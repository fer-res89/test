"""The `forge blast-species` subcommand: append best BLAST hit per
species (taxid) to each FASTA.
"""
from __future__ import annotations

import os
import re
import sys
import argparse
from pathlib import Path
from typing import List, Set, Tuple

from rich.console import Console
from rich.text import Text

from ._progress import TickingStatus

def _is_blast_db_prefix(p: Path) -> bool:
    return any(
        (p.with_suffix(sfx)).exists()
        for sfx in (".nin", ".nhr", ".nsq", ".pin", ".phr", ".psq", ".nal", ".pal")
    )

def _parse_alias_members(alias_prefix: str) -> List[str]:
    """
    Parse BLAST alias DB members from .nal/.pal.
    Returns member DB prefixes as strings.
    """
    alias_files = [Path(alias_prefix + ".nal"), Path(alias_prefix + ".pal")]
    alias_file = next((p for p in alias_files if p.exists()), None)
    if alias_file is None:
        return []

    members: List[str] = []
    in_dblist = False

    with alias_file.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("DBLIST "):
                in_dblist = True
                line = line.split(None, 1)[1]
                members.extend(line.split())
                continue

            if in_dblist:
                # stop when another KEY value section starts
                if "=" in line and not line.startswith("/") and not line.startswith("."):
                    break
                members.extend(line.split())

    return members

def _discover_db_prefixes_in_dir(db_dir: Path) -> List[str]:
    """
    Find BLAST DB prefixes in a directory.
    Excludes alias files and split-volume suffixes like .00.nin.
    """
    prefixes: Set[str] = set()

    for p in db_dir.iterdir():
        if not p.is_file():
            continue

        name = p.name

        # ignore alias files themselves
        if name.endswith(".nal") or name.endswith(".pal"):
            continue

        # nucleotide/protein main index files
        m = re.match(r"^(.*?)(?:\.(?:nhr|nin|nsq|phr|pin|psq))$", name)
        if not m:
            continue

        prefix = m.group(1)

        # skip split volume suffixes like foo.00, foo.01
        if re.search(r"\.\d{2}$", prefix):
            continue

        prefixes.add(str(db_dir / prefix))

    return sorted(prefixes)

def resolve_blast_targets(db_arg: str) -> Tuple[str, List[str]]:
    """
    Returns (db_mode, db_prefixes)
      db_mode in {"single", "alias", "directory"}
    """
    p = Path(db_arg)

    if p.is_dir():
        dbs = _discover_db_prefixes_in_dir(p)
        if not dbs:
            raise FileNotFoundError(f"No BLAST DBs found in directory: {p}")
        return "directory", dbs

    alias_members = _parse_alias_members(db_arg)
    if alias_members:
        return "alias", alias_members

    if _is_blast_db_prefix(p):
        return "single", [db_arg]

    raise FileNotFoundError(f"Could not interpret -db as a BLAST DB, alias, or directory: {db_arg}")

def run_blast_species(args: argparse.Namespace) -> None:
    """
    Append best-hit sequence(s) to each input FASTA.

    Accepted -db inputs:
      - single BLAST DB prefix
      - BLAST alias DB prefix (.nal/.pal)
      - directory containing many DB prefixes

    Behavior:
      - single DB   -> search that DB once per FASTA, using args.t threads
      - multi DB    -> create one (FASTA, DB) job per pair and run with a global executor

    Modes:
      --append full  : append full subject record (blastdbcmd -entry)
      --append range : append hit region
          - blastn vs nucl DB  : append subject DNA slice (blastdbcmd -range ...)
          - tblastn vs nucl DB : append translated AA HSP sequence (sseq from BLAST outfmt)
          - blastp/blastx vs prot DB : treated as full
      --append smart : per-hit decision
    """
    import re
    import time
    import tempfile
    import threading
    import subprocess
    import queue
    from pathlib import Path
    from datetime import datetime
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from typing import Dict, List, Set

    run_start = time.time()

    def now_stamp() -> str:
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def fmt_elapsed(seconds: float) -> str:
        total = int(round(seconds))
        h, rem = divmod(total, 3600)
        m, s = divmod(rem, 60)
        return f"{h:02d}:{m:02d}:{s:02d}"

    def estimate_eta(done: int, total: int) -> tuple[str, float]:
        if done <= 0:
            return "estimating...", 0.0
        elapsed = max(time.time() - run_start, 1e-9)
        rate = done / elapsed
        remaining = max(total - done, 0)
        eta_seconds = remaining / rate if rate > 0 else 0.0
        return fmt_elapsed(eta_seconds), rate

    console = Console(stderr=True)
    console_lock = threading.Lock()

    status: TickingStatus | None = None

    def log_console(msg: str) -> None:
        with console_lock:
            console.print(msg)

    def render_blast_status(
        *,
        started_files: int,
        total_files: int,
        total_jobs_done: int,
        total_jobs: int,
        rate: float,
        eta_str: str,
        spinner: str,
        files_in_progress: int,
        queued_for_append_count: int,
        appended_count: int,
    ) -> Text:
        parts = [
            f"{spinner} Blasting",
            f"{started_files}/{total_files} files active/started",
            f"{files_in_progress} files in progress",
            f"{total_jobs_done}/{total_jobs} jobs done",
            f"{rate:.2f} jobs/s",
            f"ETA {eta_str}",
            f"{queued_for_append_count} queued for append",
            f"{appended_count}/{total_files} appended",
        ]
        return Text(" | ".join(parts))

    def file_already_finished(out_fasta: str, log_file: str) -> bool:
        if not (os.path.exists(out_fasta) and os.path.exists(log_file)):
            return False
        try:
            with open(log_file, "r") as fh:
                log_text = fh.read()
        except OSError:
            return False
        return f"sequence(s) to {out_fasta}" in log_text

    def append_global_log_line(global_log_file: str, line: str) -> None:
        with open(global_log_file, "a") as gh:
            gh.write(line + "\n")

    fasta_files = sorted(
        [
            f for f in os.listdir(".")
            if os.path.isfile(f)
            and f.endswith(args.e)
            and not f.endswith(f".extended{args.e}")
            and ".extended" not in os.path.splitext(f)[0]
            and not f.startswith(".")
        ]
    )
    if not fasta_files:
        sys.exit(f"No files with extension '{args.e}' found in the current directory.")

    append_mode = getattr(args, "append", "smart")
    if append_mode not in ("full", "range", "smart"):
        sys.exit(f"Invalid --append mode: {append_mode}")

    max_full_bp = int(getattr(args, "max_full_bp", 200000))
    multi_db_threads = int(getattr(args, "db_threads", 1))
    if multi_db_threads < 1:
        sys.exit("--db-threads must be >= 1")

    transcript_like = re.compile(
        r"(?:^|[^A-Z0-9])(?:NM_|NR_|XM_|XR_|ENST|ENS[A-Z]{0,5}T|TRINITY_|MSTRG|mrna|transcript)(?:[^A-Z0-9]|$)",
        re.IGNORECASE,
    )

    contig_like = re.compile(
        r"(?:^|[^A-Z0-9])(?:CM\d+\.\d+|NC_\d+\.\d+|NW_\d+\.\d+|NZ_[A-Z0-9]+\.\d+|[A-Z]{4,10}\d+\.\d+|scaffold|contig|chromosome)(?:[^A-Z0-9]|$)",
        re.IGNORECASE,
    )

    def db_label(db: str) -> str:
        return Path(db).name

    def normalize_id(s: str) -> str:
        if "|" in s:
            parts = s.split("|")
            if len(parts) >= 2 and parts[1]:
                return parts[1]
        return s

    # taxid -> resolved "Genus_species" name (or a placeholder if it can't be
    # resolved). Keyed by taxid rather than by db or entry_id: a taxid always
    # means the same species regardless of which DB/entry it came from, and
    # keying this way means a single-species target DB (the common case for
    # this tool) needs at most one blastdbcmd call total, no matter how many
    # hits are appended against it.
    species_name_cache: Dict[str, str] = {}
    species_name_cache_lock = threading.Lock()

    def resolve_species_name(hit_db: str, entry_id: str, taxid: str) -> str:
        with species_name_cache_lock:
            cached = species_name_cache.get(taxid)
        if cached is not None:
            return cached

        name = None
        try:
            proc = subprocess.run(
                ["blastdbcmd", "-db", hit_db, "-entry", entry_id, "-outfmt", "%S"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
            )
            out_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
            if out_lines:
                name = out_lines[0].strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            name = None

        if not name or name.upper() in ("N/A", "NA"):
            name = f"taxid{taxid}"
            log_console(
                f"Warning: could not resolve a scientific name for taxid {taxid} "
                f"(db={db_label(hit_db)}, entry={entry_id}); using placeholder "
                f"'{name}'. Is this BLAST DB missing taxonomy info "
                f"(taxdb.bti/taxdb.btd, or built without -taxid_map)?"
            )
        else:
            name = re.sub(r"\s+", "_", name)
            name = re.sub(r"[^A-Za-z0-9_.-]", "", name)

        with species_name_cache_lock:
            species_name_cache[taxid] = name
        return name

    def clean_header(entry_id: str, species_name: str) -> str:
        return f">{entry_id}_{species_name}"

    def get_db_kind(db: str) -> str:
        try:
            proc = subprocess.run(
                ["blastdbcmd", "-db", db, "-info"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        except FileNotFoundError:
            sys.exit("Error: 'blastdbcmd' not found on PATH. Install BLAST+.")
        except subprocess.CalledProcessError as e:
            msg = (e.stdout or "").lower()
            if "protein" in msg and "nucleotide" not in msg:
                return "protein"
            if "nucleotide" in msg or "nucl" in msg:
                return "nucleotide"
            return "unknown"

        t = proc.stdout or ""
        for line in t.splitlines():
            if line.lower().strip().startswith("type:"):
                val = line.split(":", 1)[1].strip().lower()
                if "protein" in val:
                    return "protein"
                if "nucleotide" in val or "nucl" in val:
                    return "nucleotide"

        tlow = t.lower()
        if "protein" in tlow and "nucleotide" not in tlow:
            return "protein"
        if "nucleotide" in tlow or "nucl" in tlow:
            return "nucleotide"
        return "unknown"

    def fetch_full_subject(db: str, entry_id: str) -> str | None:
        try:
            proc = subprocess.run(
                ["blastdbcmd", "-db", db, "-entry", entry_id],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError:
            log_console(
                f"Warning: entry '{entry_id}' not found in DB '{db}'. Was DB built with -parse_seqids?"
            )
            return None
        return proc.stdout

    def full_bp(fasta_text: str) -> int:
        return sum(len(x.strip()) for x in fasta_text.splitlines() if x and not x.startswith(">"))

    def parse_best_hit_from_tsv(per_db_res: str, prog_name: str):
        best_for_db = None

        with open(per_db_res) as br:
            for line in br:
                parts = line.rstrip("\n").split("\t") if "\t" in line else line.strip().split()
                if len(parts) < 3:
                    continue

                raw_id = parts[0]
                sseqid = normalize_id(raw_id)
                taxids = parts[1]
                score_str = parts[2]
                taxid = (taxids.split(";", 1)[0] if taxids else "NA")
                if taxid in ("0", "NA", ""):
                    taxid = sseqid

                try:
                    score = float(score_str)
                except ValueError:
                    continue

                if prog_name == "tblastn":
                    if len(parts) < 7:
                        continue
                    try:
                        sstart = int(parts[3])
                        send = int(parts[4])
                        sframe = int(parts[5])
                    except ValueError:
                        continue
                    sseq = parts[6]
                    payload = (sseqid, score, sstart, send, sframe, sseq)

                elif prog_name == "blastn":
                    if len(parts) < 6:
                        continue
                    try:
                        sstart = int(parts[3])
                        send = int(parts[4])
                    except ValueError:
                        continue
                    sstrand = parts[5]
                    payload = (sseqid, score, sstart, send, sstrand)

                else:
                    payload = (sseqid, score)

                if best_for_db is None or score > best_for_db[1][1]:
                    best_for_db = (taxid, payload)

        return best_for_db

    def search_one_fasta_db(
        fasta: str,
        hit_db: str,
        blast_prog: str,
        blast_threads: int,
        outfmt: str,
        extra_args: list[str],
    ):
        worker_logs: List[str] = []

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv", delete=False) as tmp:
            per_db_res = tmp.name

        cmd = [
            blast_prog,
            "-query", fasta,
            "-db", hit_db,
            "-num_threads", str(blast_threads),
            "-max_target_seqs", "50",
            "-max_hsps", "1",
            "-outfmt", outfmt,
            "-out", per_db_res,
        ]
        if extra_args:
            cmd.extend(extra_args)

        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            try:
                os.remove(per_db_res)
            except OSError:
                pass
            raise RuntimeError(f"Error: '{blast_prog}' not found on PATH. Install BLAST+ and try again.")
        except subprocess.CalledProcessError as e:
            worker_logs.append(f"BLAST failed for db={db_label(hit_db)} query={fasta}: {e}")
            try:
                os.remove(per_db_res)
            except OSError:
                pass
            return fasta, hit_db, None, worker_logs

        try:
            best_for_db = parse_best_hit_from_tsv(per_db_res, os.path.basename(blast_prog).lower())
        finally:
            try:
                os.remove(per_db_res)
            except OSError:
                pass

        return fasta, hit_db, best_for_db, worker_logs

    db_mode, db_list = resolve_blast_targets(args.db)
    prog = os.path.basename(args.tl).lower()

    is_blastn = prog == "blastn"
    is_tblastn = prog == "tblastn"
    is_blastp = prog == "blastp"
    is_blastx = prog == "blastx"

    if not (is_blastn or is_tblastn or is_blastp or is_blastx):
        log_console(
            f"Warning: unrecognized BLAST program '{args.tl}'. Supported: blastn, tblastn, blastp, blastx."
        )

    if is_tblastn:
        outfmt = "6 sseqid staxids bitscore sstart send sframe sseq"
    elif is_blastn:
        outfmt = "6 sseqid staxids bitscore sstart send sstrand"
    else:
        outfmt = "6 sseqid staxids bitscore"

    db_kind_by_db: Dict[str, str] = {}
    for db in db_list:
        db_kind_by_db[db] = get_db_kind(db)

    log_console(f"[DEBUG] db_mode={db_mode} n_dbs={len(db_list)} prog={prog} outfmt='{outfmt}'")

    fasta_info: Dict[str, Dict[str, str]] = {}
    fasta_logs: Dict[str, List[str]] = {}
    best_hits_by_fasta: Dict[str, Dict[str, tuple]] = {}

    global_log_file = "blast_species.global.log"
    if not (getattr(args, "resume", False) and os.path.exists(global_log_file)):
        with open(global_log_file, "w") as gh:
            gh.write(f"[START] {now_stamp()}\n")
    else:
        with open(global_log_file, "a") as gh:
            gh.write(f"[RESUME START] {now_stamp()}\n")

    files_to_run: List[str] = []
    skipped_files: List[str] = []

    for fasta in fasta_files:
        base = fasta[:-len(args.e)] if args.e else os.path.splitext(fasta)[0]
        out_fasta = f"{base}.extended{args.e}"
        log_file = f"{base}.blast.log"

        fasta_info[fasta] = {"out_fasta": out_fasta, "log_file": log_file}
        best_hits_by_fasta[fasta] = {}

        if getattr(args, "resume", False) and file_already_finished(out_fasta, log_file):
            skipped_files.append(fasta)
            append_global_log_line(global_log_file, f'Skipping finished "{fasta}" → "{out_fasta}"')
            continue

        fasta_logs[fasta] = [
            f'Processing "{fasta}" → "{out_fasta}"',
            f"db_mode={db_mode} n_dbs={len(db_list)} prog={prog} outfmt={outfmt}",
        ]

        append_global_log_line(global_log_file, f'Processing "{fasta}" → "{out_fasta}"')

        with open(out_fasta, "w") as out_h, open(fasta) as in_h:
            out_h.write(in_h.read())

        files_to_run.append(fasta)

    state_lock = threading.Lock()
    append_state_lock = threading.Lock()

    active_db_counts: Dict[str, int] = {fasta: 0 for fasta in files_to_run}
    completed_jobs_by_fasta: Dict[str, int] = {fasta: 0 for fasta in files_to_run}
    total_files = len(files_to_run)
    total_jobs = len(files_to_run) * len(db_list)
    total_jobs_done = 0

    queued_for_append: Set[str] = set()
    appended_files: Set[str] = set()
    append_errors: List[tuple[str, str]] = []

    append_queue: queue.Queue[str | None] = queue.Queue()

    if getattr(args, "resume", False):
        log_console(
            f'{len(files_to_run)} files being processed, {len(skipped_files)} already finished and skipped'
        )
    else:
        log_console(f"{len(files_to_run)} files being processed")

    if not files_to_run:
        elapsed = time.time() - run_start
        if getattr(args, "resume", False):
            msg = f'Nothing to do: all matching files are already finished. Total time: {fmt_elapsed(elapsed)}'
        else:
            msg = f'No files with extension "{args.e}" need processing. Total time: {fmt_elapsed(elapsed)}'
        log_console(msg)
        append_global_log_line(global_log_file, f"[END] {now_stamp()} | elapsed={fmt_elapsed(elapsed)} | nothing to do")
        return

    def write_one_fasta_log(fasta: str) -> None:
        log_file = fasta_info[fasta]["log_file"]
        with append_state_lock:
            lines = list(fasta_logs[fasta])
        with open(log_file, "w") as log_h:
            for line in lines:
                log_h.write(line + "\n")

    def compute_render(spin: str) -> Text:
        """Pure(ish) renderer: reads shared state, returns a Text.
        Called by TickingStatus on its own fixed-interval schedule -- this
        is what keeps the spinner moving even when no BLAST job has
        finished in a while."""
        nonlocal total_jobs_done
        eta_str, rate = estimate_eta(total_jobs_done, total_jobs)

        with state_lock:
            started_files = sum(
                1
                for f in files_to_run
                if completed_jobs_by_fasta.get(f, 0) > 0 or active_db_counts.get(f, 0) > 0
            )

            files_in_progress = sum(
                1
                for f in files_to_run
                if 0 < completed_jobs_by_fasta.get(f, 0) < len(db_list)
            )

        with append_state_lock:
            appended_count = len(appended_files)
            queued_count = max(0, len(queued_for_append) - appended_count)

        return render_blast_status(
            started_files=started_files,
            total_files=total_files,
            total_jobs_done=total_jobs_done,
            total_jobs=total_jobs,
            rate=rate,
            eta_str=eta_str,
            spinner=spin,
            files_in_progress=files_in_progress,
            queued_for_append_count=queued_count,
            appended_count=appended_count,
        )

    def update_blast_progress(force: bool = False) -> None:
        """Ask for an immediate redraw reflecting the latest state (e.g.
        right after a job starts/finishes), on top of the TickingStatus
        background ticker that redraws at a constant rate regardless.
        `force` is accepted for call-site compatibility but no longer
        gates anything -- the ticker thread is what guarantees the display
        never goes stale, so there's nothing left to throttle here."""
        if status is not None:
            status.refresh()

    def append_one_finished_fasta(fasta: str) -> None:
        out_fasta = fasta_info[fasta]["out_fasta"]
        best_hits = best_hits_by_fasta[fasta]

        if not best_hits:
            with append_state_lock:
                fasta_logs[fasta].append(f"No BLAST hits found in {fasta}")
                fasta_logs[fasta].append(f"Appended 0 sequence(s) to {out_fasta}")
            write_one_fasta_log(fasta)
            log_console(f'Finished "{fasta}": appended 0 sequence(s)')
            return

        appended = 0
        local_log_lines: List[str] = []

        for hit_db in db_list:
            if hit_db not in best_hits:
                continue

            _db_key, taxid, payload = best_hits[hit_db]
            entry_id = payload[0]
            score = payload[1]
            db_kind = db_kind_by_db.get(hit_db, "unknown")

            mode = append_mode

            if mode == "smart":
                if is_blastp or is_blastx:
                    mode = "full"
                elif is_tblastn:
                    if transcript_like.search(entry_id):
                        mode = "full"
                    elif contig_like.search(entry_id):
                        mode = "range"
                    else:
                        mode = "range"
                elif is_blastn:
                    if transcript_like.search(entry_id):
                        mode = "full"
                    elif contig_like.search(entry_id):
                        mode = "range"
                    else:
                        mode = "full"
                else:
                    mode = "full"

                local_log_lines.append(
                    f"smart decision: db={db_label(hit_db)} entry_id={entry_id!r} "
                    f"db_kind={db_kind!r} contig_match={bool(contig_like.search(entry_id))} "
                    f"transcript_match={bool(transcript_like.search(entry_id))} mode={mode!r}"
                )

            if mode == "range":
                if is_tblastn:
                    _, _, sstart, send, sframe, sseq = payload
                    aa = sseq.replace("-", "").replace("*", "").strip()
                    if not aa:
                        log_console(f"Warning: empty HSP AA for {entry_id} taxid {taxid}")
                        continue

                    start = min(sstart, send)
                    end = max(sstart, send)
                    species_name = resolve_species_name(hit_db, entry_id, taxid)
                    header = clean_header(entry_id, species_name)

                    with open(out_fasta, "a") as out_h:
                        out_h.write(header + "\n")
                        for i in range(0, len(aa), 60):
                            out_h.write(aa[i:i+60] + "\n")

                    appended += 1
                    local_log_lines.append(
                        f"Appended RANGE-AA (tblastn HSP) for taxid {taxid} "
                        f"({species_name}): {entry_id} HSP {start}-{end} "
                        f"frame {sframe} from {db_label(hit_db)} (bitscore {score})"
                    )
                    continue

                if is_blastn:
                    _, _, sstart, send, sstrand = payload
                    start = min(sstart, send)
                    end = max(sstart, send)
                    strand = "plus" if str(sstrand).lower() in ("plus", "+") else "minus"

                    try:
                        proc = subprocess.run(
                            ["blastdbcmd", "-db", hit_db, "-entry", entry_id,
                             "-range", f"{start}-{end}", "-strand", strand],
                            check=True,
                            stdout=subprocess.PIPE,
                            text=True,
                        )
                    except subprocess.CalledProcessError:
                        log_console(f"Warning: range fetch failed for {entry_id}:{start}-{end} from {hit_db}")
                        continue

                    lines = proc.stdout.splitlines()
                    species_name = resolve_species_name(hit_db, entry_id, taxid)
                    if lines and lines[0].startswith(">"):
                        lines[0] = clean_header(entry_id, species_name)

                    with open(out_fasta, "a") as out_h:
                        out_h.write("\n".join(lines) + ("\n" if lines else ""))

                    appended += 1
                    local_log_lines.append(
                        f"Appended RANGE-DNA (blastn HSP) for taxid {taxid} "
                        f"({species_name}): {entry_id} HSP {start}-{end} "
                        f"strand {strand} from {db_label(hit_db)} (bitscore {score})"
                    )
                    continue

                mode = "full"

            if mode == "full":
                if append_mode == "smart" and is_blastn:
                    contig_hit = bool(contig_like.search(entry_id))
                    transcript_hit = bool(transcript_like.search(entry_id))
                    if contig_hit and not transcript_hit:
                        log_console(f"Warning: refusing FULL fetch for contig-like ID in smart mode: {entry_id}")
                        continue

                full = fetch_full_subject(hit_db, entry_id)
                if not full:
                    continue

                if append_mode == "smart" and is_blastn:
                    bp = full_bp(full)
                    if bp > max_full_bp:
                        if len(payload) >= 5:
                            _, _, sstart, send, sstrand = payload
                            start = min(sstart, send)
                            end = max(sstart, send)
                            strand = "plus" if str(sstrand).lower() in ("plus", "+") else "minus"

                            try:
                                proc = subprocess.run(
                                    ["blastdbcmd", "-db", hit_db, "-entry", entry_id,
                                     "-range", f"{start}-{end}", "-strand", strand],
                                    check=True,
                                    stdout=subprocess.PIPE,
                                    text=True,
                                )
                            except subprocess.CalledProcessError:
                                log_console(
                                    f"Warning: full too big ({bp} bp) and range fallback failed for {entry_id} from {hit_db}"
                                )
                                continue

                            lines = proc.stdout.splitlines()
                            species_name = resolve_species_name(hit_db, entry_id, taxid)
                            if lines and lines[0].startswith(">"):
                                lines[0] = clean_header(entry_id, species_name)

                            with open(out_fasta, "a") as out_h:
                                out_h.write("\n".join(lines) + ("\n" if lines else ""))

                            appended += 1
                            local_log_lines.append(
                                f"Appended RANGE-DNA fallback for taxid {taxid} "
                                f"({species_name}): {entry_id} HSP {start}-{end} "
                                f"strand {strand} from {db_label(hit_db)} (bitscore {score})"
                            )
                            continue
                        else:
                            log_console(
                                f"Warning: full too big ({bp} bp) for {entry_id} but no coords available; "
                                f"consider using blastn outfmt with coords."
                            )

                with open(out_fasta, "a") as out_h:
                    out_h.write(full if full.endswith("\n") else full + "\n")

                appended += 1
                local_log_lines.append(
                    f"Appended FULL for taxid {taxid}: {entry_id} "
                    f"from {db_label(hit_db)} (bitscore {score})"
                )
                continue

        local_log_lines.append(f"Appended {appended} sequence(s) to {out_fasta}")

        with append_state_lock:
            fasta_logs[fasta].extend(local_log_lines)

        write_one_fasta_log(fasta)
        log_console(f'Finished "{fasta}": appended {appended} sequence(s)')

    def queue_ready_fasta_for_append(fasta: str) -> None:
        if completed_jobs_by_fasta.get(fasta, 0) < len(db_list):
            return
        with append_state_lock:
            if fasta in queued_for_append:
                return
            queued_for_append.add(fasta)
        append_queue.put(fasta)
        update_blast_progress(force=True)

    def append_worker() -> None:
        while True:
            fasta = append_queue.get()
            try:
                if fasta is None:
                    return

                append_one_finished_fasta(fasta)

                with append_state_lock:
                    appended_files.add(fasta)

                update_blast_progress(force=True)

            except Exception as e:
                with append_state_lock:
                    append_errors.append((str(fasta), str(e)))
                    if isinstance(fasta, str) and fasta in fasta_logs:
                        fasta_logs[fasta].append(f"Append worker error: {e}")
                if isinstance(fasta, str) and fasta in fasta_info:
                    write_one_fasta_log(fasta)
                log_console(f'Error while appending "{fasta}": {e}')
                update_blast_progress(force=True)
            finally:
                append_queue.task_done()

    append_thread = threading.Thread(
        target=append_worker,
        name="blast-species-append-worker",
        daemon=True,
    )
    append_thread.start()

    try:
        with TickingStatus(
            console,
            compute_render,
            interval=0.15,
            refresh_per_second=10,
            transient=True,
        ) as ticking_status:
            status = ticking_status

            if len(db_list) == 1:
                single_db = db_list[0]
                blast_threads = args.t

                for fasta in files_to_run:
                    with state_lock:
                        active_db_counts[fasta] += 1
                    update_blast_progress(force=True)

                    fasta_name, hit_db, best_for_db, worker_logs = search_one_fasta_db(
                        fasta=fasta,
                        hit_db=single_db,
                        blast_prog=args.tl,
                        blast_threads=blast_threads,
                        outfmt=outfmt,
                        extra_args=args.extra,
                    )

                    with state_lock:
                        active_db_counts[fasta] -= 1
                        completed_jobs_by_fasta[fasta] += 1
                        total_jobs_done += 1

                    with append_state_lock:
                        fasta_logs[fasta_name].extend(worker_logs)

                    if best_for_db is not None:
                        taxid, payload = best_for_db
                        best_hits_by_fasta[fasta_name][hit_db] = (hit_db, taxid, payload)

                    update_blast_progress(force=True)
                    queue_ready_fasta_for_append(fasta)

            else:
                all_jobs = [(fasta, db) for fasta in files_to_run for db in db_list]

                if args.jobs is None:
                    cpu_n = os.cpu_count() or 1
                    max_workers = max(1, cpu_n // multi_db_threads)
                    max_workers = min(max_workers, len(all_jobs))
                else:
                    max_workers = max(1, int(args.jobs))

                log_console(
                    f"[DEBUG] multi-db mode: scheduling {len(all_jobs)} jobs "
                    f"with max_workers={max_workers}, db_threads={multi_db_threads}"
                )

                with ThreadPoolExecutor(max_workers=max_workers) as ex:
                    futures = {}

                    for fasta, db in all_jobs:
                        with state_lock:
                            active_db_counts[fasta] += 1

                        fut = ex.submit(
                            search_one_fasta_db,
                            fasta,
                            db,
                            args.tl,
                            multi_db_threads,
                            outfmt,
                            args.extra,
                        )
                        futures[fut] = (fasta, db)

                    update_blast_progress(force=True)

                    for fut in as_completed(futures):
                        fasta, db = futures[fut]
                        try:
                            fasta_name, hit_db, best_for_db, worker_logs = fut.result()
                        except Exception as e:
                            with append_state_lock:
                                fasta_logs[fasta].append(
                                    f"Worker failed for db={db_label(db)} query={fasta}: {e}"
                                )
                            with state_lock:
                                active_db_counts[fasta] -= 1
                                completed_jobs_by_fasta[fasta] += 1
                                total_jobs_done += 1
                            update_blast_progress(force=True)
                            queue_ready_fasta_for_append(fasta)
                            continue

                        with state_lock:
                            active_db_counts[fasta_name] -= 1
                            completed_jobs_by_fasta[fasta_name] += 1
                            total_jobs_done += 1

                        with append_state_lock:
                            fasta_logs[fasta_name].extend(worker_logs)

                        if best_for_db is not None:
                            taxid, payload = best_for_db
                            best_hits_by_fasta[fasta_name][hit_db] = (hit_db, taxid, payload)

                        update_blast_progress(force=True)
                        queue_ready_fasta_for_append(fasta_name)

            for fasta in files_to_run:
                queue_ready_fasta_for_append(fasta)

            while True:
                with append_state_lock:
                    done_appending = len(appended_files)
                    queued_count = len(queued_for_append)

                if done_appending >= queued_count:
                    break

                update_blast_progress(force=True)
                time.sleep(0.2)

            update_blast_progress(force=True)
            status = None

    finally:
        append_queue.put(None)
        append_thread.join()

    if append_errors:
        log_console(f"Warning: {len(append_errors)} file(s) had append-thread errors.")

    # WRITE LOG FILES (final sync pass)
    for fasta in files_to_run:
        write_one_fasta_log(fasta)

    elapsed = time.time() - run_start
    elapsed_str = fmt_elapsed(elapsed)
    log_console(f"All done in {elapsed_str}")
    append_global_log_line(global_log_file, f"[END] {now_stamp()} | elapsed={elapsed_str}")


# -------------------- GPFF → taxid map --------------------


def register(sub, formatter_class) -> None:
    """Register the `blast-species` subcommand."""
    sp_blast = sub.add_parser(
        "blast-species",
        help=("Append best BLAST hit per species (taxid) to each FASTA.\n"
              "(run 'forge blast-species -h' for options)" 
        ),
        epilog=(
            "PREREQUISITES\n"
            "-------------\n"
            "- BLAST+ (blastn, blastp, tblastn, blastdbcmd)\n"
            "- BLAST DB built with 'parse_seqids' so that sequence IDs can be retrieved.\n"
            "  e.g.:\n"
            "    * makeblastdb -in target_sequences.fasta -dbtype prot -parse_seqids -out target_db\n"
            "    * cat file1.fasta file2.fasta file3.fasta > all_seqs.fasta;\n"
            "      makeblastdb -in all_seqs.fasta -dbtype nucl -parse_seqids -out multi_species_db\n"
            "        - if database includes multiple species, database needs to be created using -taxid_map option\n"
            "        - if taxid map text file is needed, use get-taxid function (run 'forge get-taxid -h' for help)\n"
        ),
        formatter_class=formatter_class,
    )
    sp_blast.add_argument("-e", metavar="EXTENSION", default=".fasta", help="FASTA extension to scan (default: .fasta)")
    sp_blast.add_argument(
        "-db",
        metavar="DATA_BASE",
        required=True,
        help="BLAST database prefix, alias prefix, or directory containing many DBs"
    )
    sp_blast.add_argument("-tl", metavar="PROGRAM", default="blastn", help="BLAST+ program (blastn, blastp, tblastn)")
    sp_blast.add_argument("-t", metavar="THREADS", type=int, default=1, help="Threads")
    sp_blast.add_argument(
        "--jobs",
        type=int,
        default=None,
        help=("Concurrent BLAST jobs in multi-DB mode. "
              "Default: auto = min(total jobs, os.cpu_count() // --db-threads). "
              "Ignored for single-DB mode.")
    )
    sp_blast.add_argument(
        "--db-threads",
        type=int,
        default=1,
        help=("Threads per BLAST process in multi-DB mode. "
              "Single-DB mode still uses -t.")
    )
    sp_blast.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="Pass extra BLAST flags after --")
    sp_blast.add_argument(
        "--append",
        choices=["full", "range", "smart"],
        default="smart",
        help="Append full subject entries, hit ranges, or smart auto-detect (recommended for mixed DBs)."
    )
    sp_blast.add_argument(
        "--max-full-bp",
        type=int,
        default=200000,
        help="Smart mode: if full fetched nucleotide entry exceeds this many bp, fall back to range (blastn only)."
    )
    sp_blast.add_argument(
        "--resume",
        action="store_true",
        help=("Skip FASTA files whose output and log already indicate completion.")
    )
    sp_blast.set_defaults(func=run_blast_species)
