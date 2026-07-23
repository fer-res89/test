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
    import collections
    import shutil
    from pathlib import Path
    from datetime import datetime
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from typing import Dict, List, Set

    run_start = time.time()

    # Rolling window for ETA rate estimation. A pure cumulative average
    # (total_jobs_done / elapsed_since_start) gets permanently skewed by
    # whatever happened early in a long run (thread-pool warm-up, cold disk
    # cache, an unusually slow first DB) -- and for a job with hundreds of
    # thousands of units, "early" can be the first many minutes, during
    # which the cumulative average barely moves even once steady-state
    # throughput is much faster. A trailing window recovers from that within
    # RATE_WINDOW_SECONDS instead of dragging the skew for the whole run.
    RATE_WINDOW_SECONDS = 45.0
    rate_samples: "collections.deque[tuple[float, int]]" = collections.deque()
    rate_samples_lock = threading.Lock()

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

        now = time.time()
        with rate_samples_lock:
            rate_samples.append((now, done))
            cutoff = now - RATE_WINDOW_SECONDS
            while len(rate_samples) > 1 and rate_samples[0][0] < cutoff:
                rate_samples.popleft()
            oldest_t, oldest_done = rate_samples[0]

        span = now - oldest_t
        done_delta = done - oldest_done
        windowed_rate = (done_delta / span) if (span >= 3.0 and done_delta > 0) else None

        elapsed = max(now - run_start, 1e-9)
        cumulative_rate = done / elapsed

        # Use the windowed rate once we actually have a few seconds of
        # recent history to compute it from; otherwise (very start of the
        # run) fall back to the cumulative average, which is all we have.
        rate = windowed_rate if windowed_rate is not None else cumulative_rate

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
        total_jobs_done: int,
        total_jobs: int,
        rate: float,
        eta_str: str,
        elapsed_str: str,
        spinner: str,
        appended_count: "int | None" = None,
        total_files: "int | None" = None,
    ) -> Text:
        parts = [
            f"{spinner} Blasting",
            f"{total_jobs_done}/{total_jobs} jobs done",
            f"{rate:.2f} jobs/s",
            f"ETA {eta_str}",
            f"Elapsed {elapsed_str}",
        ]
        # Per-file mode: files finish (and get appended) incrementally
        # throughout the run, so "N/total files appended" is a meaningful
        # ratio. Batch mode appends individual HITS incrementally instead
        # (one per (file, db) pair, as each db finishes) with no fixed total
        # to compare against, so it gets a running count with no
        # denominator rather than a ratio that doesn't mean anything.
        if appended_count is not None and total_files is not None:
            parts.append(f"{appended_count}/{total_files} appended")
        elif appended_count is not None:
            parts.append(f"{appended_count} hits appended so far")
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

        # blastdbcmd only finds taxdb.bti/taxdb.btd via the BLASTDB env var
        # (or cwd) -- it does NOT automatically check the directory the DB
        # itself lives in. Since we already know exactly where this DB is,
        # explicitly put its directory on BLASTDB for this call, so taxonomy
        # lookups work regardless of the caller's shell environment or
        # working directory (as long as the taxdb files are colocated with
        # the DB, which is the standard/recommended setup anyway).
        env = os.environ.copy()
        db_dir = os.path.dirname(os.path.abspath(hit_db)) or "."
        existing_blastdb = env.get("BLASTDB", "")
        env["BLASTDB"] = db_dir + (os.pathsep + existing_blastdb if existing_blastdb else "")

        name = None
        try:
            proc = subprocess.run(
                ["blastdbcmd", "-db", hit_db, "-entry", entry_id, "-outfmt", "%S"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                env=env,
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
                f"'{name}'. Checked BLASTDB={env['BLASTDB']!r} for "
                f"taxdb.bti/taxdb.btd -- confirm those files are actually "
                f"there and not truncated (e.g. 'tar -tzf taxdb.tar.gz')."
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

    def parse_tsv_row(parts: List[str], prog_name: str):
        """Parse one BLAST tabular-output row (already split into fields,
        WITHOUT any qseqid column) into (taxid, score, payload), or None if
        the row is short/unparseable. Shared by both per-file mode
        (parse_best_hit_from_tsv) and batch-query mode
        (parse_best_hits_grouped_by_qseqid), so the two modes can never
        silently diverge in how a hit is interpreted."""
        if len(parts) < 3:
            return None

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
            return None

        if prog_name == "tblastn":
            if len(parts) < 7:
                return None
            try:
                sstart = int(parts[3])
                send = int(parts[4])
                sframe = int(parts[5])
            except ValueError:
                return None
            sseq = parts[6]
            payload = (sseqid, score, sstart, send, sframe, sseq)

        elif prog_name == "blastn":
            if len(parts) < 6:
                return None
            try:
                sstart = int(parts[3])
                send = int(parts[4])
            except ValueError:
                return None
            sstrand = parts[5]
            payload = (sseqid, score, sstart, send, sstrand)

        else:
            payload = (sseqid, score)

        return taxid, score, payload

    def parse_best_hit_from_tsv(per_db_res: str, prog_name: str):
        best_for_db = None

        with open(per_db_res) as br:
            for line in br:
                parts = line.rstrip("\n").split("\t") if "\t" in line else line.strip().split()
                parsed = parse_tsv_row(parts, prog_name)
                if parsed is None:
                    continue
                taxid, score, payload = parsed
                if best_for_db is None or score > best_for_db[1][1]:
                    best_for_db = (taxid, payload)

        return best_for_db

    def parse_best_hits_grouped_by_qseqid(per_db_res: str, prog_name: str) -> Dict[str, tuple]:
        """Batch-query-mode variant: the TSV has an EXTRA trailing qseqid
        column (appended, so all the existing fixed-position field indices
        parse_tsv_row relies on stay valid). Returns the best (taxid,
        payload) per qseqid, exactly like parse_best_hit_from_tsv does per
        whole file in per-file mode -- just grouped, since one combined file
        now covers every query at once."""
        best_by_qseqid: Dict[str, tuple] = {}

        with open(per_db_res) as br:
            for line in br:
                parts = line.rstrip("\n").split("\t") if "\t" in line else line.strip().split()
                if len(parts) < 2:
                    continue
                qseqid = parts[-1]
                parsed = parse_tsv_row(parts[:-1], prog_name)
                if parsed is None:
                    continue
                taxid, score, payload = parsed
                cur = best_by_qseqid.get(qseqid)
                if cur is None or score > cur[2]:
                    best_by_qseqid[qseqid] = (taxid, payload, score)

        return {q: (t, p) for q, (t, p, _score) in best_by_qseqid.items()}

    def build_combined_query_fasta(files: List[str], out_path) -> Dict[str, str]:
        """Concatenate every query FASTA in `files` into one combined
        multi-FASTA for batch-query mode, rewriting each sequence's header to
        a globally unique id (filenames are already unique within one
        directory listing, so `<stem>__<i>` cannot collide). Returns
        {qseqid: original_fasta_path} so results can be attributed back to
        the file they came from without re-parsing anything."""
        qseqid_to_fasta: Dict[str, str] = {}
        with open(out_path, "w") as out_h:
            for fasta in files:
                stem = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(fasta).stem)
                seq_idx = 0
                writing = False
                with open(fasta) as in_h:
                    for line in in_h:
                        line = line.rstrip("\n")
                        if line.startswith(">"):
                            qseqid = f"{stem}__{seq_idx}"
                            seq_idx += 1
                            qseqid_to_fasta[qseqid] = fasta
                            out_h.write(f">{qseqid}\n")
                            writing = True
                        elif writing:
                            out_h.write(line + "\n")
        return qseqid_to_fasta

    def search_one_db_batch(
        hit_db: str,
        combined_query: str,
        blast_prog: str,
        blast_threads: int,
        outfmt_batch: str,
        extra_args: list[str],
        cached_result_path: str,
        use_resume: bool,
    ):
        """Batch-query-mode search: ONE blast process for this DB, searching
        every query at once, instead of one process per (fasta, db) pair.
        This is the whole point of batch mode -- eliminating the repeated
        per-query DB-load/setup overhead that dominates when each query gets
        its own tiny process.

        Unlike per-file mode, nothing gets appended to any output file until
        EVERY db has been searched (only then do we know, for each file,
        which db had its single best hit). That means a crash partway
        through a long batch run could otherwise lose all completed search
        work with nothing to show for it. To avoid that: results are written
        to `cached_result_path` (persistent, not a deleted temp file) instead
        of being thrown away after parsing. On a subsequent --resume run, if
        that file already exists, BLAST is skipped entirely for this db and
        the existing results are just re-parsed -- so a crash only costs you
        the db(s) that hadn't finished yet, not everything that had.
        """
        worker_logs: List[str] = []
        prog_name = os.path.basename(blast_prog).lower()

        if use_resume and os.path.exists(cached_result_path):
            try:
                best_by_qseqid = parse_best_hits_grouped_by_qseqid(cached_result_path, prog_name)
                worker_logs.append(
                    f"Resumed db={db_label(hit_db)} from a previously-completed "
                    f"batch search (skipped re-running BLAST for this db)"
                )
                return hit_db, best_by_qseqid, worker_logs
            except OSError:
                pass  # cached file unreadable -- fall through and re-search

        # Write to a temp file first, then atomically rename to
        # cached_result_path ONLY after BLAST confirms success. This matters:
        # if -out pointed directly at cached_result_path, an interrupted run
        # (Ctrl+C mid-search) would leave an empty/truncated file sitting at
        # that exact path -- which os.path.exists() can't distinguish from a
        # genuinely completed result. A subsequent --resume would then treat
        # that stub as "already searched" and silently skip it forever,
        # parsing zero hits from an empty file instead of re-running BLAST.
        # The atomic rename guarantees cached_result_path only ever exists
        # once fully written -- there's no window where a partial file can
        # appear under that name.
        tmp_result_path = cached_result_path + f".tmp.{os.getpid()}.{threading.get_ident()}"
        cmd = [
            blast_prog,
            "-query", combined_query,
            "-db", hit_db,
            "-num_threads", str(blast_threads),
            "-max_target_seqs", "50",
            "-max_hsps", "1",
            "-outfmt", outfmt_batch,
            "-out", tmp_result_path,
        ]
        if extra_args:
            cmd.extend(extra_args)

        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            try:
                os.remove(tmp_result_path)
            except OSError:
                pass
            raise RuntimeError(f"Error: '{blast_prog}' not found on PATH. Install BLAST+ and try again.")
        except subprocess.CalledProcessError as e:
            worker_logs.append(f"BLAST failed for db={db_label(hit_db)} (batch query): {e}")
            try:
                os.remove(tmp_result_path)
            except OSError:
                pass
            return hit_db, None, worker_logs
        except BaseException:
            # Covers Ctrl+C (KeyboardInterrupt) and any other interruption
            # reaching this point: clean up the temp file so it can't be
            # confused with a real (never-created, since we never got here)
            # cached_result_path, then re-raise.
            try:
                os.remove(tmp_result_path)
            except OSError:
                pass
            raise

        os.replace(tmp_result_path, cached_result_path)

        best_by_qseqid = parse_best_hits_grouped_by_qseqid(cached_result_path, prog_name)
        return hit_db, best_by_qseqid, worker_logs

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
    # Running per-fasta appended-hit count. Populated incrementally in batch
    # mode (as each db's hits get appended); in per-file mode it's left empty
    # and the status line falls back to the existing appended_files-based
    # count instead, since per-file mode's append still happens in one shot
    # per fasta rather than incrementally per db.
    total_appended_by_fasta: Dict[str, int] = {}

    global_log_file = "blast_species.global.log"
    if not (getattr(args, "resume", False) and os.path.exists(global_log_file)):
        with open(global_log_file, "w") as gh:
            gh.write(f"[START] {now_stamp()}\n")
    else:
        with open(global_log_file, "a") as gh:
            gh.write(f"[RESUME START] {now_stamp()}\n")

    files_to_run: List[str] = []
    skipped_files: List[str] = []
    is_batch_query_mode_early = getattr(args, "query_mode", "per-file") == "batch"

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

        # Batch mode can have a file that's NOT fully finished yet still
        # carry genuine partial progress from a previous (crashed) run --
        # some (fasta, db) hits already durably appended, tracked via
        # [APPENDED_DB] markers in its log. Per-file mode never has this
        # state (nothing gets written until a file's full hit set is ready),
        # so this only matters in batch mode. Detect it BEFORE deciding
        # whether to reset out_fasta/the in-memory log, so we don't destroy
        # hits that were already safely on disk.
        has_partial_batch_progress = False
        if is_batch_query_mode_early and getattr(args, "resume", False):
            existing_markers: Set[str] = set()
            if os.path.exists(log_file):
                try:
                    with open(log_file) as lf:
                        for line in lf:
                            line = line.strip()
                            if line.startswith("[APPENDED_DB] "):
                                existing_markers.add(line[len("[APPENDED_DB] "):].strip())
                except OSError:
                    existing_markers = set()
            has_partial_batch_progress = bool(existing_markers) and os.path.exists(out_fasta)

        if has_partial_batch_progress:
            # Preserve on-disk content and seed the in-memory log with what's
            # already there, so write_one_fasta_log appends to history
            # instead of replacing it (which would silently drop the
            # [APPENDED_DB] markers this whole mechanism depends on).
            try:
                with open(log_file) as lf:
                    fasta_logs[fasta] = [ln.rstrip("\n") for ln in lf]
            except OSError:
                fasta_logs[fasta] = []
            fasta_logs[fasta].append(
                f'Resuming "{fasta}" -> "{out_fasta}" (preserving already-appended hits)'
            )
        else:
            fasta_logs[fasta] = [
                f'Processing "{fasta}" → "{out_fasta}"',
                f"db_mode={db_mode} n_dbs={len(db_list)} prog={prog} outfmt={outfmt}",
            ]
            with open(out_fasta, "w") as out_h, open(fasta) as in_h:
                out_h.write(in_h.read())

        append_global_log_line(global_log_file, f'Processing "{fasta}" → "{out_fasta}"')

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

    def load_already_appended_dbs(fasta: str) -> Set[str]:
        """Batch-mode resume support: each successful append writes a
        machine-parseable '[APPENDED_DB] <label>' marker line (in addition
        to its descriptive log line). On --resume, read a fasta's existing
        .blast.log (if any) and collect those markers, so hits from dbs
        already appended in a prior (possibly crashed) run are never
        appended a second time."""
        log_file = fasta_info.get(fasta, {}).get("log_file")
        if not log_file or not os.path.exists(log_file):
            return set()
        done: Set[str] = set()
        try:
            with open(log_file) as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith("[APPENDED_DB] "):
                        done.add(line[len("[APPENDED_DB] "):].strip())
        except OSError:
            return set()
        return done

    def compute_render(spin: str) -> Text:
        """Pure(ish) renderer: reads shared state, returns a Text.
        Called by TickingStatus on its own fixed-interval schedule -- this
        is what keeps the spinner moving even when no BLAST job has
        finished in a while."""
        nonlocal total_jobs_done
        eta_str, rate = estimate_eta(total_jobs_done, total_jobs)
        elapsed_str = fmt_elapsed(time.time() - run_start)

        is_batch = getattr(args, "query_mode", "per-file") == "batch"
        appended_count = None
        show_total_files = None
        if is_batch:
            with append_state_lock:
                appended_count = sum(total_appended_by_fasta.values())
            # No natural fixed denominator for "hits appended" in batch mode
            # (unlike per-file mode's N/total-files), so show a running
            # count only.
        else:
            with append_state_lock:
                appended_count = len(appended_files)
            show_total_files = total_files

        return render_blast_status(
            total_jobs_done=total_jobs_done,
            total_jobs=total_jobs,
            rate=rate,
            eta_str=eta_str,
            elapsed_str=elapsed_str,
            spinner=spin,
            appended_count=appended_count,
            total_files=show_total_files,
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

    fasta_append_locks: Dict[str, threading.Lock] = {fasta: threading.Lock() for fasta in files_to_run}

    def format_and_append_one_hit(fasta: str, hit_db: str, taxid: str, payload: tuple) -> tuple[List[str], bool]:
        """Append ONE (fasta, db) hit to fasta's out_fasta file right now --
        a true, permanent file append, not part of a larger batched rewrite.
        Returns (log_lines, appended_ok). Safe to call concurrently for
        DIFFERENT fastas; calls for the SAME fasta (which can happen in
        batch mode, where multiple DBs' results land around the same time)
        are serialized via fasta_append_locks so two hits never interleave
        their writes to the same file.

        This is the exact per-hit logic previously inlined in
        append_one_finished_fasta's per-db loop, factored out so it can be
        invoked immediately per (fasta, db) pair -- in per-file mode once a
        fasta's own jobs are all done (unchanged from before), and in batch
        mode the moment EACH db's combined search returns, without waiting
        for every other db too.
        """
        out_fasta = fasta_info[fasta]["out_fasta"]
        entry_id = payload[0]
        score = payload[1]
        db_kind = db_kind_by_db.get(hit_db, "unknown")
        local_log_lines: List[str] = []
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

        lock = fasta_append_locks.get(fasta)
        if lock is None:
            lock = threading.Lock()
            fasta_append_locks[fasta] = lock

        with lock:
            if mode == "range":
                if is_tblastn:
                    _, _, sstart, send, sframe, sseq = payload
                    aa = sseq.replace("-", "").replace("*", "").strip()
                    if not aa:
                        log_console(f"Warning: empty HSP AA for {entry_id} taxid {taxid}")
                        return local_log_lines, False

                    start = min(sstart, send)
                    end = max(sstart, send)
                    species_name = resolve_species_name(hit_db, entry_id, taxid)
                    header = clean_header(entry_id, species_name)

                    with open(out_fasta, "a") as out_h:
                        out_h.write(header + "\n")
                        for i in range(0, len(aa), 60):
                            out_h.write(aa[i:i+60] + "\n")

                    local_log_lines.append(
                        f"Appended RANGE-AA (tblastn HSP) for taxid {taxid} "
                        f"({species_name}): {entry_id} HSP {start}-{end} "
                        f"frame {sframe} from {db_label(hit_db)} (bitscore {score})"
                    )
                    local_log_lines.append(f"[APPENDED_DB] {db_label(hit_db)}")
                    return local_log_lines, True

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
                        return local_log_lines, False

                    lines = proc.stdout.splitlines()
                    species_name = resolve_species_name(hit_db, entry_id, taxid)
                    if lines and lines[0].startswith(">"):
                        lines[0] = clean_header(entry_id, species_name)

                    with open(out_fasta, "a") as out_h:
                        out_h.write("\n".join(lines) + ("\n" if lines else ""))

                    local_log_lines.append(
                        f"Appended RANGE-DNA (blastn HSP) for taxid {taxid} "
                        f"({species_name}): {entry_id} HSP {start}-{end} "
                        f"strand {strand} from {db_label(hit_db)} (bitscore {score})"
                    )
                    local_log_lines.append(f"[APPENDED_DB] {db_label(hit_db)}")
                    return local_log_lines, True

                mode = "full"

            if mode == "full":
                if append_mode == "smart" and is_blastn:
                    contig_hit = bool(contig_like.search(entry_id))
                    transcript_hit = bool(transcript_like.search(entry_id))
                    if contig_hit and not transcript_hit:
                        log_console(f"Warning: refusing FULL fetch for contig-like ID in smart mode: {entry_id}")
                        return local_log_lines, False

                full = fetch_full_subject(hit_db, entry_id)
                if not full:
                    return local_log_lines, False

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
                                return local_log_lines, False

                            lines = proc.stdout.splitlines()
                            species_name = resolve_species_name(hit_db, entry_id, taxid)
                            if lines and lines[0].startswith(">"):
                                lines[0] = clean_header(entry_id, species_name)

                            with open(out_fasta, "a") as out_h:
                                out_h.write("\n".join(lines) + ("\n" if lines else ""))

                            local_log_lines.append(
                                f"Appended RANGE-DNA fallback for taxid {taxid} "
                                f"({species_name}): {entry_id} HSP {start}-{end} "
                                f"strand {strand} from {db_label(hit_db)} (bitscore {score})"
                            )
                            local_log_lines.append(f"[APPENDED_DB] {db_label(hit_db)}")
                            return local_log_lines, True
                        else:
                            log_console(
                                f"Warning: full too big ({bp} bp) for {entry_id} but no coords available; "
                                f"consider using blastn outfmt with coords."
                            )

                with open(out_fasta, "a") as out_h:
                    out_h.write(full if full.endswith("\n") else full + "\n")

                local_log_lines.append(
                    f"Appended FULL for taxid {taxid}: {entry_id} "
                    f"from {db_label(hit_db)} (bitscore {score})"
                )
                local_log_lines.append(f"[APPENDED_DB] {db_label(hit_db)}")
                return local_log_lines, True

        return local_log_lines, False

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
            hit_log_lines, ok = format_and_append_one_hit(fasta, hit_db, taxid, payload)
            local_log_lines.extend(hit_log_lines)
            if ok:
                appended += 1

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

    # append_one_finished_fasta is safe to run concurrently across DIFFERENT
    # fasta files: each writes to its own out_fasta/log_file, and all shared
    # state it touches (fasta_logs, appended_files, queued_for_append,
    # append_errors, species_name_cache) is already protected by locks. A
    # single append worker becomes a bottleneck once the search phase is
    # well-parallelized (e.g. --jobs 64): every appended hit does at least
    # one blastdbcmd subprocess call (a range/full fetch, plus a cached
    # taxid->species lookup), and funneling that through one thread can leave
    # most cores idle once the search phase races ahead of a single-threaded
    # append queue. --append-workers lets that stage scale independently.
    n_append_workers = max(1, int(getattr(args, "append_workers", 8)))
    is_batch_query_mode = getattr(args, "query_mode", "per-file") == "batch"

    append_threads: List[threading.Thread] = []
    if not is_batch_query_mode:
        append_threads = [
            threading.Thread(
                target=append_worker,
                name=f"blast-species-append-worker-{i}",
                daemon=True,
            )
            for i in range(n_append_workers)
        ]
        for t in append_threads:
            t.start()

    try:
        with TickingStatus(
            console,
            compute_render,
            interval=0.15,
            refresh_per_second=10,
            transient=True,
        ) as ticking_status:
            status = ticking_status

            if getattr(args, "query_mode", "per-file") == "batch":
                batch_work_dir = tempfile.mkdtemp(prefix="forge_batch_query_")
                combined_query_path = os.path.join(batch_work_dir, "combined_queries.fasta")

                # Persistent (NOT cleaned up automatically) per-db SEARCH
                # results, so a crash mid-run doesn't lose already-completed
                # BLAST searches. --resume reuses any db's result here
                # instead of re-running BLAST for it.
                batch_results_dir = "forge_batch_results"
                os.makedirs(batch_results_dir, exist_ok=True)
                use_batch_resume = bool(getattr(args, "resume", False))

                log_console(
                    f"[DEBUG] batch-query mode: combining {len(files_to_run)} "
                    f"query file(s) into one multi-FASTA for {len(db_list)} db job(s)"
                )
                qseqid_to_fasta = build_combined_query_fasta(files_to_run, combined_query_path)

                outfmt_batch = outfmt + " qseqid"

                if args.jobs is None:
                    max_workers = max(1, min(len(db_list), os.cpu_count() or 1))
                else:
                    max_workers = max(1, int(args.jobs))

                # Each db's hits get appended to their files the moment that
                # db's search returns -- there's no need to wait for every
                # other db first, since each db is an independent genome and
                # every db's best hit for a file is appended regardless of
                # what any other db found (this is exactly what --append
                # range already does in per-file mode: one appended sequence
                # per db, not a single cross-db winner). Per-(fasta, db)
                # append tracking (not just per-db SEARCH tracking) means a
                # crash loses only whatever hadn't been appended yet, down to
                # the individual (file, db) pair -- not a whole run's worth
                # of completed searches, and not even a whole db's worth of
                # appends.
                for fasta in files_to_run:
                    total_appended_by_fasta[fasta] = 0
                already_appended_dbs_by_fasta: Dict[str, Set[str]] = {}
                if use_batch_resume:
                    n_partially_done = 0
                    for fasta in files_to_run:
                        done_dbs = load_already_appended_dbs(fasta)
                        if done_dbs:
                            already_appended_dbs_by_fasta[fasta] = done_dbs
                            total_appended_by_fasta[fasta] = len(done_dbs)
                            n_partially_done += 1
                    if n_partially_done:
                        log_console(
                            f"[DEBUG] batch-query mode: {n_partially_done} file(s) have "
                            f"some hits already appended from a previous run and will "
                            f"only get NEW dbs appended, not re-appended"
                        )

                    n_cached = sum(
                        1 for db in db_list
                        if os.path.exists(os.path.join(batch_results_dir, f"{db_label(db)}.tsv"))
                    )
                    log_console(
                        f"[DEBUG] batch-query mode: {n_cached}/{len(db_list)} db(s) "
                        f"already have a cached search result and will be re-parsed, not re-searched"
                    )

                log_console(
                    f"[DEBUG] batch-query mode: scheduling {len(db_list)} db job(s) "
                    f"with max_workers={max_workers}, db_threads={multi_db_threads}"
                )

                # Dedicated append pool for batch mode: hits get queued as
                # each db finishes (not once at the very end), so this needs
                # to be running throughout the search phase, not started
                # after it like per-file mode's queue_ready_fasta_for_append.
                batch_append_queue: queue.Queue = queue.Queue()

                def batch_append_worker() -> None:
                    while True:
                        item = batch_append_queue.get()
                        try:
                            if item is None:
                                return
                            b_fasta, b_hit_db, b_taxid, b_payload = item
                            hit_log_lines, ok = format_and_append_one_hit(
                                b_fasta, b_hit_db, b_taxid, b_payload
                            )
                            with append_state_lock:
                                fasta_logs[b_fasta].extend(hit_log_lines)
                                if ok:
                                    total_appended_by_fasta[b_fasta] = (
                                        total_appended_by_fasta.get(b_fasta, 0) + 1
                                    )
                            write_one_fasta_log(b_fasta)
                        except Exception as e:
                            with append_state_lock:
                                append_errors.append((str(item), str(e)))
                            log_console(f"Error while appending batch hit {item}: {e}")
                        finally:
                            batch_append_queue.task_done()

                batch_append_threads = [
                    threading.Thread(
                        target=batch_append_worker,
                        name=f"blast-species-batch-append-worker-{i}",
                        daemon=True,
                    )
                    for i in range(n_append_workers)
                ]
                for t in batch_append_threads:
                    t.start()

                with ThreadPoolExecutor(max_workers=max_workers) as ex:
                    futures = {
                        ex.submit(
                            search_one_db_batch,
                            db,
                            combined_query_path,
                            args.tl,
                            multi_db_threads,
                            outfmt_batch,
                            args.extra,
                            os.path.join(batch_results_dir, f"{db_label(db)}.tsv"),
                            use_batch_resume,
                        ): db
                        for db in db_list
                    }

                    # Batch mode has no meaningful per-file granularity until
                    # a DB's single combined search returns, so mark every
                    # file "active" against every DB up front.
                    with state_lock:
                        for fasta in files_to_run:
                            active_db_counts[fasta] += len(db_list)
                    update_blast_progress(force=True)

                    for fut in as_completed(futures):
                        hit_db = futures[fut]
                        try:
                            _, best_by_qseqid, worker_logs = fut.result()
                        except Exception as e:
                            log_console(f"Batch search failed for db={db_label(hit_db)}: {e}")
                            with state_lock:
                                for fasta in files_to_run:
                                    active_db_counts[fasta] -= 1
                                    completed_jobs_by_fasta[fasta] += 1
                                    total_jobs_done += 1
                            update_blast_progress(force=True)
                            continue

                        # Reduce per-qseqid hits down to per-fasta: a file may
                        # contribute more than one query sequence, so keep the
                        # single best -- this matches per-file mode's
                        # semantics of "best hit across all of this file's
                        # queries against this db" exactly.
                        best_by_fasta_for_this_db: Dict[str, tuple] = {}
                        if best_by_qseqid:
                            for qseqid, (taxid, payload) in best_by_qseqid.items():
                                owning_fasta = qseqid_to_fasta.get(qseqid)
                                if owning_fasta is None:
                                    continue
                                score = payload[1]
                                cur = best_by_fasta_for_this_db.get(owning_fasta)
                                if cur is None or score > cur[1][1]:
                                    best_by_fasta_for_this_db[owning_fasta] = (taxid, payload)

                        with append_state_lock:
                            for fasta in files_to_run:
                                fasta_logs[fasta].extend(worker_logs)

                        with state_lock:
                            for fasta in files_to_run:
                                active_db_counts[fasta] -= 1
                                completed_jobs_by_fasta[fasta] += 1
                                total_jobs_done += 1

                        # Queue this db's hits for immediate appending --
                        # skipping any (fasta, db) pair already appended in a
                        # prior run, per that fasta's own log manifest.
                        db_label_str = db_label(hit_db)
                        for fasta, (taxid, payload) in best_by_fasta_for_this_db.items():
                            if db_label_str in already_appended_dbs_by_fasta.get(fasta, ()):
                                continue
                            batch_append_queue.put((fasta, hit_db, taxid, payload))

                        update_blast_progress(force=True)

                # Drain all queued hit-appends before declaring any file
                # finished -- otherwise the final summary below could run
                # ahead of appends still in flight.
                batch_append_queue.join()
                for _ in batch_append_threads:
                    batch_append_queue.put(None)
                for t in batch_append_threads:
                    t.join()

                shutil.rmtree(batch_work_dir, ignore_errors=True)

                # Every file's hits are now fully and permanently appended
                # (whether from this run or reused from a previous one).
                # Write the true completion marker per file -- the same
                # phrase file_already_finished() looks for -- so a FUTURE
                # --resume run correctly recognizes these as done.
                for fasta in files_to_run:
                    n = total_appended_by_fasta.get(fasta, 0)
                    out_fasta = fasta_info[fasta]["out_fasta"]
                    with append_state_lock:
                        if n == 0:
                            fasta_logs[fasta].append(f"No BLAST hits found in {fasta}")
                        fasta_logs[fasta].append(f"Appended {n} sequence(s) to {out_fasta}")
                    write_one_fasta_log(fasta)
                    log_console(f'Finished "{fasta}": appended {n} sequence(s)')

            elif len(db_list) == 1:
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

            if not is_batch_query_mode:
                # Per-file mode only: queue_ready_fasta_for_append/append_queue/
                # append_threads are the per-file append mechanism, drained here.
                # Batch mode has its own separate append queue+workers, already
                # fully drained and joined inside the batch branch above --
                # redundantly touching queued_for_append/appended_files here
                # would just add fastas to queued_for_append with nothing ever
                # marking them in appended_files (since no append_worker runs
                # in batch mode), hanging this wait loop forever.
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
        for _ in append_threads:
            append_queue.put(None)
        for t in append_threads:
            t.join()

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
    sp_blast.add_argument(
        "--query-mode",
        choices=["per-file", "batch"],
        default="per-file",
        help=("'per-file' (default): one BLAST process per (query, DB) pair "
              "-- simple, but for many queries against few DBs, most of the "
              "time goes to repeated per-query DB-load/setup overhead. "
              "'batch': combine ALL query files into one multi-FASTA and run "
              "just ONE BLAST process per DB, searching every query at once. "
              "Much less process/DB-load overhead for large query counts "
              "against a handful of DBs (the common case for this tool), at "
              "the cost of not seeing per-file completion until every DB "
              "has finished (results are only attributed back to individual "
              "files after each DB's single combined search returns). "
              "In batch mode, --jobs means concurrent DB searches (not "
              "concurrent query-DB pairs), and --db-threads means threads "
              "for that one big per-DB search (so it can usefully be much "
              "higher than in per-file mode).")
    )
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
    sp_blast.add_argument(
        "--append-workers",
        type=int,
        default=8,
        help=("Concurrent workers for the post-search append step (fetching "
              "hit sequences and resolving species names via blastdbcmd). "
              "This runs independently of --jobs/--db-threads, which only "
              "control the BLAST search itself. Increase this if CPU "
              "utilization during a run stays low despite a high --jobs "
              "value -- a single-threaded append stage can bottleneck an "
              "otherwise well-parallelized search phase. Default: 8.")
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
