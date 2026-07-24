"""The `forge miniprot-recover` subcommand: recover missing taxa into
per-orthogroup FASTAs by aligning all-taxa OG protein queries against raw
genome assemblies with miniprot (splice-aware protein-to-genome alignment),
instead of forge blast-species' tblastn route.

Why this exists alongside blast-species: tblastn's -max_hsps 1 caps a hit
at a single HSP, which on an unannotated multi-exon gene is a single exon
-- typically a small fraction of the full protein, and a DIFFERENT
fraction per taxon depending on which exon happened to score best. If your
downstream analysis measures gene tree discordance (e.g. to distinguish
hemiplasy from true homoplasy), that taxon-dependent fragmentation can
distort branch-level concordance estimates. miniprot reconstructs the full
spliced coding sequence in one alignment, avoiding this.

Output convention matches blast-species: original OG FASTAs are never
modified. A copy is made at '<OG>.extended<ext>' on first append, and
recovered sequences are appended there -- so both tools can accumulate
species into the SAME extended file across separate runs.
"""
from __future__ import annotations

import argparse
import collections
import os
import shutil
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

from rich.console import Console
from rich.text import Text

from ._progress import TickingStatus

GENOME_EXTS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
DONE_MARKER = "##MINIPROT_RECOVER_DONE"
RATE_WINDOW_SECONDS = 45.0


def _fmt_elapsed(seconds: float) -> str:
    total = int(round(seconds))
    h, rem = divmod(total, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _find_genomes(genomes_dir: str):
    found = []
    for p in sorted(Path(genomes_dir).iterdir()):
        if not p.is_file():
            continue
        for ext in GENOME_EXTS:
            if p.name.endswith(ext):
                found.append((p.name[: -len(ext)], p))
                break
    return found


def _build_og_query_fasta(ogseq_dir: str, ext: str, out_path: str):
    """Concatenate every '<OG><ext>' file in ogseq_dir into one miniprot
    query FASTA, prefixing each header with its OG id ('OG0001234|<id>')
    so hits can be attributed back to an orthogroup after searching. Skips
    any file whose stem already contains '.extended', so a re-run doesn't
    fold previously-recovered sequences back in as new queries.
    """
    files = sorted(Path(ogseq_dir).glob(f"*{ext}"))
    n_ogs = n_seqs = 0
    with open(out_path, "w") as out_h:
        for path in files:
            stem = path.name[: -len(ext)] if path.name.endswith(ext) else path.stem
            if ".extended" in stem:
                continue
            wrote = False
            with open(path, errors="replace") as in_h:
                for line in in_h:
                    if line.startswith(">"):
                        orig_id = line[1:].strip().split()[0]
                        out_h.write(f">{stem}|{orig_id}\n")
                        wrote = True
                        n_seqs += 1
                    else:
                        out_h.write(line)
            if wrote:
                n_ogs += 1
    return n_ogs, n_seqs


def _parse_paf_with_trans(path):
    """miniprot --trans emits one '##STA <protein>' line per alignment,
    interleaved with PAF lines in matching order."""
    paf_rows, sta_seqs = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##STA"):
                parts = line.split(None, 1)
                sta_seqs.append(parts[1].strip() if len(parts) > 1 else "")
                continue
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 12:
                continue
            paf_rows.append(f)

    if not paf_rows:
        return
    if len(sta_seqs) != len(paf_rows):
        print(
            f"WARNING: {path}: {len(paf_rows)} PAF record(s) but {len(sta_seqs)} "
            f"##STA line(s). Was miniprot run with --trans? Skipping.",
            file=sys.stderr,
        )
        return

    for f, seq in zip(paf_rows, sta_seqs):
        try:
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            tname, tstart, tend = f[5], int(f[7]), int(f[8])
            nmatch, blocklen = int(f[9]), int(f[10])
        except ValueError:
            continue

        score = 0
        for tag in f[12:]:
            if tag.startswith("AS:i:"):
                try:
                    score = int(tag[5:])
                except ValueError:
                    pass
                break

        yield {
            "qname": f[0],
            "og": f[0].split("|")[0],
            "cov": (qend - qstart) / qlen if qlen else 0.0,
            "ident": nmatch / blocklen if blocklen else 0.0,
            "tname": tname,
            "tstart": min(tstart, tend),
            "tend": max(tstart, tend),
            "score": score,
            "seq": seq.replace("*", "").replace("-", "").strip(),
        }


def _overlaps(a, b, min_frac=0.5) -> bool:
    if a["tname"] != b["tname"]:
        return False
    lo, hi = max(a["tstart"], b["tstart"]), min(a["tend"], b["tend"])
    if hi <= lo:
        return False
    shorter = min(a["tend"] - a["tstart"], b["tend"] - b["tstart"])
    return shorter > 0 and (hi - lo) / shorter >= min_frac


def _wrapped(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


class _ExtendedWriter:
    """Manages '<OG>.extended<ext>' files shared across genomes and worker
    threads: creates one from the original OG FASTA on first write, and
    tracks which species have already been appended into it via a
    companion '<OG>.miniprot_recover.log' -- so --resume never re-appends
    the same species twice, and two genomes racing to append into the
    SAME orthogroup (a real case: any species can be each other's best
    locus-conflict winner) never interleave their writes.
    """

    def __init__(self, ogseq_dir: str, ext: str):
        self.dir = Path(ogseq_dir)
        self.ext = ext
        self._locks: dict = defaultdict(threading.Lock)
        self._locks_guard = threading.Lock()

    def _lock_for(self, og: str) -> threading.Lock:
        with self._locks_guard:
            return self._locks[og]

    def _log_path(self, og: str) -> Path:
        return self.dir / f"{og}.miniprot_recover.log"

    def _extended_path(self, og: str) -> Path:
        return self.dir / f"{og}.extended{self.ext}"

    def _orig_path(self, og: str) -> Path:
        return self.dir / f"{og}{self.ext}"

    def already_appended(self, og: str, species: str) -> bool:
        lp = self._log_path(og)
        if not lp.exists():
            return False
        marker = f"[APPENDED_SPECIES] {species}\n"
        try:
            return marker in lp.read_text()
        except OSError:
            return False

    def append(self, og: str, species: str, header: str, seq: str):
        lock = self._lock_for(og)
        with lock:
            ext_path = self._extended_path(og)
            if not ext_path.exists():
                orig = self._orig_path(og)
                if not orig.exists():
                    return False, "original OG file not found"
                shutil.copy(orig, ext_path)
            with open(ext_path, "a") as fh:
                fh.write(f"{header}\n{_wrapped(seq)}\n")
            with open(self._log_path(og), "a") as lh:
                lh.write(f"[APPENDED_SPECIES] {species}\n")
            return True, "ok"


def _process_one_genome(
    stem, genome_path, query_path, work_dir, ogseq_dir, ext,
    threads, outn, min_cov, min_ident, overlap_frac,
    miniprot_bin, resume, writer, log_lock, log_path,
):
    def log(msg: str) -> None:
        line = f"[{_now()}] {stem}: {msg}"
        with log_lock:
            with open(log_path, "a") as lh:
                lh.write(line + "\n")

    paf_dir = Path(work_dir) / "paf"
    paf_dir.mkdir(parents=True, exist_ok=True)
    index_dir = Path(work_dir) / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    paf_path = paf_dir / f"{stem}.paf"
    index_path = index_dir / f"{stem}.mpi"

    def paf_done(p: Path) -> bool:
        if not p.exists():
            return False
        with open(p, "rb") as fh:
            fh.seek(max(0, os.path.getsize(p) - 128))
            return DONE_MARKER in fh.read().decode(errors="replace")

    if not (resume and paf_done(paf_path)):
        if not index_path.exists():
            log(f"indexing ({genome_path.name}) -> {index_path.name}")
            try:
                subprocess.run(
                    [miniprot_bin, "-t", str(threads), "-d", str(index_path), str(genome_path)],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                )
            except FileNotFoundError:
                log("ERROR: miniprot not found on PATH")
                return stem, False, 0
            except subprocess.CalledProcessError as e:
                msg = e.stderr.decode(errors="replace")[:500] if e.stderr else str(e)
                log(f"ERROR: indexing failed: {msg}")
                return stem, False, 0

        log(f"searching ({threads} threads, outn={outn})")
        tmp = paf_path.with_name(paf_path.name + f".tmp.{os.getpid()}")
        cmd = [miniprot_bin, "-t", str(threads), "-I", "--trans", "--outn", str(outn),
               str(index_path), query_path]
        try:
            with open(tmp, "w") as out_h:
                subprocess.run(cmd, check=True, stdout=out_h, stderr=subprocess.PIPE)
        except FileNotFoundError:
            log("ERROR: miniprot not found on PATH")
            tmp.unlink(missing_ok=True)
            return stem, False, 0
        except subprocess.CalledProcessError as e:
            msg = e.stderr.decode(errors="replace")[:500] if e.stderr else str(e)
            log(f"ERROR: search failed: {msg}")
            tmp.unlink(missing_ok=True)
            return stem, False, 0
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        with open(tmp, "a") as out_h:
            out_h.write(f"{DONE_MARKER} {_now()}\n")
        os.replace(tmp, paf_path)
    else:
        log("PAF already complete, reusing")

    hits = [
        h for h in _parse_paf_with_trans(paf_path)
        if h["seq"] and h["cov"] >= min_cov and h["ident"] >= min_ident
    ]

    best_by_og = {}
    for h in hits:
        cur = best_by_og.get(h["og"])
        if cur is None or h["score"] > cur["score"]:
            best_by_og[h["og"]] = h

    by_contig = defaultdict(list)
    for h in best_by_og.values():
        by_contig[h["tname"]].append(h)

    kept, n_conflict = [], 0
    for group in by_contig.values():
        group.sort(key=lambda h: -h["score"])
        claimed = []
        for h in group:
            if any(_overlaps(h, c, overlap_frac) for c in claimed):
                n_conflict += 1
                continue
            claimed.append(h)
        kept.extend(claimed)

    n_appended = n_skipped_resume = n_missing_og = 0
    for h in kept:
        if resume and writer.already_appended(h["og"], stem):
            n_skipped_resume += 1
            continue
        header = f">{h['tname']}_{h['tstart']}_{stem}"
        ok, _status = writer.append(h["og"], stem, header, h["seq"])
        if ok:
            n_appended += 1
        else:
            n_missing_og += 1

    log(
        f"done: {len(hits)} alignments passed filters, {len(best_by_og)} OG(s), "
        f"{n_conflict} dropped as locus conflicts, {n_appended} appended, "
        f"{n_skipped_resume} already-done (resume), {n_missing_og} missing OG file"
    )
    return stem, True, n_appended


def run_miniprot_recover(args: argparse.Namespace) -> None:
    console = Console(stderr=True)
    run_start = time.time()

    ogseq_dir = args.ogseq_dir
    ext = args.e
    work_dir = args.work_dir
    os.makedirs(work_dir, exist_ok=True)

    query_path = args.query or os.path.join(work_dir, "combined_queries.faa")
    if args.query and not os.path.exists(args.query):
        sys.exit(f"--query given but not found: {args.query}")
    if not args.query:
        console.print(f"Building combined OG query FASTA from {ogseq_dir} -> {query_path}")
        n_ogs, n_seqs = _build_og_query_fasta(ogseq_dir, ext, query_path)
        console.print(f"  {n_ogs} orthogroup(s), {n_seqs} query sequence(s)")
        if n_ogs == 0:
            sys.exit(f"No '*{ext}' files found in {ogseq_dir} to build a query from.")

    genomes = _find_genomes(args.genomes_dir)
    if args.only:
        wanted = set(args.only)
        genomes = [(s, p) for s, p in genomes if s in wanted]
        missing = wanted - {s for s, _ in genomes}
        if missing:
            sys.exit(f"--only requested stem(s) not found as genome files: {sorted(missing)}")
    if not genomes:
        sys.exit(f"No genome FASTA files found in {args.genomes_dir} (extensions: {GENOME_EXTS})")

    cpu_n = os.cpu_count() or 1
    jobs = args.jobs if args.jobs is not None else max(1, cpu_n // max(1, args.t))
    if jobs * args.t > cpu_n:
        console.print(
            f"[yellow]Warning:[/yellow] --jobs {jobs} x -t {args.t} = {jobs * args.t} "
            f"exceeds {cpu_n} detected CPUs; runs will oversubscribe."
        )

    def paf_done(p: Path) -> bool:
        if not p.exists():
            return False
        with open(p, "rb") as fh:
            fh.seek(max(0, os.path.getsize(p) - 128))
            return DONE_MARKER in fh.read().decode(errors="replace")

    paf_dir = Path(work_dir) / "paf"
    n_already_done = 0
    if args.resume:
        for stem, _path in genomes:
            if paf_done(paf_dir / f"{stem}.paf"):
                n_already_done += 1

    writer = _ExtendedWriter(ogseq_dir, ext)
    log_lock = threading.Lock()
    batch_log = os.path.join(work_dir, "miniprot_recover.log")

    console.print(
        f"[DEBUG] n_genomes={len(genomes)} threads={args.t} jobs={jobs} "
        f"outn={args.outn} min_cov={args.min_cov} min_ident={args.min_ident}"
    )
    console.print(
        f"{len(genomes)} genome(s) being processed"
        + (f", {n_already_done} already finished and skipped" if args.resume else "")
    )

    results = []
    done_count = 0
    total_appended_so_far = 0
    state_lock = threading.Lock()
    rate_samples: "collections.deque[tuple[float, int]]" = collections.deque()

    def estimate_eta(done: int, total: int) -> tuple[str, float]:
        if done <= 0:
            return "estimating...", 0.0
        now = time.time()
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
        rate = windowed_rate if windowed_rate is not None else cumulative_rate
        remaining = max(total - done, 0)
        eta_seconds = remaining / rate if rate > 0 else 0.0
        return _fmt_elapsed(eta_seconds), rate

    def render(spin: str) -> Text:
        with state_lock:
            d, appended = done_count, total_appended_so_far
        eta_str, rate = estimate_eta(d, len(genomes))
        elapsed_str = _fmt_elapsed(time.time() - run_start)
        return Text(
            f"{spin} miniprot-recover | {d}/{len(genomes)} genomes done | "
            f"{rate:.2f} genomes/s | ETA {eta_str} | Elapsed {elapsed_str} | "
            f"{appended} sequence(s) appended so far"
        )

    with TickingStatus(console, render, interval=0.15, refresh_per_second=10, transient=True):
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {
                ex.submit(
                    _process_one_genome, stem, path, query_path, work_dir, ogseq_dir, ext,
                    args.t, args.outn, args.min_cov, args.min_ident, args.overlap_frac,
                    args.miniprot_bin, args.resume, writer, log_lock, batch_log,
                ): stem
                for stem, path in genomes
            }
            for fut in as_completed(futures):
                stem = futures[fut]
                try:
                    result = fut.result()
                except Exception as e:
                    result = (stem, False, 0)
                    with open(batch_log, "a") as lh:
                        lh.write(f"[{_now()}] {stem}: unhandled exception: {e}\n")
                results.append(result)
                with state_lock:
                    done_count += 1
                    if result[1]:
                        total_appended_so_far += result[2]

    elapsed = time.time() - run_start
    n_ok = sum(1 for r in results if r[1])
    n_fail = len(results) - n_ok
    total_appended = sum(r[2] for r in results if r[1])
    console.print(
        f"Finished {len(results)} genome(s) in {elapsed / 60:.1f} min: "
        f"{n_ok} ok, {n_fail} failed, {total_appended} sequence(s) appended overall"
    )
    console.print(f"Per-genome detail: {batch_log}")
    if n_fail:
        console.print("[red]Failed:[/red]")
        for r in results:
            if not r[1]:
                console.print(f"  {r[0]}")
        sys.exit(1)


def register(sub, formatter_class) -> None:
    """Register the `miniprot-recover` subcommand."""
    sp = sub.add_parser(
        "mp-recover",
        help=(
            "Recover missing taxa into per-OG FASTAs via miniprot "
            "(splice-aware alternative to blast-species for unannotated genomes).\n"
            "(run 'forge mp-recover -h' for options)"
        ),
        epilog=(
            "PREREQUISITES\n"
            "-------------\n"
            "- miniprot on PATH (https://github.com/lh3/miniprot)\n"
            "- --ogseq-dir containing one all-taxa FASTA per orthogroup "
            "(e.g. Orthogroup_Sequences/all_taxa_fas from 'forge filter')\n"
            "- --genomes-dir containing one raw genome FASTA per taxon to recover, "
            "named as the species label you want in headers (e.g. Gekko_gekko.fa)\n\n"
            "OUTPUT\n"
            "------\n"
            "Original OG FASTAs are never modified. Recovered sequences are appended "
            "into '<OG>.extended<ext>' (created from the original on first write), "
            "the same convention forge blast-species uses -- so both tools can "
            "accumulate species into the same extended file across separate runs.\n"
        ),
        formatter_class=formatter_class,
    )
    sp.add_argument("-e", metavar="EXTENSION", default=".fa",
                     help="OG FASTA extension in --ogseq-dir (default: .fa)")
    sp.add_argument("--ogseq-dir", required=True, metavar="DIR",
                     help="directory with one all-taxa FASTA per orthogroup")
    sp.add_argument("--genomes-dir", required=True, metavar="DIR",
                     help="directory of raw genome FASTA files, one per species to recover")
    sp.add_argument("--work-dir", default="forge_miniprot_work", metavar="DIR",
                     help="working directory for indexes, PAFs, and logs")
    sp.add_argument("--query", default=None, metavar="FASTA",
                     help="use this pre-built combined query FASTA instead of building one "
                          "from --ogseq-dir (must have 'OGid|seqid' headers)")
    sp.add_argument("--miniprot-bin", default="miniprot", help="path to miniprot binary")
    sp.add_argument("-t", metavar="THREADS", type=int, default=32,
                     help="threads per miniprot process (default: 32)")
    sp.add_argument("--jobs", type=int, default=None,
                     help="concurrent genomes to process (default: cpu_count() // -t)")
    sp.add_argument("--outn", type=int, default=3, help="miniprot --outn (default: 3)")
    sp.add_argument("--min-cov", type=float, default=0.5,
                     help="minimum fraction of the query protein aligned (default: 0.5)")
    sp.add_argument("--min-ident", type=float, default=0.4,
                     help="minimum residue identity (default: 0.4)")
    sp.add_argument("--overlap-frac", type=float, default=0.5,
                     help="genomic overlap fraction treated as the same locus, for "
                          "resolving orthogroups that both claim one region (default: 0.5)")
    sp.add_argument("--resume", action="store_true",
                     help="skip genomes whose PAF is already complete, and skip "
                          "(OG, species) pairs already appended into an extended file")
    sp.add_argument("--only", nargs="*", default=None, metavar="STEM",
                     help="restrict to these genome stems -- handy for a first test run")
    sp.set_defaults(func=run_miniprot_recover)
