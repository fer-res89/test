"""The `forge filter` subcommand: copy all-taxa gene trees, QC
iso/paralog monophyly, copy matching FASTAs, with optional dedupe and/or
label cleaning.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from rich.console import Console
from rich.text import Text

try:
    from Bio import Phylo  # type: ignore
except Exception:
    Phylo = None  # type: ignore

from phylotools_common.fileutils import (
    ensure_dir,
    newly_created_dir_was,
    describe_and_print_new_dirs,
)
from ._common import (
    TREE_FILE_RE,
    TREE_OG_RE,
    FASTA_GLOBS,
    count_tree_files,
    count_fasta_files,
    materialize_trees_if_needed,
    species_from_species_tree,
)
from ._progress import TickingStatus
from .dedupe import run_dedupe_dir_inplace
from .clean_labels import load_taxon_filter, clean_fasta, clean_newick

_console = Console(stderr=True)


def _fmt_elapsed(seconds: float) -> str:
    total = int(round(seconds))
    h, rem = divmod(total, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def _estimate_eta(done: int, total: int, start_time: float) -> tuple[str, float]:
    if done <= 0:
        return "estimating...", 0.0
    elapsed = max(time.time() - start_time, 1e-9)
    rate = done / elapsed
    remaining = max(total - done, 0)
    eta_seconds = remaining / rate if rate > 0 else 0.0
    return _fmt_elapsed(eta_seconds), rate

def verify_expected_inputs(
    og_counts_tsv: Path,
    trees_path: Path,          # dir OR file
    ogseq_dir: Path,
    species_tree: Path | None = None,
    need_species_tree: bool = False,
) -> None:
    problems = []

    if not og_counts_tsv.is_file():
        problems.append(f"Missing Orthogroups.GeneCount.tsv: {og_counts_tsv}")

    if not trees_path.exists():
        problems.append(f"Missing trees input: {trees_path}")
    else:
        if not (trees_path.is_dir() or trees_path.is_file()):
            problems.append(f"Trees path is not a file/dir: {trees_path}")

    if not ogseq_dir.is_dir():
        problems.append(f"Missing Orthogroup_Sequences dir: {ogseq_dir}")

    if need_species_tree:
        if species_tree is None or not species_tree.is_file():
            problems.append(f"Missing species tree file (required): {species_tree}")

    if problems:
        msg = "Required inputs not found:\n  - " + "\n  - ".join(problems)
        raise FileNotFoundError(msg)

def list_all_taxa_orthogroups(og_counts_tsv: Path) -> List[str]:
    keep: List[str] = []
    with og_counts_tsv.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        _ = next(reader, None)
        for row in reader:
            if not row:
                continue
            og_id = row[0].strip()
            if not og_id or og_id.lower().startswith("total"):
                continue
            counts = row[1:]
            def is_zero(x: str) -> bool:
                x = x.strip()
                if x in {"", "NA", "NaN", "nan"}:
                    return True
                try:
                    return float(x) == 0.0
                except ValueError:
                    return True
            if any(is_zero(x) for x in counts):
                continue
            keep.append(og_id)
    return keep

def copy_gene_trees_for_ids(tree_dir: Path,
                            ids: Iterable[str],
                            dst_dir: Path) -> List[Path]:
    ensure_dir(dst_dir)
    copied: List[Path] = []

    # Build an index of OG -> list of candidate tree files
    og_to_files: Dict[str, List[Path]] = {}
    for p in tree_dir.iterdir():
        if not (p.is_file() and TREE_FILE_RE.match(p.name)):
            continue
        m = TREE_OG_RE.match(p.name)
        if not m:
            continue
        og = m.group("og")
        og_to_files.setdefault(og, []).append(p)

    # determinism: pick "first" by sorted filename
    for og, files in og_to_files.items():
        og_to_files[og] = sorted(files, key=lambda x: x.name)

    ids = list(ids)
    total = len(ids)
    done = 0
    start_time = time.time()

    def render(spin: str) -> Text:
        eta_str, rate = _estimate_eta(done, total, start_time)
        return Text(
            f"{spin} Copying gene trees | {done}/{total} | "
            f"{rate:.1f}/s | ETA {eta_str}"
        )

    with TickingStatus(_console, render, interval=0.15):
        for tid in ids:
            files = og_to_files.get(tid, [])
            if files:
                src = files[0]  # "pick the first"
                dst = dst_dir / src.name
                if not dst.exists():
                    shutil.copy2(src, dst)
                copied.append(dst)
            done += 1

    return copied

def open_tree_file(tree_file: Path):
    with tree_file.open("r") as fh:
        text = fh.read()
    if Phylo is None:
        sys.exit("BioPython is required for tree QC. Please install biopython.")
    tree = Phylo.read(str(tree_file), "newick")
    return text, tree

def tip_species(label: str) -> str:
    return (label or "").split("|", 1)[0].strip()

def is_species_multiplicity_monophyletic(tree_obj, species: str) -> bool:
    leaves = [t for t in tree_obj.get_terminals() if t.name]
    targets = [t for t in leaves if tip_species(t.name) == species]

    if len(targets) <= 1:
        return True

    mrca = tree_obj.common_ancestor(targets)
    descendants = list(mrca.get_terminals())
    return all(d.name and tip_species(d.name) == species for d in descendants)

def iter_tree_files(tree_dir: Path):
    for p in tree_dir.iterdir():
        if p.is_file() and TREE_FILE_RE.match(p.name):
            yield p

def qc_gene_trees_for_isos(tree_dir: Path, species_list: List[str]) -> Tuple[List[str], List[str]]:
    files = sorted(iter_tree_files(tree_dir))
    total = len(files)
    good: List[str] = []
    bad: List[str] = []
    done = 0
    current_name = ""
    start_time = time.time()

    # Single-threaded loop: `done`/`current_name` are only ever written by
    # this thread. TickingStatus's ticker thread only reads them for display,
    # so no lock is needed here (unlike blast-species, which has multiple
    # worker threads mutating shared counters).
    def render(spin: str) -> Text:
        eta_str, rate = _estimate_eta(done, total, start_time)
        parts = [
            f"{spin} QC'ing gene trees",
            f"{done}/{total} trees",
            f"{rate:.1f} trees/s",
            f"ETA {eta_str}",
        ]
        if current_name:
            parts.append(f"current: {current_name}")
        return Text(" | ".join(parts))

    with TickingStatus(_console, render, interval=0.15) as status:
        for p in files:
            current_name = p.name
            try:
                txt, t = open_tree_file(p)
            except Exception:
                bad.append(p.name)
            else:
                flags = [is_species_multiplicity_monophyletic(t, sp) for sp in species_list]
                (good if all(flags) else bad).append(p.name)
            done += 1

    elapsed = time.time() - start_time
    _console.print(
        f"QC'd {total} gene tree(s) in {_fmt_elapsed(elapsed)}: "
        f"{len(good)} good, {len(bad)} bad"
    )

    (tree_dir / "good_tree_files.list").write_text("\n".join(good) + "\n")
    (tree_dir / "bad_tree_files.list").write_text("\n".join(bad) + "\n")
    return good, bad

def extract_tree_id_from_filename(name: str) -> str:
    stem = Path(name).stem  # drops .txt
    if stem.endswith("_tree"):
        stem = stem[:-5]
    return stem

def copy_fastas_for_good_trees(ogseq_dir: Path,
                               good_tree_fnames: Iterable[str],
                               dst_dir: Path) -> List[Path]:
    ensure_dir(dst_dir)

    copied: List[Path] = []
    missing: List[str] = []

    tree_names = sorted(good_tree_fnames)
    total = len(tree_names)
    done = 0
    start_time = time.time()

    def render(spin: str) -> Text:
        eta_str, rate = _estimate_eta(done, total, start_time)
        return Text(
            f"{spin} Copying FASTAs | {done}/{total} | "
            f"{rate:.1f}/s | ETA {eta_str}"
        )

    with TickingStatus(_console, render, interval=0.15):
        for tree_name in tree_names:
            tid = extract_tree_id_from_filename(tree_name)

            found = False
            for ext in (".fa", ".fas", ".fasta", ".FA", ".FAS", ".FASTA"):
                cand = ogseq_dir / f"{tid}{ext}"
                if cand.is_file():
                    dst = dst_dir / cand.name
                    dst.write_bytes(cand.read_bytes())
                    copied.append(dst)
                    found = True
                    break

            if not found:
                missing.append(tid)
            done += 1

    rpt = dst_dir / "missing_fastas_for_good_trees.txt"
    if missing:
        rpt.write_text("\n".join(missing) + "\n")
        print(f"[WARN] {len(missing)} good trees had no matching FASTA. See: {rpt}")
    elif rpt.exists():
        rpt.unlink()

    return copied

def _pick_default_trees(base_dir: Path) -> Path | None:
    gene_trees_dir = base_dir / "Gene_Trees"
    resolved_txt   = base_dir / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt"
    resolved_dir   = base_dir / "Resolved_Gene_Trees"

    # Prefer the resolved multi-tree file if present
    if resolved_txt.is_file():
        return resolved_txt

    # Next best: a directory that actually has *_tree*.txt files
    if gene_trees_dir.is_dir() and count_tree_files(gene_trees_dir) > 0:
        return gene_trees_dir

    if resolved_dir.is_dir() and count_tree_files(resolved_dir) > 0:
        return resolved_dir

    # Last resort: if dirs exist but are empty, keep old behavior (maybe user will materialize elsewhere)
    if gene_trees_dir.exists():
        return gene_trees_dir
    if resolved_dir.exists():
        return resolved_dir

    return None

def _pick_default_species_tree(base_dir: Path) -> Path | None:
    cands = [
        base_dir / "Species_Tree" / "SpeciesTree_rooted.txt",
        base_dir / "Species_Tree" / "SpeciesTree.txt",
    ]
    for c in cands:
        if c.is_file():
            return c
    return None

def run_filter_pipeline(args: argparse.Namespace) -> None:
    base_dir: Path = args.base_dir

    og_counts = (args.gene_counts if args.gene_counts else
                 (base_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"))

    trees_path = (args.trees if args.trees else _pick_default_trees(base_dir))
    if trees_path is None:
        raise FileNotFoundError(f"Could not find gene trees. Use --trees PATH. (base-dir={base_dir.resolve()})")

    ogseq_dir = (args.ogseq if args.ogseq else (base_dir / "Orthogroup_Sequences"))
    if not ogseq_dir.is_dir():
        raise FileNotFoundError(f"Could not find Orthogroup_Sequences. Use --ogseq PATH. (got {ogseq_dir})")

    work_dir = (args.work_dir if args.work_dir else (base_dir / ".forge_work"))
    ensure_dir(work_dir)

    # Only require species tree if dedupe and no taxa-file
    need_species_tree = bool(args.dedupe) and not args.taxa_file
    species_tree_path = args.species_tree if args.species_tree else _pick_default_species_tree(base_dir)
    if need_species_tree and (species_tree_path is None or not species_tree_path.is_file()):
        raise FileNotFoundError("Dedupe requires --species-tree when --taxa-file is not provided "
                                "(or a default species tree under base-dir).")

    verify_expected_inputs(
        og_counts_tsv=Path(og_counts),
        trees_path=Path(trees_path),
        ogseq_dir=Path(ogseq_dir),
        species_tree=Path(species_tree_path) if species_tree_path else None,
        need_species_tree=need_species_tree,
    )

    # Materialize trees into a directory (handles Resolved_Gene_Trees.txt)
    trees_dir = materialize_trees_if_needed(
        Path(trees_path),
        work_dir=work_dir / "trees_split"
    )

    tree_dir_rel = Path("Gene_Trees/all_taxa_orthogroups")
    tree_dir = base_dir / tree_dir_rel
    fasta_out_rel = Path("Orthogroup_Sequences/all_taxa_fas")
    fasta_dir = base_dir / fasta_out_rel

    tree_dir_was_new = newly_created_dir_was(tree_dir)
    fasta_dir_was_new = newly_created_dir_was(fasta_dir)

    og_ids = list_all_taxa_orthogroups(og_counts)

    # Copy trees from trees_dir (not base_dir/Gene_Trees)
    copied_trees = copy_gene_trees_for_ids(
        tree_dir=trees_dir,
        ids=og_ids,
        dst_dir=tree_dir
    )

    species_for_qc = species_from_species_tree(species_tree_path) if species_tree_path else []
    # After you determine species_for_qc from species tree or TSV header...
    if args.taxa_file and os.path.isfile(args.taxa_file):
        with open(args.taxa_file) as f:
            species_for_qc = [ln.strip() for ln in f if ln.strip()]
        print(f"[filter] QC species list loaded from --taxa-file ({len(species_for_qc)} taxa)", file=sys.stderr)
    if not species_for_qc:
        with og_counts.open(newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None) or []
        species_for_qc = [h for h in header[1:] if h]

    good, bad = qc_gene_trees_for_isos(tree_dir, species_for_qc)
    copied_fastas = copy_fastas_for_good_trees(ogseq_dir, good, dst_dir=fasta_dir)

    processed_msg_parts: List[str] = []

    # Optional dedupe (FASTA only; in-place)
    if args.dedupe:
        taxa_list: Optional[List[str]] = None
        if args.taxa_file and os.path.isfile(args.taxa_file):
            with open(args.taxa_file) as f:
                taxa_list = [line.strip() for line in f if line.strip()]
        species_fallback = species_from_species_tree(species_tree_path) if species_tree_path else []
        n_files = run_dedupe_dir_inplace(fasta_dir, taxa_list, species_fallback)
        processed_msg_parts.append(f"deduped {n_files} FASTA file(s)")

    # Optional clean-labels for FASTA **and** GENE TREES
    if args.clean_labels:
        taxon_list, taxon_regex = load_taxon_filter(args.clean_labels)

        # FASTA files in output dir
        fasta_pats = [str(fasta_dir / g) for g in FASTA_GLOBS]
        fasta_files: List[str] = []
        for pat in fasta_pats:
            fasta_files.extend(glob.glob(pat))
        for fp in fasta_files:
            base, ext = os.path.splitext(fp)
            out = fp if args.inplace else f"{base}_cleaned{ext}"
            clean_fasta(fp, taxon_list, taxon_regex, args.keep_accession, out)

        # Gene trees in the copied tree dir
        tree_files: List[str] = [str(p) for p in tree_dir.glob("*_tree*.txt") if p.is_file()]
        for tp in tree_files:
            base, ext = os.path.splitext(tp)
            out = tp if args.inplace else f"{base}_cleaned{ext}"
            clean_newick(tp, taxon_list, taxon_regex, args.keep_accession, out)

        processed_msg_parts.append(
            f"clean-labels on FASTA + trees{' (in-place)' if args.inplace else ''}"
        )

        # ---- report & new-directories summary ----
    print("DONE.")
    print(f"Copied trees: {len(copied_trees)}")
    print(f"Good trees:   {len(good)} | Bad trees: {len(bad)}")
    print(f"FASTAs copied: {len(copied_fastas)}")
    if processed_msg_parts:
        print("Post-processing:", ", ".join(processed_msg_parts))
    else:
        print("Post-processing: none")

    new_dirs: List[Tuple[Path, str]] = []

    if tree_dir_was_new and tree_dir.exists():
        n_tree = count_tree_files(tree_dir)
        desc = f"{n_tree} all-taxa gene tree files; QC lists: good_tree_files.list, bad_tree_files.list"
        new_dirs.append((tree_dir_rel, desc))

    if fasta_dir_was_new and fasta_dir.exists():
        n_fa = count_fasta_files(fasta_dir)
        miss = fasta_dir / "missing_fastas_for_good_trees.txt"
        miss_extra = ""
        if miss.exists():
            try:
                miss_n = sum(1 for _ in miss.read_text().splitlines() if _.strip())
                miss_extra = f"; missing list: {miss.name} ({miss_n} entries)"
            except Exception:
                miss_extra = f"; missing list: {miss.name}"
        desc = f"{n_fa} FASTA files (good OGs){miss_extra}"
        new_dirs.append((fasta_out_rel, desc))

    describe_and_print_new_dirs(new_dirs)



def register(sub, formatter_class) -> None:
    """Register the `filter` subcommand."""
    sp_filter = sub.add_parser(
        "filter",
        help=("Copy all-taxa trees, QC iso/paralog monophyly, copy matching FASTAs; "
            "optional dedupe and/or clean-labels.\n"
            "(run 'forge filter -h' for options)"),
        formatter_class=formatter_class,
    )
    # Keep base-dir for convenience/defaults (backwards compatible)
    sp_filter.add_argument(
        "--base-dir", type=Path, default=Path("."),
        help=("OrthoFinder results dir (used for defaults). If you pass --trees/--ogseq/--species-tree "
            "these can be outside base-dir or absolute paths.")
    )
    #explicit inputs (dir or file)
    sp_filter.add_argument(
        "--trees", type=Path, default=None,
        help=("Gene trees input. Can be a directory of per-OG tree files (e.g. Gene_Trees/) "
            "OR a multi-tree file (e.g. Resolved_Gene_Trees/Resolved_Gene_Trees.txt). "
            "If omitted, defaults are searched under --base-dir.")
    )
    sp_filter.add_argument(
        "--ogseq", type=Path, default=None,
        help=("Orthogroup sequences directory (e.g. Orthogroup_Sequences/). "
            "If omitted, defaults are searched under --base-dir.")
    )
    sp_filter.add_argument(
        "--gene-counts", type=Path, default=None,
        help=("Orthogroups.GeneCount.tsv path. If omitted, defaults to "
            "--base-dir/Orthogroups/Orthogroups.GeneCount.tsv.")
    )
    sp_filter.add_argument(
        "--species-tree", type=Path, default=None,
        help=("Species tree file path. Only required when using --dedupe without --taxa-file. "
            "If omitted, defaults are searched under --base-dir.")
    )
    # Optional: where to write split trees if --trees is a multi-tree file
    sp_filter.add_argument(
        "--work-dir", type=Path, default=None,
        help=("Work directory for generated intermediates (e.g., split trees). "
            "Default: --base-dir/.forge_work/")
    )
    sp_filter.add_argument(
        "--dedupe", action="store_true",
        help="Deduplicate FASTAs in-place (uses --taxa-file if given; else species tree)."
    )
    sp_filter.add_argument(
        "--taxa-file", type=str, default=None,
        help="Text file with taxa (one per line) for dedupe."
    )
    sp_filter.add_argument(
        "--clean-labels", type=str, default=None,
        help="Clean labels on output FASTAs **and gene trees** using this taxa file or regex."
    )
    sp_filter.add_argument(
        "--keep-accession", action="store_true",
        help="When cleaning, keep accession prefix in IDs (e.g., XP_12345_taxon)."
    )
    sp_filter.add_argument(
        "--inplace", action="store_true",
        help="When cleaning, overwrite files instead of creating _cleaned copies."
    )
    sp_filter.set_defaults(func=run_filter_pipeline)
