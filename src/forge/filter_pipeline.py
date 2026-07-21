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

# Sentinels for the single-pass monophyly scan. _MIXED marks a subtree
# containing more than one species; _UNNAMED is a pseudo-species for tips
# with no usable label (they behave like a distinct species, so a named
# species' clade that contains one is not monophyletic).
_MIXED = object()
_UNNAMED = object()


def nonmonophyletic_species(tree_obj, species_of_interest: Iterable[str]) -> Set[str]:
    """Return the subset of `species_of_interest` whose tips are NOT
    rooted-monophyletic in this tree, computed in a SINGLE post-order pass.

    Semantics are identical to calling is_species_multiplicity_monophyletic()
    once per species (a species with <= 1 tip is always fine; otherwise all
    its tips must form a rooted clade containing no tip of any other species,
    including unnamed tips) -- but O(N) per tree instead of
    O(species x N x targets), since it avoids recomputing get_terminals() and
    an MRCA path-trace for every species.

    Correctness argument: for a species s with total tip count n_s >= 2, s is
    rooted-monophyletic iff some node's subtree is "pure s" (contains only
    s-tips) AND holds all n_s of them. We track, per node, (subtree tip count,
    pure species | _MIXED), and record the largest pure-s clade size seen;
    s is monophyletic iff that largest size equals n_s.
    """
    interest = set(species_of_interest)

    total_count: Dict[object, int] = {}     # species -> total tips in tree
    max_pure_count: Dict[object, int] = {}  # species -> largest pure clade size

    count: Dict[int, int] = {}       # id(node) -> tips in subtree
    pure: Dict[int, object] = {}     # id(node) -> species name / _UNNAMED / _MIXED

    # Explicit iterative post-order (avoids Python recursion limits on very
    # large trees and BioPython's slower recursive helpers).
    root = tree_obj.root
    stack = [(root, False)]
    while stack:
        node, processed = stack.pop()
        children = node.clades
        if not processed:
            stack.append((node, True))
            for ch in children:
                stack.append((ch, False))
            continue

        if not children:
            # terminal
            sp = tip_species(node.name) if node.name else ""
            sp_key: object = sp if sp else _UNNAMED
            count[id(node)] = 1
            pure[id(node)] = sp_key
            if sp_key is not _UNNAMED:
                total_count[sp_key] = total_count.get(sp_key, 0) + 1
            continue

        c = 0
        subtree_species: object = None
        seen_first = False
        mixed = False
        for ch in children:
            c += count.pop(id(ch))
            ps = pure.pop(id(ch))
            if mixed:
                continue
            if ps is _MIXED:
                mixed = True
            elif not seen_first:
                subtree_species = ps
                seen_first = True
            elif ps != subtree_species:
                # != is value-equality for species-name strings and identity
                # for the _UNNAMED sentinel (a distinct object), which is
                # exactly what we want.
                mixed = True

        count[id(node)] = c
        if mixed:
            pure[id(node)] = _MIXED
        else:
            pure[id(node)] = subtree_species
            if subtree_species is not _UNNAMED and subtree_species is not None:
                if max_pure_count.get(subtree_species, 0) < c:
                    max_pure_count[subtree_species] = c

    bad: Set[str] = set()
    for sp in interest:
        n = total_count.get(sp, 0)
        if n <= 1:
            continue
        if max_pure_count.get(sp, 0) != n:
            bad.add(sp)
    return bad


def is_species_multiplicity_monophyletic(tree_obj, species: str) -> bool:
    """Backward-compatible single-species check. Prefer
    nonmonophyletic_species() for QC over many species -- this shim exists
    for any external callers and simply asks the batch scan about one
    species."""
    return species not in nonmonophyletic_species(tree_obj, (species,))

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
                # Single post-order pass over the whole tree, checking every
                # species at once, instead of one MRCA path-trace per species.
                violations = nonmonophyletic_species(t, species_list)
                (good if not violations else bad).append(p.name)
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

def resolve_fasta_out_dir(
    value: Optional[str],
    base_dir: Path,
    default_parent: str = "Orthogroup_Sequences",
    default_name: str = "all_taxa_fas",
) -> Path:
    """Resolve the --fasta-out value into an output directory.

    One flag, three behaviors:
      - None (default): base_dir / Orthogroup_Sequences / all_taxa_fas
      - a BARE NAME with no path separator (e.g. "my_fastas"): treated as a
        rename of the folder inside the default parent, i.e.
        base_dir / Orthogroup_Sequences / my_fastas
      - a PATH with directory components:
          * relative (e.g. "out/fastas", "../elsewhere"): resolved against
            base_dir
          * absolute (e.g. "/data/fastas"): used as-is

    This lets a user either just rename the output folder or relocate it
    entirely, without a second flag.
    """
    if value is None:
        return base_dir / default_parent / default_name

    p = Path(value)
    if p.is_absolute():
        return p
    if p.parent == Path("."):
        # bare name -> rename inside the default parent
        return base_dir / default_parent / p
    # relative path with directory components -> relative to base_dir
    return base_dir / p


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
    fasta_dir = resolve_fasta_out_dir(getattr(args, "fasta_out", None), base_dir)
    # For the "new directories" report, show a path relative to base_dir when
    # the output lives under it, otherwise show the resolved path as-is.
    try:
        fasta_out_rel = fasta_dir.relative_to(base_dir)
    except ValueError:
        fasta_out_rel = fasta_dir

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
    # Optional clean-labels for FASTA **and** GENE TREES
    if args.clean_labels:
        _phase("Cleaning labels on FASTAs and gene trees...")
        taxon_list, taxon_regex = load_taxon_filter(args.clean_labels)

        # FASTA files in output dir
        fasta_pats = [str(fasta_dir / g) for g in FASTA_GLOBS]
        fasta_files: List[str] = []
        for pat in fasta_pats:
            fasta_files.extend(glob.glob(pat))
        tree_files: List[str] = [str(p) for p in tree_dir.glob("*_tree*.txt") if p.is_file()]

        total_clean = len(fasta_files) + len(tree_files)
        done_clean = 0
        clean_start = time.time()
        clean_kind = "FASTAs"

        def _clean_render(spin: str) -> Text:
            eta_str, rate = _estimate_eta(done_clean, total_clean, clean_start)
            return Text(
                f"{spin} Cleaning labels ({clean_kind}) | "
                f"{done_clean}/{total_clean} | {rate:.1f}/s | ETA {eta_str}"
            )

        with TickingStatus(_console, _clean_render, interval=0.15):
            for fp in fasta_files:
                base, ext = os.path.splitext(fp)
                out = fp if args.inplace else f"{base}_cleaned{ext}"
                clean_fasta(fp, taxon_list, taxon_regex, args.keep_accession, out)
                done_clean += 1

            clean_kind = "gene trees"
            for tp in tree_files:
                base, ext = os.path.splitext(tp)
                out = tp if args.inplace else f"{base}_cleaned{ext}"
                clean_newick(tp, taxon_list, taxon_regex, args.keep_accession, out)
                done_clean += 1

        _console.print(
            f"[filter] Cleaned {len(fasta_files)} FASTA(s) + {len(tree_files)} tree(s)."
        )
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
    sp_filter.add_argument(
        "--fasta-out", metavar="NAME_OR_PATH", type=str, default=None,
        help=("Output directory for the copied good-OG FASTAs. Accepts EITHER "
            "a bare name (e.g. 'my_fastas'), which just renames the folder "
            "inside the default Orthogroup_Sequences/ parent, OR a path "
            "(relative to --base-dir, or absolute) to place the output "
            "anywhere. Default: Orthogroup_Sequences/all_taxa_fas")
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
