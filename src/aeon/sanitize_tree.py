"""The `sanitize-tree` subcommand: normalize a Newick file (collapse unary
nodes, drop numeric internal support labels, strip branch lengths, etc.)
and optionally merge calibration labels from a separate calibration tree.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from io import StringIO
from pathlib import Path
from typing import Optional

import dendropy
from Bio import Phylo

def _mc_tip_labels(t: dendropy.Tree) -> set[str]:
    return {
        lf.taxon.label.strip()
        for lf in t.leaf_node_iter()
        if lf.taxon and lf.taxon.label
    }

def _mc_prune_to(t: dendropy.Tree, keep: set[str]) -> None:
    """Prune leaves not in 'keep' and suppress unifurcations (compat with distro DendroPy)."""
    keep = set(keep)
    for lf in list(t.leaf_node_iter()):
        lab = lf.taxon.label.strip() if (lf.taxon and lf.taxon.label) else None
        if not lab or lab not in keep:
            parent = lf.parent_node
            if parent is not None:
                parent.remove_child(lf)
    try:
        t.suppress_unifurcations(update_splits=False)
    except TypeError:
        t.suppress_unifurcations()

def _mc_leafset(node: dendropy.Node) -> frozenset[str]:
    return frozenset(
        lf.taxon.label.strip()
        for lf in node.leaf_nodes()
        if lf.taxon and lf.taxon.label
    )

def _mc_smallest_covering_node(tree: dendropy.Tree, groups: tuple[set[str], ...]):
    """Return the smallest internal node whose leafset intersects ALL groups."""
    groups = tuple({x for x in g} for g in groups if g)
    if not groups:
        return None

    best = None
    best_sz = None
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        L = _mc_leafset(nd)
        if len(L) < 2:
            continue
        if all(L & G for G in groups):
            sz = len(L)
            if best is None or sz < best_sz:
                best, best_sz = nd, sz
    return best

def _mc_internal_clade_map(t: dendropy.Tree) -> dict[frozenset[str], str]:
    """Map internal clades (frozenset of tip labels) → node.label (non-empty only)."""
    m: dict[frozenset[str], str] = {}
    for nd in t.postorder_node_iter():
        if nd.is_leaf():
            continue
        leaves = _mc_leafset(nd)
        if not leaves:
            continue
        lab = (nd.label or "").strip()
        if lab:
            m[leaves] = lab
    return m

def _mc_load_quartet_sets(tsv_path: str, name: str) -> tuple[set[str], set[str], set[str], set[str]]:
    import csv
    A = B = C = D = set()
    with open(tsv_path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if (row.get("name") or "").strip() == name:
                def s(col: str) -> set[str]:
                    return {t for t in (row.get(col) or "").split(",") if t}
                A, B, C, D = s("A_set"), s("B_set"), s("C_set"), s("D_set")
                break
    return A, B, C, D

def _mc_merge_calibrations(
    target_tree: dendropy.Tree,
    cal_tree: dendropy.Tree,
    no_overwrite: bool,
    label_prefix: str = "",
    quartet_sets: tuple[set[str], set[str], set[str], set[str]] | None = None,
) -> int:
    """
    Copy labels from cal_tree into target_tree by exact clade match (after pruning cal_tree
    to shared taxa). If quartet_sets=(A,B,C,D) is provided, *skip* copying labels for the
    AB and ABC clades (so those calibrations never appear on the gene tree).
    """

    # 1) Shared taxa; prune ONLY the calibration tree
    tgt_tips = _mc_tip_labels(target_tree)
    cal_tips = _mc_tip_labels(cal_tree)
    common = tgt_tips & cal_tips
    if not common:
        return 0

    # Important: we do NOT prune/simplify target_tree.
    _mc_prune_to(cal_tree, common)

    # 2) Build clade → label map from the pruned calibration tree
    cal_map = _mc_internal_clade_map(cal_tree)

    # 3) If quartet info is given, identify AB and ABC clades and exclude them
    if quartet_sets is not None:
        A, B, C, D = quartet_sets
        present = common  # tips actually present in this gene tree

        A1 = {x.strip() for x in (A or set()) if x.strip() in present}
        B1 = {x.strip() for x in (B or set()) if x.strip() in present}
        C1 = {x.strip() for x in (C or set()) if x.strip() in present}

        excluded: set[frozenset[str]] = set()
        if A1 and B1:
            excluded.add(frozenset(A1 | B1))          # AB
        if A1 and B1 and C1:
            excluded.add(frozenset(A1 | B1 | C1))     # ABC

        if excluded:
            cal_map = {
                leaves: lab
                for leaves, lab in cal_map.items()
                if leaves not in excluded
            }

    # 4) Copy remaining calibrations onto target_tree by exact clade match
    merged = 0
    for nd in target_tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        leaves = _mc_leafset(nd)
        if len(leaves) < 2:
            continue
        lab = cal_map.get(leaves)
        if not lab:
            continue
        if no_overwrite and (nd.label or "").strip():
            continue
        nd.label = f"{label_prefix}{lab}" if label_prefix else lab
        merged += 1

    return merged


# batch mode for sanitiza-tree helpers


def _st_load_quartet_names(tsv_path: Path) -> set[str]:
    import csv
    names: set[str] = set()
    with tsv_path.open(newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if "name" not in (rdr.fieldnames or []):
            raise SystemExit(f"[error] --quartets-tsv {tsv_path} is missing a 'name' column.")
        for row in rdr:
            nm = (row.get("name") or "").strip()
            if nm:
                names.add(nm)
    return names

def _st_find_matching_dirs(root: Path, quartet_names: set[str]) -> list[Path]:
    return [p for p in root.iterdir() if p.is_dir() and p.name in quartet_names]

def _st_iter_tree_files(base_dir: Path, ext: str):
    ext = ext if ext.startswith(".") else f".{ext}"
    for dirpath, _dirnames, filenames in os.walk(base_dir):
        dp = Path(dirpath)
        for fn in filenames:
            if fn.endswith(ext):
                yield dp / fn


# ---------- Rooted filter helpers ----------


NUMERIC_RE = re.compile(r'^[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?$')
LENGTHS_RE = re.compile(r":[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def parse_sanitize_letters(letters: str | None):
    """
    Map a compact letter string into sanitizer options.

    Letters:
      i = sanitize input   (rearrange-triplet only)
      o = sanitize output  (rearrange-triplet only)
      l = keep lengths
      s = keep numeric internal labels (support)
      p = allow polytomies
      h = add '<num_taxa> <num_trees>' header
      c = NO collapse (i.e., keep unary internal nodes)
    """
    s = (letters or "")
    return {
        "do_input":        ("i" in s),
        "do_output":       ("o" in s),
        "collapse":        ("c" not in s),
        "keep_support":    ("s" in s),
        "keep_lengths":    ("l" in s),
        "allow_polytomies":("p" in s),
        "add_header":      ("h" in s),
        "quiet":           True,
    }

def _collapse_unary(clade):
    """Collapse internal nodes with a single child."""
    while not clade.is_terminal() and len(clade.clades) == 1:
        child = clade.clades[0]
        clade.clades = child.clades
        clade.name = child.name
    for c in clade.clades:
        _collapse_unary(c)

def _drop_numeric_internal_labels(tree):
    """Remove numeric internal node labels (support values)."""
    for cl in tree.find_clades(order="preorder"):
        if cl.name and (not cl.is_terminal()) and NUMERIC_RE.match(str(cl.name).strip()):
            cl.name = None

def _sanitize_newick_text(nwk_text, strip_lengths=True, strip_trailing_numeric_root=True):
    """Strip branch lengths & trailing numeric root labels; clean whitespace."""
    if strip_lengths:
        nwk_text = LENGTHS_RE.sub("", nwk_text)
    nwk_text = re.sub(r"\s+", "", nwk_text)
    if strip_trailing_numeric_root:
        nwk_text = re.sub(r"\)(?:\s*[-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*;$", ");", nwk_text)
    return nwk_text


def sanitize_tree_file(
    in_path: Path,
    out_path: Path,
    collapse: bool = True,
    keep_support: bool = False,
    keep_lengths: bool = False,
    allow_polytomies: bool = False,
    add_header: bool = False,
    quiet: bool = True,
    # optional calibration arguments
    cal_tree: Path | None = None,
    quartets_tsv: Path | None = None,
    quartet_name: str | None = None,
    no_overwrite_cal: bool = False,
    label_prefix_cal: str = "",
):
    """
    Sanitize a Newick file on disk and write to out_path.

    If cal_tree is provided, calibration labels are merged *after* sanitization
    using the DendroPy-based _mc_merge_calibrations() helper. Any AB/ABC
    exclusions (from quartet_sets) are handled entirely inside _mc_merge_calibrations.
    """

    # --- Biopython-based sanitizing (unchanged logic) ---
    tree = Phylo.read(str(in_path), "newick")
    if collapse:
        _collapse_unary(tree.root)
    if not keep_support:
        _drop_numeric_internal_labels(tree)

    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    nwk = buf.getvalue()
    nwk = _sanitize_newick_text(nwk, strip_lengths=not keep_lengths)

    # Second pass with Biopython for bifurcation warning and header taxon count
    t2 = Phylo.read(StringIO(nwk), "newick")
    if (not t2.is_bifurcating()) and (not allow_polytomies) and (not quiet):
        print(
            "WARNING: Tree not fully bifurcating. Use -p/--allow-polytomies to allow.",
            file=sys.stderr,
        )

    # === Case 1: no calibration tree, behave exactly as before ===
    if cal_tree is None:
        if add_header:
            num_taxa = len(t2.get_terminals())
            nwk_out = f"{num_taxa} 1\n" + nwk
        else:
            nwk_out = nwk
        out_path.write_text(nwk_out)
        if not quiet:
            print(f"Done. Tips={len(t2.get_terminals())}", file=sys.stderr)
        return out_path

    # === Case 2: merge calibrations AFTER sanitization ===

    # Parse sanitized tree with DendroPy (keep whatever rooting DendroPy infers)
    target = dendropy.Tree.get(
        data=nwk,
        schema="newick",
        preserve_underscores=True,
    )

    # Read calibration tree
    cal = dendropy.Tree.get(
        path=str(cal_tree),
        schema="newick",
        preserve_underscores=True,
    )

    # Load quartet sets if provided
    quartet_sets = None
    if quartets_tsv is not None and quartet_name:
        A, B, C, D = _mc_load_quartet_sets(str(quartets_tsv), quartet_name)
        quartet_sets = (A, B, C, D)

    merged = _mc_merge_calibrations(
        target_tree=target,
        cal_tree=cal,
        no_overwrite=no_overwrite_cal,
        label_prefix=label_prefix_cal,
        quartet_sets=quartet_sets,
    )

    if not quiet:
        print(f"[sanitize] merged {merged} calibration label(s)", file=sys.stderr)

    # Write final tree from DendroPy.
    # Branch lengths were already stripped in 'nwk', so target has no lengths
    # unless keep_lengths=True.
    base_kwargs = dict(
        schema="newick",
        suppress_leaf_node_labels=False,
        suppress_internal_node_labels=False,
        suppress_edge_lengths=not keep_lengths,
    )
    try:
        nwk2 = target.as_string(unquoted_underscores=True, **base_kwargs).strip()
    except TypeError:
        nwk2 = target.as_string(**base_kwargs).strip()

    if add_header:
        num_taxa = len(target.leaf_nodes())
        nwk2 = f"{num_taxa} 1\n{nwk2}"

    out_path.write_text(nwk2)
    if not quiet:
        print(f"Done. Tips={len(target.leaf_nodes())}", file=sys.stderr)

    return out_path


# =========================
# pipeline core
# =========================


def cmd_sanitize_tree(args):
    in_arg  = args.input
    out_arg = args.output

    # Decide "extension batch" mode first:
    ext_batch = in_arg.startswith(".") and out_arg.startswith(".")
    do_batch  = bool(getattr(args, "batch", False) or ext_batch)

    # Common helpers
    def _merge_opts():
        return dict(
            collapse = not getattr(args, "no_collapse", False),
            keep_support = getattr(args, "keep_support", False),
            keep_lengths = getattr(args, "keep_lengths", False),
            allow_polytomies = getattr(args, "allow_polytomies", False),
            add_header = getattr(args, "add_header", False),
            quiet = False,
        )

    # ---------- BATCH PATH ----------
    if do_batch:
        # Required: -i and -o as extensions
        if not (in_arg.startswith(".") and out_arg.startswith(".")):
            sys.exit("[error] In batch mode, -i and -o must both be extensions (e.g., -i .rtd -o .san).")

        ext    = in_arg
        suffix = out_arg

        # Optional calibration inputs
        cal_tree = Path(args.cal_tree).resolve() if getattr(args, "cal_tree", None) else None
        if cal_tree and not cal_tree.is_file():
            sys.exit(f"[error] calibration tree not found: {cal_tree}")

        quartets_tsv = Path(args.quartets_tsv).resolve() if getattr(args, "quartets_tsv", None) else None
        if quartets_tsv and not quartets_tsv.is_file():
            sys.exit(f"[error] quartets TSV not found: {quartets_tsv}")

        # Mode A: with --quartets-tsv → walk subdirs named by TSV 'name'
        if quartets_tsv:
            quartet_names = _st_load_quartet_names(quartets_tsv)
            root = Path(".").resolve()
            hits = _st_find_matching_dirs(root, quartet_names)
            if not hits:
                print(f"[info] no subdirectories under {root} match quartet names in {quartets_tsv}")
                return

            print(f"[info] batch sanitize (quartets): ext={ext} suffix={suffix} root={root}")
            for qdir in hits:
                qname = qdir.name
                for in_path in _st_iter_tree_files(qdir, ext):
                    out_path = in_path.with_name(in_path.name + suffix)
                    print(f"[run] quartet={qname}  {in_path} -> {out_path}")
                    if getattr(args, "dry_run", False):
                        continue
                    sanitize_tree_file(
                        in_path=in_path,
                        out_path=out_path,
                        **_merge_opts(),
                        # merge calibrations AFTER sanitize:
                        cal_tree=cal_tree,
                        quartets_tsv=quartets_tsv,
                        quartet_name=qname if getattr(args, "quartet", None) is None else args.quartet,
                        no_overwrite_cal=getattr(args, "no_overwrite", False),
                        label_prefix_cal=getattr(args, "label_prefix", ""),
                    )
            return

        # Mode B: no TSV → operate flat in CWD
        cwd = Path(".").resolve()
        files = [p for p in cwd.iterdir() if p.is_file() and p.name.endswith(ext)]
        if not files:
            print(f"[info] no files ending with '{ext}' in {cwd}")
            return

        print(f"[info] batch sanitize (cwd): ext={ext} suffix={suffix} cwd={cwd}")
        for in_path in sorted(files):
            out_path = in_path.with_name(in_path.name + suffix)
            print(f"[run] {in_path} -> {out_path}")
            if getattr(args, "dry_run", False):
                continue
            sanitize_tree_file(
                in_path=in_path,
                out_path=out_path,
                **_merge_opts(),
                cal_tree=cal_tree,
                quartets_tsv=None,
                quartet_name=None,
                no_overwrite_cal=getattr(args, "no_overwrite", False),
                label_prefix_cal=getattr(args, "label_prefix", ""),
            )
        return

    # ---------- SINGLE-FILE PATH ----------
    # Only here do we treat -i/-o as actual file paths
    sanitize_tree_file(
        Path(in_arg), Path(out_arg),
        **_merge_opts(),
        cal_tree = Path(args.cal_tree) if getattr(args, "cal_tree", None) else None,
        quartets_tsv = Path(args.quartets_tsv) if getattr(args, "quartets_tsv", None) else None,
        quartet_name = getattr(args, "quartet", None),
        no_overwrite_cal = getattr(args, "no_overwrite", False),
        label_prefix_cal = getattr(args, "label_prefix", ""),
    )
    print(f"Wrote sanitized tree → {args.output}")



# =========================
# CLI
# =========================



def register(sub) -> None:
    """Register the `sanitize-tree` subcommand."""
    # sanitize trees
    p6 = sub.add_parser(
        'sanitize-tree',
        help='Sanitize a Newick file (collapse unary, drop numeric internal labels, strip lengths, etc.)'
    )
    p6.add_argument('-i', '--input', required=True, help="Input tree file (single-file mode). In batch mode, may be an extension like '.rtd'.")
    p6.add_argument('-o', '--output', required=True, help="Output file (single-file mode). In batch mode, may be a suffix like '.san'.")
    
    # short flags:
    p6.add_argument('-l', '--keep-lengths', action='store_true', help='Keep branch lengths')
    p6.add_argument('-s', '--keep-support', action='store_true', help='Keep numeric internal node labels')
    p6.add_argument('-p', '--allow-polytomies', action='store_true', help='Allow non-bifurcating trees')
    # -h is help, so use -H for header
    p6.add_argument('-H', '--add-header', action='store_true', help="Add '<num_taxa> <num_trees>' header")
    p6.add_argument('-c', '--no-collapse', action='store_true', help='Do NOT collapse unary internal nodes')
    p6.add_argument("--cal-tree",
                     help="Optional calibration tree with internal node labels to merge "
                          "AFTER sanitizing.")
    p6.add_argument("--quartets-tsv",
                     help="Quartet TSV with columns: name, A_set, B_set, C_set, D_set "
                          "to define A,B,C,D for MRCA scrubbing.")
    p6.add_argument("--quartet",
                     help="Name in --quartets-tsv specifying which quartet to use.")
    p6.add_argument("--no-overwrite", action="store_true",
                     help="Do not overwrite existing (non-empty) internal labels on the sanitized tree.")
    p6.add_argument("--label-prefix", default="",
                     help="Optional prefix to add to merged calibration labels (e.g. 'cal:').")
    p6.add_argument("--batch", action="store_true",
                 help="Batch mode: treat -i as an extension and -o as a suffix; "
                      "walk subdirectories whose names match quartets in --quartets-tsv.")
    p6.add_argument("--root", default=".", help="Root directory to search (default: current directory).")
    p6.add_argument("--tree-ext", help="Tree extension for batch mode (e.g., .rtd). If omitted and -i starts with '.', -i is used.")
    p6.add_argument("--suffix", help="Output suffix for batch mode (e.g., .san). If omitted and -o starts with '.', -o is used; else defaults to .san.")
    p6.add_argument("--dry-run", action="store_true", help="Batch: show what would run, do not write files.")

    p6.set_defaults(func=cmd_sanitize_tree)


