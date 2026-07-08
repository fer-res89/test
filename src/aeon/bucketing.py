"""The `bucket-files` subcommand: bucket gene trees by
quartet topology support and optionally copy matching FASTA sequences.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import dendropy

def _qb_parse_quartet_arg(qstr: str):
    # "Q1:A1,A2;B1;C1,C2;D1,D2"
    name, rest = qstr.split(":", 1)
    parts = rest.split(";")
    if len(parts) != 4:
        raise ValueError("Quartet must have four groups separated by ';'")
    def to_set(s): return set([t for t in s.split(",") if t])
    return name, to_set(parts[0]), to_set(parts[1]), to_set(parts[2]), to_set(parts[3])

def _qb_load_quartets(quartet_str: str | None, quartets_tsv: str | None):
    """
    Returns a list of dicts, each:
      {
        "name": str,
        "A": set[str], "B": set[str], "C": set[str], "D": set[str],
        "labels": {"A": str, "B": str, "C": str, "D": str} | None
      }

    TSV mode supports a companion row per quartet named '<name>_names' with columns
    A_name, B_name, C_name, D_name (or A_set/B_set/... if you prefer), e.g.:

      name    A_set           B_set     C_set   D_set
      Q1      aves1,aves2     croc1     turt1   archo1,archo2
      Q1_names A_name=aves    B_name=croc   C_name=turt   D_name=archo
    """
    out = []

    def parse_set(csv_cell: str) -> set[str]:
        return {x for x in (t.strip() for t in (csv_cell or "").split(",")) if x}

    def parse_single(qstr: str):
        nm, rest = qstr.split(":", 1)
        parts = rest.split(";")
        if len(parts) != 4:
            raise SystemExit("Quartet must have four groups separated by ';' (A;B;C;D)")
        def toks(p): return {t for t in (x.strip() for x in p.split(",")) if t}
        A, B, C, D = toks(parts[0]), toks(parts[1]), toks(parts[2]), toks(parts[3])
        return {"name": nm, "A": A, "B": B, "C": C, "D": D, "labels": None}

    if quartet_str:
        out.append(parse_single(quartet_str))

    if quartets_tsv:
        rows = []
        with open(quartets_tsv, newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            rows = list(rdr)

        # Index by name (exact match)
        by_name: dict[str, dict] = {}
        pending_labels: dict[str, dict] = {}

        for r in rows:
            name = (r.get("name") or "").strip()
            if not name:
                continue

            # companion labels row like "Q1_names"
            if name.endswith("_names"):
                base = name[:-6]
                # accept either A_name/B_name/... OR A_set/B_set/... as labels, pick *_name first
                labels = {
                    "A": (r.get("A_name") or r.get("A_set") or "").strip(),
                    "B": (r.get("B_name") or r.get("B_set") or "").strip(),
                    "C": (r.get("C_name") or r.get("C_set") or "").strip(),
                    "D": (r.get("D_name") or r.get("D_set") or "").strip(),
                }
                # keep only non-empty
                labels = {k: v for k, v in labels.items() if v}
                pending_labels[base] = labels
                continue

            # regular quartet row (expects *_set)
            need = ["A_set", "B_set", "C_set", "D_set"]
            if not all(k in r for k in need):
                raise SystemExit(f"[ERROR] TSV row for '{name}' missing columns: {need}")
            q = {
                "name": name,
                "A": parse_set(r["A_set"]),
                "B": parse_set(r["B_set"]),
                "C": parse_set(r["C_set"]),
                "D": parse_set(r["D_set"]),
                "labels": None,
            }
            by_name[name] = q

        # attach labels if we saw a <name>_names row
        for nm, q in by_name.items():
            if nm in pending_labels and pending_labels[nm]:
                q["labels"] = pending_labels[nm]
            out.append(q)

    if not out:
        raise SystemExit("[ERROR] No quartets provided.")
    return out


def _qb_collapse_by_support(tree: dendropy.Tree, thresh: float):
    """Contract internal edges whose child node label < thresh (non-numeric → 0)."""
    for edge in list(tree.postorder_edge_iter()):
        head = edge.head_node
        if head is None or head.is_leaf():
            continue
        try:
            val = float(head.label) if head.label is not None else 0.0
        except Exception:
            val = 0.0
        if val < thresh:
            parent = edge.tail_node
            if parent is None:
                continue
            for ch in list(head.child_node_iter()):
                head.remove_child(ch)
                parent.add_child(ch)
            parent.remove_child(head)

def _qb_encode_splits(tree: dendropy.Tree):
    """Return (splits, present_taxa). splits is list of (frozenset(left), frozenset(right))."""
    present = {lf.taxon.label.strip()
               for lf in tree.leaf_node_iter()
               if lf.taxon and lf.taxon.label}
    splits = []
    for e in tree.postorder_edge_iter():
        if not e.is_internal():
            continue
        head = e.head_node
        if head is None or head.is_leaf():
            continue
        left = {lf.taxon.label.strip()
                for lf in head.leaf_nodes()
                if lf.taxon and lf.taxon.label}
        right = present - left
        if left and right:
            splits.append((frozenset(left), frozenset(right)))
    return splits, present

def _qb_has_bipartition_for(top_left, top_right, splits, present_taxa):
    """True if a split separates top_left from top_right (after intersecting with present_taxa)."""
    Lwant = (set(x.strip() for x in top_left) & present_taxa)
    Rwant = (set(x.strip() for x in top_right) & present_taxa)
    if not Lwant or not Rwant:
        return False
    for L, R in splits:
        if Lwant.issubset(L) and Rwant.issubset(R):
            return True
        if Lwant.issubset(R) and Rwant.issubset(L):
            return True
    return False

def _qb_classify_tree_for_quartet(tree: dendropy.Tree, quartet, bs_thresh: float, min_per_side: int):
    """Return one of: 'AB|CD', 'AC|BD', 'AD|BC', 'unresolved', 'insufficient'."""
    name, A, B, C, D = quartet
    t = tree.clone(depth=2)
    if bs_thresh and bs_thresh > 0:
        _qb_collapse_by_support(t, bs_thresh)
    splits, present = _qb_encode_splits(t)

    for grp in (A, B, C, D):
        if len(grp & present) < min_per_side:
            return "insufficient"

    if _qb_has_bipartition_for(A|B, C|D, splits, present):
        return "AB|CD"
    if _qb_has_bipartition_for(A|C, B|D, splits, present):
        return "AC|BD"
    if _qb_has_bipartition_for(A|D, B|C, splits, present):
        return "AD|BC"
    return "unresolved"

def _qb_classify_quartet_bipartition(
    tree: dendropy.Tree,
    quartet: tuple[str, set[str], set[str], set[str], set[str]],
    bs_thresh: float,
    min_per_side: int,
) -> str:
    """
    Thin wrapper around the existing bipartition-based classification.

    quartet = (name, A, B, C, D).
    """
    return _qb_classify_tree_for_quartet(
        tree,
        quartet,
        bs_thresh,
        min_per_side,
    )


def _qb_classify_quartet_lenient(
    tree: dendropy.Tree,
    quartet: tuple[str, set[str], set[str], set[str], set[str]],
    bs_thresh: float,
    min_per_side: int,
) -> str:
    """
    Lenient clade-aware mode:

      - First classify by bipartitions (AB|CD, AC|BD, AD|BC, unresolved, insufficient)
      - If outcome is not one of the 3 topologies, return it as-is.
      - If outcome is AB|CD / AC|BD / AD|BC, additionally require:

          * A, B, C are each monophyletic relative to the union (A∪B∪C∪D)
            in this tree.

      - If monophyly fails, downgrade to "unresolved".
    """
    name, A, B, C, D = quartet
    outcome = _qb_classify_quartet_bipartition(
        tree, quartet, bs_thresh, min_per_side
    )
    if outcome not in ("AB|CD", "AC|BD", "AD|BC"):
        return outcome

    present = _qb_present_labels(tree)
    universe = (A | B | C | D) & present

    for G in (A, B, C):
        if not _qb_is_monophyletic_wrt(tree, G & present, universe):
            return "unresolved"

    return outcome


def _qb_classify_quartet_strict(
    tree: dendropy.Tree,
    quartet: tuple[str, set[str], set[str], set[str], set[str]],
    bs_thresh: float,
    min_per_side: int,
) -> str:
    """
    Strict clade-aware mode:

      1) Start from bipartition classification (as above).
      2) Require:
         - A, B, C are monophyletic relative to union(A,B,C,D)
         - The MRCA of A∪B∪C contains **no taxa** outside A∪B∪C∪D
           (i.e. ingroup = A,B,C plus possibly D; no extra outsiders inside).

      3) If either condition fails, return "unresolved".
    """
    name, A, B, C, D = quartet
    outcome = _qb_classify_quartet_bipartition(
        tree, quartet, bs_thresh, min_per_side
    )
    if outcome not in ("AB|CD", "AC|BD", "AD|BC"):
        return outcome

    present = _qb_present_labels(tree)
    universe = (A | B | C | D) & present

    # A, B, C each monophyletic wrt universe
    for G in (A, B, C):
        if not _qb_is_monophyletic_wrt(tree, G & present, universe):
            return "unresolved"

    # Ingroup clade = MRCA(A∪B∪C)
    ingroup = (A | B | C) & present
    nd_abc = _qb_mrca_or_none(tree, ingroup)
    if nd_abc is None:
        return "unresolved"

    leaves_abc = _qb_leafset(nd_abc)
    extraneous = leaves_abc - (A | B | C | D)
    if extraneous:
        # There are taxa inside the A/B/C clade that are not A/B/C/D
        return "unresolved"

    return outcome


def _qb_iter_gene_trees(gene_trees: str | None, tree_glob: list[str] | None):
    """
    Yield (gene_id, tree, src_path_or_none).
    For 'glob' mode, you can pass multiple patterns/paths (nargs='+').
    """
    if tree_glob:
        # expand every token into files; allow explicit paths
        expanded: list[str] = []
        for tok in tree_glob:
            matches = glob.glob(tok)
            if matches:
                expanded.extend(matches)
            else:
                p = Path(tok)
                if p.exists():
                    expanded.append(str(p))
                else:
                    print(f"[warn] no matches for pattern/path: {tok}", file=sys.stderr)

        files = sorted(set(expanded))
        if not files:
            print("[warn] no input files found for --tree-glob", file=sys.stderr)
            return

        for path in files:
            p = Path(path)
            try:
                tl = dendropy.TreeList.get(path=str(p), schema="newick",
                                           rooting="force-unrooted", preserve_underscores=True)
            except Exception as e:
                print(f"[warn] failed to read {path}: {e}", file=sys.stderr)
                continue
            if len(tl) == 1:
                yield p.stem, tl[0], p
            else:
                for i, t in enumerate(tl, 1):
                    yield f"{p.stem}_T{str(i).zfill(3)}", t, p
    else:
        tl = dendropy.TreeList.get(path=gene_trees, schema="newick",
                                   rooting="force-unrooted", preserve_underscores=True)
        for i, t in enumerate(tl, 1):
            yield f"T{str(i).zfill(6)}", t, None


def _qb_bucket_dir_map(q: dict, focal: str | None = None) -> dict[str, str]:
    """
    Given a quartet dict (from _qb_load_quartets), return mapping from outcomes to
    human-friendly directory names. Outcomes: 'AB|CD', 'AC|BD', 'AD|BC',
    plus 'unresolved' and 'insufficient'.

    If labels are provided (q['labels'] with keys A/B/C/D), and a focal group is given
    ('A'|'B'|'C'|'D' or the label value), we always name the directory with the pair that
    CONTAINS the focal group (e.g., '<other>_<focal>').
    If no labels, fall back to 'AB_CD', etc.

    Example:
      labels: {'A':'aves','B':'croc','C':'turt','D':'archo'}, focal='C'
      AB|CD -> choose CD -> 'archo_turt'
      AC|BD -> choose AC -> 'aves_turt'
      AD|BC -> choose BC -> 'croc_turt'
    """
    # default generic names
    generic = {
        "AB|CD": "AB_CD",
        "AC|BD": "AC_BD",
        "AD|BC": "AD_BC",
        "unresolved": "unresolved",
        "insufficient": "insufficient",
    }

    labels = q.get("labels")
    if not labels:
        return generic

    # normalize label lookup by either key ('A'/'B'/'C'/'D') or value
    def key_of(x: str) -> str | None:
        x = x.strip().lower()
        for k, v in labels.items():
            if k.lower() == x or (v and v.lower() == x):
                return k
        return None

    # pick focal; default to 'C' if not provided
    k_focal = key_of(focal) if focal else "C"
    if k_focal not in ("A","B","C","D"):
        k_focal = "C"

    def name_pair(k1, k2):
        # produce '<label_other>_<label_focal>' lowercased
        lab1 = labels.get(k1, k1).lower()
        lab2 = labels.get(k2, k2).lower()
        # keep focal second
        if k2 == k_focal:
            other, foc = lab1, lab2
        else:
            other, foc = lab2, lab1
        return f"{other}_{foc}"

    # For each topology X|Y, choose the pair that contains the focal
    # Pairs: AB|CD, AC|BD, AD|BC
    mapping = {}
    # AB|CD
    mapping["AB|CD"] = (
        name_pair("C","D") if k_focal in ("C","D") else name_pair("A","B")
    )
    # AC|BD
    mapping["AC|BD"] = (
        name_pair("A","C") if k_focal in ("A","C") else name_pair("B","D")
    )
    # AD|BC
    mapping["AD|BC"] = (
        name_pair("B","C") if k_focal in ("B","C") else name_pair("A","D")
    )
    # standard others
    mapping["unresolved"] = "unresolved"
    mapping["insufficient"] = "insufficient"
    return mapping


def _qb_find_matching_seq_by_suffix(
    tree_name: str,
    seq_dir: Path,
    seq_ext: str,
    strip_suffixes: Optional[Iterable[str]] = None,
) -> Optional[Path]:
    """
    Tries to find a matching sequence file by stripping known suffixes from tree_name.

    Example:
      tree_name = 'OG001234_aln.fa.tre.rtd'
      strip_suffixes = ['.fa.tre.rtd', '.tre.rtd']
      seq_ext = '.fas'
      -> will try OG001234_aln.fas, OG001234_aln.fa.fas, etc.
    """
    # 1) Try explicit suffixes, longest first
    if strip_suffixes:
        for suf in sorted(strip_suffixes, key=len, reverse=True):
            if tree_name.endswith(suf):
                base = tree_name[: -len(suf)]
                cand = seq_dir / f"{base}{seq_ext}"
                if cand.exists():
                    return cand

    # 2) Fallback: generic progressive stripping of dot-segments
    base = tree_name
    while True:
        cand = seq_dir / f"{base}{seq_ext}"
        if cand.exists():
            return cand
        if "." not in base:
            break
        base = base.rsplit(".", 1)[0]

    # 3) Final fallback: pure stem
    stem = Path(tree_name).stem
    cand = seq_dir / f"{stem}{seq_ext}"
    if cand.exists():
        return cand

    return None


def _qb_present_labels(tree: dendropy.Tree) -> set[str]:
    """All leaf labels in a DendroPy tree, stripped."""
    return {
        lf.taxon.label.strip()
        for lf in tree.leaf_node_iter()
        if lf.taxon and lf.taxon.label
    }


def _qb_leafset(node: dendropy.Node) -> set[str]:
    """Leaf labels below a node."""
    return {
        lf.taxon.label.strip()
        for lf in node.leaf_nodes()
        if lf.taxon and lf.taxon.label
    }


def _qb_mrca_or_none(tree: dendropy.Tree, labels: set[str]) -> dendropy.Node | None:
    """MRCA for a set of labels, or None if undefined / <2 tips."""
    labels = {x for x in labels if x}
    if len(labels) < 2:
        return None
    try:
        return tree.mrca(taxon_labels=sorted(labels))
    except Exception:
        return None


def _qb_is_monophyletic_wrt(
    tree: dendropy.Tree,
    group: set[str],
    universe: set[str],
) -> bool:
    """
    Test if 'group' is monophyletic **relative to** 'universe'.

    I.e. find MRCA(group) and check that, among labels in 'universe',
    that clade contains group and no other 'universe' labels.
    """
    group = group & universe
    if len(group) < 1:
        return False
    nd = _qb_mrca_or_none(tree, group)
    if nd is None:
        return False
    leaves = _qb_leafset(nd) & universe
    # group must be subset, and no "foreign" labels from universe
    return group.issubset(leaves) and not (leaves - group)


# ===== Calibration merge helpers (used by sanitize-tree) =====


def _rf_leafset(node: dendropy.Node) -> frozenset[str]:
    """Return the set of leaf labels under a node (stripped, no Nones)."""
    leaves = set()
    for lf in node.leaf_nodes():
        if lf.taxon and lf.taxon.label:
            leaves.add(lf.taxon.label.strip())
    return frozenset(leaves)


def _rf_present_labels(tree: dendropy.Tree) -> set[str]:
    return {
        lf.taxon.label.strip()
        for lf in tree.leaf_node_iter()
        if lf.taxon and lf.taxon.label
    }


def _rf_is_monophyletic(tree: dendropy.Tree, group: set[str]) -> bool:
    """
    True iff:
      - group ∩ present is non-empty, and
      - if |group∩present| == 1: trivially monophyletic
      - if |group∩present| >= 2: there exists an internal node whose leafset
        is exactly that set.
    """
    present = _rf_present_labels(tree)
    G = {x for x in group if x in present}
    if not G:
        # nothing from this group in the tree → cannot assess monophyly
        return False

    # singleton group is trivially monophyletic
    if len(G) == 1:
        return True

    # logic for ≥2 taxa
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        leaves = _rf_leafset(nd)
        if leaves and leaves <= present and leaves == G:
            return True

    return False

def _rf_smallest_covering_node(tree: dendropy.Tree, sets: list[set[str]]) -> dendropy.Node | None:
    """
    Smallest node whose leafset contains all the taxa in the union of 'sets'.
    Returns None if union is empty or cannot be covered (shouldn't happen if leaves exist).
    """
    present = _rf_present_labels(tree)
    target = set()
    for S in sets:
        target |= {x for x in S if x in present}
    if not target:
        return None

    best = None
    best_size = None

    for nd in tree.postorder_node_iter():
        leaves = _rf_leafset(nd)
        if not leaves:
            continue
        if not target.issubset(leaves):
            continue
        if best is None or len(leaves) < best_size:
            best = nd
            best_size = len(leaves)

    return best

def _rf_check_rooted_constraints(
    tree: dendropy.Tree,
    topology: str,
    A: set[str],
    B: set[str],
    C: set[str],
    D: set[str],
    mode: str,
    max_strays: int = 2,
) -> bool:
    """
    Rooted second-stage filter.

    We assume 'tree' is already the gene tree, possibly rooted by mammals. We only
    use rooted structure to check that A,B,C behave like coherent clades under a
    species-tree-aware view.

    Returns True if the tree passes the rooted filter for the given 'topology'
    (AB|CD, AC|BD, AD|BC), False if it should be demoted to 'root_fail'.
    """
    if mode == "bipartition":
        # No rooted filtering; accept any unrooted informative classification.
        return True

    present = _rf_present_labels(tree)
    A1 = {x for x in A if x in present}
    B1 = {x for x in B if x in present}
    C1 = {x for x in C if x in present}
    D1 = {x for x in D if x in present}

    # If we don't have any of A,B,C, classification is meaningless anyway.
    if not (A1 or B1 or C1):
        return False

    # 1) Ensure A,B,C are each monophyletic (within present tips), if any members are present.
    for label, G in (("A", A1), ("B", B1), ("C", C1)):
        if G and not _rf_is_monophyletic(tree, G):
            return False

    # 2) Find the smallest clade covering A∪B∪C.
    union_ABC = A1 | B1 | C1
    if not union_ABC:
        return False

    nd_abc = _rf_smallest_covering_node(tree, [union_ABC])
    if nd_abc is None:
        return False

    leaves_abc = _rf_leafset(nd_abc)

    # Focal tips inside ABC:
    focal_in_abc = leaves_abc & (A1 | B1 | C1)
    nonfocal_in_abc = leaves_abc - focal_in_abc

    # 3) In lenient/strict modes, allow at most max_strays non-focal taxa inside ABC.
    if len(nonfocal_in_abc) > max_strays:
        return False

    if mode == "lenient":
        # In lenient mode we don't care what those strays are (could be D or others),
        # as long as they are few enough.
        return True

    # 4) STRICT mode: disallow any taxon in ABC that is not in A∪B∪C∪D.
    allowed = A1 | B1 | C1 | D1
    bad = [t for t in leaves_abc if t not in allowed]
    if bad:
        return False

    return True


# ===== extract-dates: quartet helpers =====



def cmd_bucket_files(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    quartets = _qb_load_quartets(args.quartet, args.quartets_tsv)

    modes = set(args.out_mode)                     # e.g. {"copy_seqs"}
    do_lists = "lists" in modes
    want_copy_trees = "copy_trees" in modes
    want_copy_seqs  = "copy_seqs"  in modes

    per_file_input = bool(args.tree_glob)          # --tree-glob (nargs="+") present?
    if (want_copy_trees or want_copy_seqs) and not per_file_input:
        print("[info] copy_* modes require --tree-glob; disabling copy modes.", file=sys.stderr)
        want_copy_trees = False
        want_copy_seqs  = False

    # seq_dir first, then derive effective flags
    seq_dir = Path(args.seq_dir) if args.seq_dir else None
    do_copy_trees = want_copy_trees
    do_copy_seqs  = want_copy_seqs and (seq_dir is not None)
    if want_copy_seqs and seq_dir is None:
        print("[info] --out-mode copy_seqs requested but --seq-dir not provided; disabling copy_seqs.",
              file=sys.stderr)
        do_copy_seqs = False

    # ---- writers + bucket dirs ----
    writers: dict[str, dict[str, TextIO]] = {}
    buckets: dict[str, dict[str, Path]] = {}

    for q in quartets:
        qname = q["name"]
        qdir = outdir / qname
        qdir.mkdir(parents=True, exist_ok=True)

        # folder names for the 3 informative outcomes
        dir_map = _qb_bucket_dir_map(q, getattr(args, "focal", None)) if "_qb_bucket_dir_map" in globals() else {
            "AB|CD": "AB_CD",
            "AC|BD": "AC_BD",
            "AD|BC": "AD_BC",
            # we will not create dirs for 'unresolved', 'insufficient', 'root_fail'
        }

        # --- Writers ---
        writers[qname] = {}

        # Always open lists for UNRESOLVED / INSUFFICIENT / ROOT_FAIL
        for k in ("unresolved", "insufficient", "root_fail"):
            writers[qname][k] = open(qdir / f"{qname}.{k}.txt", "w")

        # Lists for informative outcomes only if requested
        if do_lists:
            for k, fname in (("AB|CD", "AB_CD"),
                             ("AC|BD", "AC_BD"),
                             ("AD|BC", "AD_BC")):
                writers[qname][k] = open(qdir / f"{qname}.{fname}.txt", "w")

        # --- Bucket dirs (copy destinations) ---
        if do_copy_trees or do_copy_seqs:
            buckets[qname] = {
                "AB|CD": qdir / dir_map["AB|CD"],
                "AC|BD": qdir / dir_map["AC|BD"],
                "AD|BC": qdir / dir_map["AD|BC"],
                # no dirs for unresolved/insufficient/root_fail
            }
            for p in buckets[qname].values():
                p.mkdir(parents=True, exist_ok=True)

    # initialize counts (including root_fail)
    counts = {q["name"]: defaultdict(int) for q in quartets}

    # ---- iterate inputs ----
    for gid, tree, src_path in _qb_iter_gene_trees(args.gene_trees, args.tree_glob):
        for q in quartets:
            qname = q["name"]
            A, B, C, D = q["A"], q["B"], q["C"], q["D"]

            # STEP 1: unrooted bipartition-based classification
            base_outcome = _qb_classify_tree_for_quartet(
                tree,
                (qname, A, B, C, D),
                args.bs_thresh,
                args.min_per_side,
            )
            outcome = base_outcome

            # STEP 2: optional rooted filter for informative outcomes
            if base_outcome in ("AB|CD", "AC|BD", "AD|BC") and args.bucket_mode != "bipartition":
                passes = _rf_check_rooted_constraints(
                    tree=tree,
                    topology=base_outcome,
                    A=A,
                    B=B,
                    C=C,
                    D=D,
                    mode=args.bucket_mode,
                    max_strays=args.max_strays,
                )
                if not passes:
                    outcome = "root_fail"

            counts[qname][outcome] += 1

            # Always record unresolved / insufficient / root_fail
            if outcome in ("unresolved", "insufficient", "root_fail"):
                writers[qname][outcome].write(gid + "\n")
                continue  # never copy for these

            # For informative outcomes: write lists only if requested
            if do_lists:
                writers[qname][outcome].write(gid + "\n")

            # Copy trees/seqs only if requested and we have a source file
            if (do_copy_trees or do_copy_seqs) and src_path is not None:
                dest_dir = buckets[qname][outcome]

                if do_copy_trees:
                    shutil.copy2(src_path, dest_dir / src_path.name)

                if do_copy_seqs and seq_dir is not None:
                    cand = _qb_find_matching_seq_by_suffix(
                        tree_name=src_path.name,
                        seq_dir=seq_dir,
                        seq_ext=args.seq_ext,
                        strip_suffixes=getattr(args, "seq_strip_suffix", []),
                    )
                    if cand is not None:
                        try:
                            shutil.copy2(cand, dest_dir / cand.name)
                        except Exception as e:
                            print(f"[warn] copy seq failed: {cand} → {dest_dir} ({e})",
                                  file=sys.stderr)
                    else:
                        print(
                            f"[info] no matching seq for {src_path.name} "
                            f"(checked {args.seq_strip_suffix or 'default stripping'})",
                            file=sys.stderr,
                        )

    # ---- close and summarize ----
    if do_lists:
        for qname in writers:
            for f in writers[qname].values():
                try:
                    f.close()
                except Exception:
                    pass

    with open(outdir / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "quartet",
            "AB|CD",
            "AC|BD",
            "AD|BC",
            "unresolved",
            "insufficient",
            "root_fail",
            "total",
        ])
        for qname, c in counts.items():
            w.writerow([
                qname,
                c["AB|CD"],
                c["AC|BD"],
                c["AD|BC"],
                c["unresolved"],
                c["insufficient"],
                c["root_fail"],
                sum(c.values()),
            ])

    print(f"[bucket-files] wrote results under: {outdir}")


# =========================
# Sanitizer utilities
# =========================



def register(sub) -> None:
    """Register the `bucket-files` subcommand."""
    # bucket-seqs
    pBQ = sub.add_parser(
        'bucket-files',
        help='Bucket gene trees by quartet support under different quartet models.'
    )
    src = pBQ.add_mutually_exclusive_group(required=True)
    src.add_argument("--gene-trees", help="Concatenated Newick file (one or more trees).")
    src.add_argument("--tree-glob", nargs="+", help="Glob pattern(s) for tree files, e.g. trees/*.tre")
    pBQ.add_argument(
        "--bucket-mode",
        choices=["bipartition", "lenient", "strict"],
        default="bipartition",
        help=(
            "How to interpret quartet support:\n"
            "  bipartition: unrooted, split-based (MSC-style; current default behavior).\n"
            "  lenient: first classify by bipartition, then require some rooted constraints\n"
            "          (e.g. focal clades roughly monophyletic) for an informative label.\n"
            "  strict: as lenient, but with stronger rooted clade requirements.\n"
            "All rooted modes demote failures to a separate 'root_fail' bucket."
        )
    )
    pBQ.add_argument(
    "--max-strays",
    type=int,
    default=2,
    help=("For rooted 'strict' mode: maximum number of focal-taxon strays allowed "
          "outside the main A/B/C clades before the gene tree is demoted to root_fail "
          "(default 2). Ignored in 'bipartition' mode.")
    )
    pBQ.add_argument("--quartet", help='Single quartet like NAME:A1,A2;B1;C1,C2;D1')
    pBQ.add_argument("--quartets-tsv", help="TSV with columns: name, A_set, B_set, C_set, D_set")
    pBQ.add_argument("--bs-thresh", type=float, default=0.0,
                     help="Collapse edges with support < threshold before evaluating (default 0).")
    pBQ.add_argument("--min-per-side", type=int, default=1,
                     help="Minimum taxa required in each of A,B,C,D (default 1).")
    pBQ.add_argument("--outdir", required=True, help="Output root directory.")
    pBQ.add_argument(
        "--out-mode",
        nargs="+",
        choices=["lists", "copy_trees", "copy_seqs"],
        default=["lists"],
        help=("One or more output actions: "
              "'lists' (write per-topology ID lists), "
              "'copy_trees' (copy tree files into bucket dirs), "
              "'copy_seqs' (copy sequence files by basename). "
              "Example: --out-mode lists copy_seqs")
    )
    # Optional copy of FASTAs by basename (only meaningful with --tree-glob)
    pBQ.add_argument("--seq-dir", help="Directory with matching sequence files (for --out-mode copy).")
    pBQ.add_argument("--seq-ext", default=".fas", help="Sequence filename extension (default .fas).")
    pBQ.add_argument(
        "--seq-strip-suffix",
        action="append",
        default=[],
        help="Optional suffix to strip from tree filenames when searching for matching sequence files. "
             "Can be given multiple times. Example: --seq-strip-suffix .fa.tre.rtd --seq-strip-suffix .tre"
    )
    pBQ.set_defaults(func=cmd_bucket_files)

