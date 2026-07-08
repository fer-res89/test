"""The `rearrange-triplet` subcommand: from a 3-group taxa list, produce
(g1,g3)|g2 and (g2,g3)|g1 Newick topologies while preserving the original
root split/order.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from ete3 import Tree as EteTree

from .sanitize_tree import parse_sanitize_letters, sanitize_tree_file

BLEN_RE = re.compile(r':-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?')

def _strip_blens(s: str) -> str:
    return BLEN_RE.sub('', s)

def _get_root_sides(newick: str):
    t = EteTree(newick, format=0)
    r = t.get_tree_root()
    if len(r.children) < 2:
        s = {lf.name for lf in t.get_leaves()}
        return s, set()
    left, right = r.children[0], r.children[1]
    left_set  = {lf.name for lf in left.get_leaves()}
    right_set = {lf.name for lf in right.get_leaves()}
    return left_set, right_set

def _reroot_to_original(newick: str, left_set: set, right_set: set) -> str:
    t = EteTree(newick, format=0)
    present_right = [n for n in right_set if t.search_nodes(name=n)]
    if present_right:
        mrca_right = t.get_common_ancestor(present_right)
        if mrca_right is not t.get_tree_root():
            t.set_outgroup(mrca_right)
    r = t.get_tree_root()
    if len(r.children) == 2:
        c0, c1 = r.children
        def has_any(node, names):
            return bool(set(leaf.name for leaf in node.iter_leaves()) & names)
        if not has_any(c0, left_set) and has_any(c1, left_set):
            c0.detach(); c1.detach()
            r.add_child(c1); r.add_child(c0)
    return _strip_blens(t.write(format=1, format_root_node=True))

SINGLETON_RE = re.compile(r'\(([A-Za-z0-9_.-]+)\)')

def _collapse_singletons(newick: str) -> str:
    #Remove parentheses around single taxa, repeat until stable.
    prev = None
    while prev != newick:
        prev = newick
        newick = SINGLETON_RE.sub(r'\1', newick)
    return newick

def _rearrange_triplet(input_newick, sisterA, sisterB, outgroup):
    tree = EteTree(input_newick, format=0)
    setA, setB, setO = set(sisterA), set(sisterB), set(outgroup)
    if (setA & setB) or (setA & setO) or (setB & setO):
        overlap = sorted(list((setA & setB) | (setA & setO) | (setB & setO)))
        raise ValueError("Taxa appear in multiple groups: " + ", ".join(overlap))
    lcaA = tree.get_common_ancestor(sisterA)
    lcaB = tree.get_common_ancestor(sisterB)
    lcaO = tree.get_common_ancestor(outgroup)
    mrca_live = tree.get_common_ancestor(lcaA, lcaB, lcaO)
    sub = mrca_live.copy(method="deepcopy")
    sub_names = {lf.name for lf in sub.get_leaves()}
    A_names = [t for t in sisterA if t in sub_names]
    B_names = [t for t in sisterB if t in sub_names]
    O_names = [t for t in outgroup if t in sub_names]
    if not A_names or not B_names or not O_names:
        missing = []
        if not A_names: missing += [f"A:{x}" for x in sisterA if x not in sub_names]
        if not B_names: missing += [f"B:{x}" for x in sisterB if x not in sub_names]
        if not O_names: missing += [f"O:{x}" for x in outgroup if x not in sub_names]
        raise ValueError("Some taxa not found under MRCA: " + ", ".join(missing))
    union = set(A_names) | set(B_names) | set(O_names)

    def pruned_copy(leafnames):
        c = sub.copy(method="deepcopy")
        c.prune(list(leafnames), preserve_branch_length=True)
        return c

    A_copy = pruned_copy(A_names)
    B_copy = pruned_copy(B_names)
    O_copy = pruned_copy(O_names)

    inner = EteTree(); inner.add_child(A_copy); inner.add_child(B_copy)
    triple = EteTree(); triple.add_child(inner); triple.add_child(O_copy)

    children = list(sub.children)
    def leafset(node): return set(node.get_leaf_names())
    child_sets = [leafset(ch) for ch in children]
    covered_idx = [i for i,S in enumerate(child_sets) if S & union]
    leftovers = []
    for i in covered_idx:
        rem = child_sets[i] - union
        leftovers.append(pruned_copy(rem) if rem else None)

    R = EteTree(); R.add_child(triple)
    for node in leftovers:
        if node is not None:
            R.add_child(node)

    first_idx = covered_idx[0]
    rebuilt = EteTree(); inserted = False
    for i, ch in enumerate(children):
        if i == first_idx and not inserted:
            rebuilt.add_child(R); inserted = True
        if i in covered_idx:
            continue
        rebuilt.add_child(ch.detach())
    if not inserted:
        rebuilt.add_child(R)

    if mrca_live.up is None:
        result_tree = rebuilt
    else:
        up = mrca_live.up
        mrca_live.detach()
        up.add_child(rebuilt)
        result_tree = tree
    raw = result_tree.write(format=1, format_root_node=True)
    return re.sub(BLEN_RE, '', raw)


def cmd_rearrange_triplet(args):
    in_path = Path(args.input)
    san = parse_sanitize_letters(getattr(args, 'sanitize', None))
    if san["do_input"]:
        tmp_in = in_path.with_suffix('.sanitized.input.tmp.tre')
        sanitize_tree_file(
            in_path, tmp_in,
            collapse=san["collapse"],
            keep_support=san["keep_support"],
            keep_lengths=san["keep_lengths"],
            allow_polytomies=san["allow_polytomies"],
            add_header=san["add_header"],
            quiet=san["quiet"]
        )
        in_path = tmp_in
    taxa_file = Path(args.taxa_file)
    if not in_path.is_file():
        sys.exit(f"Error: Tree file '{in_path}' not found.")
    if not taxa_file.is_file():
        sys.exit(f"Error: Taxa file '{taxa_file}' not found.")

    newick = in_path.read_text().strip()
    left_set, right_set = _get_root_sides(newick)

    lines = [line.strip() for line in taxa_file.read_text().splitlines() if line.strip()]
    if len(lines) < 3:
        sys.exit("Error: Taxa file must have at least 3 non-empty lines.")
    group1 = [t.strip() for t in lines[0].split(',') if t.strip()]
    group2 = [t.strip() for t in lines[1].split(',') if t.strip()]
    group3 = [t.strip() for t in lines[2].split(',') if t.strip()]

    outA_newick = _rearrange_triplet(newick, group1, group3, group2)  # (g1,g3)|g2
    outB_newick = _rearrange_triplet(newick, group2, group3, group1)  # (g2,g3)|g1
    outA_newick = _reroot_to_original(outA_newick, left_set, right_set)
    outB_newick = _reroot_to_original(outB_newick, left_set, right_set)
    outA_newick = _collapse_singletons(outA_newick)
    outB_newick = _collapse_singletons(outB_newick)

    if args.output_names:
        names = re.split(r'[,\s]+', args.output_names.strip())
        if len(names) != 2:
            sys.exit("Error: -o/--output-names must provide exactly TWO names.")
        outA_path = Path(names[0] + ".tre")
        outB_path = Path(names[1] + ".tre")
    else:
        base, ext = os.path.splitext(in_path.name)
        ext = ext if ext else ".tre"
        outA_path = Path(f"{base}_g1g3_sis_g2out{ext}")
        outB_path = Path(f"{base}_g2g3_sis_g1out{ext}")

    if san["do_output"]:
        tmpA = outA_path.with_suffix('.sanitized.output.tmp.tre'); tmpA.write_text(outA_newick + "\n")
        tmpB = outB_path.with_suffix('.sanitized.output.tmp.tre'); tmpB.write_text(outB_newick + "\n")
    
        sanitize_tree_file(
            tmpA, outA_path,
            collapse=san["collapse"],
            keep_support=san["keep_support"],
            keep_lengths=san["keep_lengths"],
            allow_polytomies=san["allow_polytomies"],
            add_header=san["add_header"],
            quiet=san["quiet"]
        )
        sanitize_tree_file(
            tmpB, outB_path,
            collapse=san["collapse"],
            keep_support=san["keep_support"],
            keep_lengths=san["keep_lengths"],
            allow_polytomies=san["allow_polytomies"],
            add_header=san["add_header"],
            quiet=san["quiet"]
        )
        for t in (tmpA, tmpB):
            try: t.unlink()
            except: pass
    else:
        outA_path.write_text(outA_newick + "\n")
        outB_path.write_text(outB_newick + "\n")
    
    print(f"Wrote {outA_path}")
    print(f"Wrote {outB_path}")




# =========================
# Sanitize trees
# =========================



def register(sub) -> None:
    """Register the `rearrange-triplet` subcommand."""
    # rearrange-triplet
    p5 = sub.add_parser('rearrange-triplet',
                        help='From one taxa-list (3 lines), output (g1,g3)|g2 and (g2,g3)|g1; preserve root split/order')
    p5.add_argument('-i', '--input', required=True, help='Input Newick tree file')
    p5.add_argument('-t', '--taxa-file', required=True,
                    help='3 lines: group1, group2, group3 (comma-separated taxa)')
    p5.add_argument('-o', '--output-names',
                    help='Two names (comma/space-separated); files written as <name>.tre')
    p5.add_argument(
    '--sanitize',
    nargs='?',           # zero or one argument
    const='',            # `--sanitize` alone => defaults (collapse, drop support, strip lengths)
    help=("Sanitize trees with letter options: "
          "i=input, o=output, l=keep-lengths, s=keep-support, "
          "p=allow-polytomies, h=add-header, c=no-collapse. "
          "Examples: '--sanitize iop', '--sanitize s', '--sanitize'")
    )
    p5.set_defaults(func=cmd_rearrange_triplet)

