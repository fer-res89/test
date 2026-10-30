"""Shared code for the gene-tree dating ILS test (gene_tree_dating_ILS_test.md).

Contents
  * a small dependency-free Newick/NEXUS parser that keeps labels, comments
    ([&95%HPD=...], ASTRAL q-annotations) and literal underscores
  * rerooting on an outgroup leaf
  * species-tree branch definitions (Section 9.1) and the polytomy node (Section 13)
  * closed-form MSC expectations (Section 3.3 / Appendix A)
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Newick parsing
# --------------------------------------------------------------------------


class Node:
    __slots__ = ("children", "parent", "name", "length", "comment", "label", "idx")

    def __init__(self):
        self.children: list[Node] = []
        self.parent: Node | None = None
        self.name: str | None = None      # tip name or internal label
        self.length: float | None = None
        self.comment: str = ""            # concatenated [...] comments
        self.label: str = ""              # raw internal label (quoted or not)
        self.idx: int = -1                # preorder index of internal nodes

    def is_leaf(self) -> bool:
        return not self.children

    def preorder(self):
        stack = [self]
        while stack:
            nd = stack.pop()
            yield nd
            stack.extend(reversed(nd.children))

    def postorder(self):
        out = list(self.preorder())
        return reversed(out)

    def leaves(self):
        return [n for n in self.preorder() if n.is_leaf()]

    def leaf_names(self) -> frozenset[str]:
        return frozenset(n.name for n in self.leaves())


def parse_newick(s: str) -> Node:
    """Parse one Newick string. Internal nodes get .idx in the order their
    opening parenthesis appears (= preorder = MCMCTree/PAML node numbering
    offset). Quoted labels and [comments] are preserved verbatim."""
    s = s.strip()
    j = _end_of_tree(s)
    if j >= 0:
        s = s[:j]
    i, n = 0, len(s)
    root = None
    cur: Node | None = None
    internal_counter = 0
    stack: list[Node] = []

    def read_comment(j):
        depth, k = 0, j
        while k < n:
            if s[k] == "[":
                depth += 1
            elif s[k] == "]":
                depth -= 1
                if depth == 0:
                    return s[j + 1:k], k + 1
            k += 1
        raise ValueError("unterminated comment")

    def read_label(j):
        if j < n and s[j] == "'":
            k = j + 1
            buf = []
            while k < n:
                if s[k] == "'":
                    if k + 1 < n and s[k + 1] == "'":
                        buf.append("'"); k += 2; continue
                    return "".join(buf), k + 1, True
                buf.append(s[k]); k += 1
            raise ValueError("unterminated quote")
        k = j
        while k < n and s[k] not in ",():;[" and not s[k].isspace():
            k += 1
        return s[j:k], k, False

    def skip_ws(j):
        while j < n and s[j].isspace():
            j += 1
        return j

    def read_tail(nd: Node, j):
        """label, comments, :length, comments after a node."""
        j = skip_ws(j)
        while j < n and s[j] == "[":
            c, j = read_comment(j); nd.comment += c; j = skip_ws(j)
        if j < n and s[j] not in ",():;[":
            lab, j, _ = read_label(j)
            if nd.is_leaf():
                nd.name = lab
            else:
                nd.label = lab
            j = skip_ws(j)
        while j < n and s[j] == "[":
            c, j = read_comment(j); nd.comment += c; j = skip_ws(j)
        if j < n and s[j] == ":":
            j = skip_ws(j + 1)
            m = re.match(r"[-+0-9.eE]+", s[j:])
            if m:
                nd.length = float(m.group(0)); j += m.end()
            j = skip_ws(j)
        while j < n and s[j] == "[":
            c, j = read_comment(j); nd.comment += c; j = skip_ws(j)
        return j

    i = skip_ws(i)
    while i < n:
        ch = s[i]
        if ch == "(":
            nd = Node(); nd.idx = internal_counter; internal_counter += 1
            if stack:
                nd.parent = stack[-1]; stack[-1].children.append(nd)
            else:
                root = nd
            stack.append(nd); i += 1
            # a leaf may follow
            i = skip_ws(i)
            if i < n and s[i] not in "(),":
                lf = Node(); lf.parent = nd; nd.children.append(lf)
                lab, i, _ = read_label(i); lf.name = lab
                i = read_tail(lf, i)
        elif ch == ",":
            i = skip_ws(i + 1)
            if i < n and s[i] not in "(),":
                lf = Node(); lf.parent = stack[-1]; stack[-1].children.append(lf)
                lab, i, _ = read_label(i); lf.name = lab
                i = read_tail(lf, i)
        elif ch == ")":
            nd = stack.pop()
            i = read_tail(nd, i + 1)
        elif ch.isspace():
            i += 1
        else:
            # single-leaf tree or junk
            if root is None:
                root = Node(); lab, i, _ = read_label(i); root.name = lab
                i = read_tail(root, i)
            else:
                raise ValueError(f"unexpected character {ch!r} at {i}")
    if stack:
        raise ValueError("unbalanced parentheses")
    return root


def read_tree_file(path: str | Path) -> Node:
    """Read the first tree from a Newick file, an MCMCTree tree file (with
    '<ntax> <ntree>' header) or a simple NEXUS file (MCMCTree/BPP FigTree.tre;
    no Translate block)."""
    txt = Path(path).read_text()
    if txt.lstrip().upper().startswith("#NEXUS"):
        if re.search(r"(?i)\btranslate\b", txt):
            raise ValueError(f"{path}: NEXUS Translate blocks are not supported")
        m = re.search(r"(?is)\b[UR]?TREE\s+[^=;]*=\s*(?:\[&[RU]\]\s*)?(.*)", txt)
        if not m:
            raise ValueError(f"{path}: no tree found in NEXUS file")
        body = m.group(1)
        return parse_newick(body[: _end_of_tree(body) + 1])
    lines = [ln for ln in txt.splitlines() if ln.strip()]
    if lines and re.fullmatch(r"\s*\d+\s+\d+\s*", lines[0]):
        lines = lines[1:]
    body = " ".join(lines)
    j = _end_of_tree(body)
    return parse_newick(body[: j + 1] if j >= 0 else body)


def _end_of_tree(s: str) -> int:
    """Index of the ';' that ends the first tree, ignoring quotes and [comments]."""
    q = False; depth = 0
    for i, ch in enumerate(s):
        if q:
            if ch == "'":
                q = False
        elif ch == "'":
            q = True
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        elif ch == ";" and depth == 0:
            return i
    return -1


HPD_RE = re.compile(r"95%(?:HPD)?\s*=\s*\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\}")


def node_ages(root: Node) -> dict[Node, float]:
    """Ages (height above the tips) for an ultrametric tree with lengths.
    Uses the maximum root-to-tip depth so rounding does not produce negative
    tip ages."""
    depth = {root: 0.0}
    for nd in root.preorder():
        for c in nd.children:
            depth[c] = depth[nd] + (c.length or 0.0)
    h = max(depth[lf] for lf in root.leaves())
    return {nd: h - d for nd, d in depth.items()}


def to_newick(root: Node, lengths=True, internal_labels=False) -> str:
    def rec(nd):
        if nd.is_leaf():
            s = nd.name
        else:
            s = "(" + ",".join(rec(c) for c in nd.children) + ")"
            if internal_labels and nd.label:
                s += nd.label
        if lengths and nd.length is not None and nd is not root:
            s += f":{nd.length:.8g}"
        return s
    return rec(root) + ";"


# --------------------------------------------------------------------------
# Tree manipulation
# --------------------------------------------------------------------------

def _reroot_on_node(root: Node, og: Node) -> Node:
    """Return a new rooted tree with `outgroup` as one child of the root.
    Edge lengths and internal labels/comments stay attached to the node below
    each edge, as far as that is well defined (labels on edges whose direction
    flips are moved with the edge)."""
    if og.parent is root and len(root.children) == 2:
        return root
    # Build undirected adjacency with edge data keyed by the child node
    adj: dict[Node, list[tuple[Node, Node]]] = {}
    for nd in root.preorder():
        adj.setdefault(nd, [])
        for c in nd.children:
            adj[nd].append((c, c))      # (neighbor, edge-owner)
            adj.setdefault(c, []).append((nd, c))
    new_root = Node()
    other = og.parent
    edge_len = og.length or 0.0

    def build(nd: Node, frm: Node | None, edge_owner: Node | None, length):
        cp = Node()
        cp.name = nd.name
        cp.length = length
        if edge_owner is not None:
            cp.label = edge_owner.label
            cp.comment = edge_owner.comment
        for nb, owner in adj[nd]:
            if nb is frm:
                continue
            ch = build(nb, nd, owner, owner.length)
            ch.parent = cp
            cp.children.append(ch)
        # suppress unary nodes (old root of degree 2)
        if len(cp.children) == 1 and not nd.is_leaf():
            only = cp.children[0]
            only.length = (only.length or 0.0) + (cp.length or 0.0)
            only.parent = None
            return only
        return cp

    a = build(og, other, og, edge_len)
    if og.is_leaf():
        a.children = []
        a.name = og.name
    b = build(other, og, og, 0.0)
    for x in (a, b):
        x.parent = new_root
        new_root.children.append(x)
    a.label = a.comment = ""
    renumber(new_root)
    return new_root


def reroot_on_leaf(root: Node, outgroup: str) -> Node:
    """Root on the edge above one leaf."""
    leaves = {lf.name: lf for lf in root.leaves()}
    if outgroup not in leaves:
        raise ValueError(f"outgroup {outgroup} not in tree")
    return _reroot_on_node(root, leaves[outgroup])


def parse_outgroup(spec) -> frozenset:
    """'Gallus_gallus' or 'Mus_musculus,Homo_sapiens,...' -> set of taxa."""
    if isinstance(spec, (set, frozenset, list, tuple)):
        return frozenset(spec)
    return frozenset(t.strip() for t in str(spec).split(",") if t.strip())


class OutgroupError(ValueError):
    pass


def reroot_on_outgroup(root: Node, outgroup) -> Node:
    """Root on the edge separating the outgroup taxa (those present in the
    tree) from everything else. The outgroup may be one leaf or a clade; if
    the present outgroup taxa are not a clade on either side of some edge,
    OutgroupError is raised."""
    og = parse_outgroup(outgroup) & root.leaf_names()
    if not og:
        raise OutgroupError("no outgroup taxon in tree")
    if len(og) == 1:
        return reroot_on_leaf(root, next(iter(og)))
    allt = root.leaf_names()
    if len(root.children) == 2 and any(c.leaf_names() == og for c in root.children):
        return root
    rest = allt - og
    for nd in root.preorder():
        if nd is root:
            continue
        ls = nd.leaf_names()
        if ls == og:
            return _reroot_on_node(root, nd)
        if ls == rest:                     # outgroup is on the root side of this edge
            return _reroot_on_node(root, nd)
    raise OutgroupError("outgroup taxa are not monophyletic")


def outgroup_at_root(root: Node, outgroup) -> bool:
    og = parse_outgroup(outgroup) & root.leaf_names()
    return bool(og) and len(root.children) == 2 and any(c.leaf_names() == og for c in root.children)


def renumber(root: Node):
    k = 0
    for nd in root.preorder():
        if not nd.is_leaf():
            nd.idx = k; k += 1


def prune_to(root: Node, keep: set[str]) -> Node | None:
    """Prune leaves not in keep; suppress unary nodes. Returns new root."""
    def rec(nd):
        if nd.is_leaf():
            return nd if nd.name in keep else None
        kids = [k for k in (rec(c) for c in nd.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            k = kids[0]
            k.length = (k.length or 0.0) + (nd.length or 0.0)
            return k
        nd.children = kids
        for k in kids:
            k.parent = nd
        return nd
    r = rec(root)
    if r is not None:
        r.parent = None
        r.length = None
        renumber(r)
    return r


# --------------------------------------------------------------------------
# MSC expectations (Section 3.3, Appendix A)
# --------------------------------------------------------------------------

def msc_expectations(T: float) -> dict[str, float]:
    """Closed forms for a branch of length T coalescent units."""
    if T is None or not np.isfinite(T):
        return dict(P_C=np.nan, P_D_each=np.nan, dTin_CU=np.nan, dTout_CU=np.nan,
                    frac_C_older=np.nan)
    e = math.exp(-T)
    PC = 1 - (2 / 3) * e
    w1 = 1 - e
    w2 = e / 3
    if T > 0:
        Es = 1 - T * e / (1 - e)
    else:
        Es = 0.0
    E_in_C = (w1 * (Es - T) + w2 / 3) / PC
    E_out_C = (w1 * 1 + w2 * 4 / 3) / PC
    return dict(P_C=PC, P_D_each=e / 3,
                dTin_CU=1 / 3 - E_in_C,
                dTout_CU=4 / 3 - E_out_C,
                frac_C_older=w2 / PC)


# --------------------------------------------------------------------------
# Species-tree branches (Section 9.1) and the polytomy node (Section 13)
# --------------------------------------------------------------------------

@dataclass
class Branch:
    bid: str
    kind: str                       # "binary" or "polytomy"
    groups: list[list[str]]         # binary: [K1, K2, S]; polytomy: [G1, G2, G3]
    group_names: list[str]
    T: float = float("nan")
    classes: list[str] = field(default_factory=list)
    q_astral: dict = field(default_factory=dict)   # class -> q (if mappable)
    q_note: str = ""

    @property
    def n_triplets(self):
        a, b, c = (len(g) for g in self.groups)
        return a * b * c

    def label(self):
        return " | ".join(self.group_names)


def _clade_name(taxa: list[str], clade_names: dict[frozenset, str]):
    fs = frozenset(taxa)
    if fs in clade_names:
        return clade_names[fs]
    if len(taxa) == 1:
        return taxa[0]
    return f"{sorted(taxa)[0]}+{len(taxa) - 1}"


def parse_q(nd: Node) -> dict[str, float]:
    txt = (nd.label or "") + " " + (nd.comment or "")
    out = {}
    for k in ("q1", "q2", "q3"):
        m = re.search(rf"\b{k}\s*=\s*([0-9.eE+-]+)", txt)
        if m:
            out[k] = float(m.group(1))
    return out


def _astral_q_for_branch(orig_root: Node, K1, K2, S, O) -> tuple[dict, str]:
    """Map ASTRAL q1/q2/q3 (as annotated on the species tree in the file's own
    orientation) onto C / D1 / D2.

    Convention assumed (ASTRAL-III -t 2): at node v with children L1, L2,
    sibling R1 and remainder R2: q1 = L1L2|R1R2, q2 = L1R1|L2R2,
    q3 = L1R2|L2R1. The four groups are matched by taxon sets, so the mapping
    does not depend on how the tree was rerooted afterwards -- but it is only
    valid if the annotations were written in this orientation.
    Here: C = K1K2|SO, D1 = K1S|K2O, D2 = K2S|K1O.
    """
    sets = {"K1": frozenset(K1), "K2": frozenset(K2), "S": frozenset(S), "O": frozenset(O)}
    allt = orig_root.leaf_names()
    for v in orig_root.preorder():
        if v.is_leaf() or v.parent is None:
            continue
        q = parse_q(v)
        if len(q) < 3 or len(v.children) != 2:
            continue
        L1, L2 = (c.leaf_names() for c in v.children)
        sibs = [c for c in v.parent.children if c is not v]
        if len(sibs) != 1:
            continue  # sibling ambiguous at a multifurcating parent
        R1 = sibs[0].leaf_names()
        R2 = allt - L1 - L2 - R1
        quad = {"L1": L1, "L2": L2, "R1": R1, "R2": R2}
        if set(quad.values()) != set(sets.values()):
            continue

        def pairing(a, b):
            """frozenset of two frozensets: the unrooted split {a∪b | rest}."""
            return frozenset([a | b, allt - (a | b)])
        qtop = {"q1": pairing(L1, L2), "q2": pairing(L1, R1), "q3": pairing(L1, R2)}
        ours = {"C": pairing(sets["K1"], sets["K2"]),
                "D1": pairing(sets["K1"], sets["S"]),
                "D2": pairing(sets["K2"], sets["S"])}
        out = {}
        for cls, split in ours.items():
            for qk, qs in qtop.items():
                if qs == split:
                    out[cls] = q[qk]
        if len(out) == 3:
            return out, "mapped"
    return {}, "no matching annotated node (or sibling ambiguous)"


def define_branches(species_tree: Node, outgroup: str, ingroup: set[str],
                    polytomy: list[tuple[str, list[str]]] | None = None,
                    clade_names: dict[frozenset, str] | None = None,
                    polytomy_T: float = float("nan")) -> tuple[list[Branch], list[str]]:
    """Branch-specific triplet groups for every testable internal branch.

    species_tree: as read from file (ASTRAL output). It is rerooted on the
    outgroup; q-annotations are read in the original orientation.
    Returns (branches, warnings)."""
    clade_names = clade_names or {}
    warnings = []
    orig = species_tree
    rt = reroot_on_outgroup(orig, outgroup)
    allt = rt.leaf_names()
    poly_sets = None
    if polytomy:
        poly_sets = {frozenset(t) for _, t in polytomy}
    branches = []
    k = 0
    for K in rt.preorder():
        if K.is_leaf() or K.parent is None:
            continue
        N = K.parent
        if N.parent is None:
            continue  # N is the root: S would be the outgroup
        if len(K.children) != 2 or len(N.children) != 2:
            warnings.append(f"multifurcation at clade {_clade_name(sorted(K.leaf_names()), clade_names)}: "
                            f"branch skipped (use --polytomy)")
            continue
        K1 = sorted(K.children[0].leaf_names())
        K2 = sorted(K.children[1].leaf_names())
        S = sorted([c for c in N.children if c is not K][0].leaf_names())
        if not (set(K1 + K2 + S) <= ingroup):
            continue
        if poly_sets and {frozenset(K1), frozenset(K2), frozenset(S)} == poly_sets:
            warnings.append("species-tree branch matching the --polytomy groups replaced by the 3-way node")
            continue
        O = allt - set(K1) - set(K2) - set(S)
        T = K.length if K.length is not None else float("nan")
        k += 1
        br = Branch(bid=f"b{k:02d}", kind="binary", groups=[K1, K2, S],
                    group_names=[_clade_name(K1, clade_names), _clade_name(K2, clade_names),
                                 _clade_name(S, clade_names)],
                    T=T, classes=["C", "D1", "D2"])
        br.q_astral, br.q_note = _astral_q_for_branch(orig, K1, K2, S, O)
        branches.append(br)
    if polytomy:
        if len(polytomy) != 3:
            raise ValueError("--polytomy needs exactly three groups")
        names = [g for g, _ in polytomy]
        groups = [sorted(t) for _, t in polytomy]
        for g, t in polytomy:
            miss = set(t) - allt
            if miss:
                raise ValueError(f"polytomy group {g}: taxa not in species tree: {sorted(miss)}")
        cls = [f"{names[0]}+{names[1]}", f"{names[0]}+{names[2]}", f"{names[1]}+{names[2]}"]
        branches.append(Branch(bid="poly", kind="polytomy", groups=groups, group_names=names,
                               T=polytomy_T, classes=cls))
    return branches, warnings


def species_node_for_branch(br: Branch):
    """Taxon sets whose MRCAs are the branch's species nodes:
    tau_in = MRCA(K1 ∪ K2), tau_out = MRCA(K1 ∪ K2 ∪ S)."""
    a, b, c = br.groups
    return sorted(a + b), sorted(a + b + c)


def mrca_age(root: Node, ages: dict, taxa) -> float:
    taxa = set(taxa)
    best = None
    for nd in root.preorder():
        if taxa <= nd.leaf_names():
            best = nd
    if best is None:
        return float("nan")
    return ages[best]


def read_clade_names(path: str | None) -> dict[frozenset, str]:
    out = {}
    if not path:
        return out
    for ln in Path(path).read_text().splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        name, taxa = ln.split("\t")[:2]
        out[frozenset(t.strip() for t in taxa.split(",") if t.strip())] = name.strip()
    return out


def read_polytomy(path: str | None):
    if not path:
        return None
    rows = []
    for ln in Path(path).read_text().splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        name, taxa = ln.split("\t")[:2]
        rows.append((name.strip(), [t.strip() for t in taxa.split(",") if t.strip()]))
    return rows
