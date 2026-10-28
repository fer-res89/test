"""Simulate gene trees under the multispecies coalescent on the ASTRAL species
tree and write them as MCMCTree-style outputs (<gene>.FigTree.tre +
<gene>.mcmc.txt) for a data run and a prior-only run.

Two uses:
  1. Validation of `aeon ils-extract` / `aeon ils-analyze`: with --post-sd 0
     --gene-sd 0 --shrink 1, the observed F(c)/sum(F) must match the MSC
     topology probabilities and the observed Delta t_in / Delta t_out (Ma) must
     match the closed-form predictions x 2Ne*g within bootstrap error.
  2. Power / bias analysis: add realistic posterior noise (--post-sd),
     gene-level dating error (--gene-sd), missing taxa (--occupancy) and
     shrinkage of posterior ages toward a topology-dependent prior (--shrink),
     then ask whether the planned statistics still recover the truth.

Sequence mode (--sequences) is the stronger check: it evolves amino-acid
alignments along the simulated gene trees (WAG + Gamma4, independent
lognormal branch rates, gene-to-gene rate variation) and writes

  seqs/<gene>.phy     alignments to date with `aeon pipe-all --root-cal ...`
  trees/<gene>.tre    the true rooted topologies (dates on the true topology,
                      so dating error is separated from gene-tree error)
  truth/<gene>/       noise-free FigTree.tre + mcmc.txt with the TRUE ages

Run `aeon ils-extract` + `aeon ils-analyze` on the MCMCTree output and on
truth/, and compare: any difference is what MCMCTree (rates, priors, root
calibration, finite data) does to the MSC signal.

The "posterior" in the default mode is a caricature: centre = shrink * true age + (1 - shrink)
* prior centre, times multiplicative gene and sample noise, with ages forced
to be ordered. The prior centre depends only on topology (node depth from the
root and longest path to a tip), mimicking how a birth-death prior with a
calibrated root spaces nodes. It is a device for checking the statistics, not
a model of MCMCTree; sequence simulation through the real dating pipeline is
the stronger check.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ._common import WAG_DAT_EMBED
from .ils_common import Node, read_tree_file, renumber, reroot_on_leaf

AA = "ARNDCQEGHILKMFPSTWYV"


class AAModel:
    """WAG + discrete Gamma(4) substitution model, normalised to 1 sub/site."""

    def __init__(self, alpha: float, ncat: int = 4):
        from scipy.stats import gamma as G
        nums = [float(x) for x in WAG_DAT_EMBED.split()]
        exch, pi = nums[:190], np.array(nums[190:210])
        pi = pi / pi.sum()
        S = np.zeros((20, 20)); k = 0
        for i in range(1, 20):
            for j in range(i):
                S[i, j] = S[j, i] = exch[k]; k += 1
        Q = S * pi[None, :]
        np.fill_diagonal(Q, -Q.sum(1))
        Q /= -(pi * np.diag(Q)).sum()
        sq = np.sqrt(pi)
        B = (sq[:, None] * Q) / sq[None, :]
        lam, U = np.linalg.eigh((B + B.T) / 2)
        self.pi, self.lam, self.U, self.sq = pi, lam, U, sq
        # category means of a mean-1 Gamma(alpha) (equal-probability bins)
        edges = G.ppf(np.linspace(0, 1, ncat + 1), alpha, scale=1 / alpha)
        from scipy.special import gammainc
        cdf1 = gammainc(alpha + 1, edges * alpha)
        self.rates = (np.diff(cdf1) * ncat)

    def P(self, t):
        M = (self.U * np.exp(self.lam * t)) @ self.U.T
        P = M * (self.sq[None, :] / self.sq[:, None])
        P = np.clip(P, 0, None)
        return P / P.sum(1, keepdims=True)


def evolve(tree_g, model: "AAModel", nsites: int, mu: float, sigma: float, rng):
    """Evolve along a GNode tree (ages in Myr). Returns {tip name: sequence}."""
    cats = rng.integers(0, len(model.rates), nsites)
    root = rng.choice(20, nsites, p=model.pi)
    out = {}

    def rec(x, states):
        if not x.children:
            out[x.name] = "".join(AA[i] for i in states); return
        for c in x.children:
            r = mu * np.exp(rng.normal(-sigma ** 2 / 2, sigma))
            b = r * (x.age - c.age)
            new = np.empty_like(states)
            for k, rk in enumerate(model.rates):
                m = cats == k
                if not m.any():
                    continue
                cum = np.cumsum(model.P(b * rk)[states[m]], 1)
                new[m] = (rng.random(m.sum())[:, None] > cum).sum(1).clip(0, 19)
            rec(c, new)
    rec(tree_g, root)
    return out


def species_heights(sp: Node, tip_cu: float, root_branch_cu: float, outgroup: str):
    """Heights (CU above present) for an outgroup-rooted species tree that keep
    every internal branch at exactly its ASTRAL length: internal nodes are
    placed top-down from the crown, and terminal branches absorb the rest
    (their length does not matter with one sampled lineage per species)."""
    crown = [c for c in sp.children if not (c.is_leaf() and c.name == outgroup)][0]
    depth = {crown: 0.0}
    for nd in crown.preorder():
        for c in nd.children:
            if not c.is_leaf():
                if c.length is None or not np.isfinite(c.length):
                    raise SystemExit("species tree: internal branch without CU length")
                depth[c] = depth[nd] + c.length
    H = max(depth.values()) + tip_cu
    h = {sp: H + root_branch_cu}
    for nd in sp.preorder():
        if nd is sp:
            continue
        h[nd] = 0.0 if nd.is_leaf() else H - depth[nd]
    return h


class GNode:
    __slots__ = ("children", "age", "name")

    def __init__(self, age, name=None, children=None):
        self.age = age; self.name = name; self.children = children or []


def sim_gene_tree(sp: Node, h, present: set, rng):
    lineages = {}
    for nd in sp.postorder():
        if nd.is_leaf():
            lineages[nd] = [GNode(0.0, nd.name)] if nd.name in present else []
            continue
        pool = []
        for c in nd.children:
            L = list(lineages[c])
            t, t_end = h[c], h[nd]
            while len(L) > 1:
                k = len(L)
                t += rng.exponential(2.0 / (k * (k - 1)))
                if t >= t_end:
                    break
                i, j = rng.choice(k, 2, replace=False)
                new = GNode(t, children=[L[i], L[j]])
                L = [x for m, x in enumerate(L) if m not in (i, j)] + [new]
            pool.extend(L)
        lineages[nd] = pool
    L = lineages[sp]
    t = h[sp]
    while len(L) > 1:
        k = len(L)
        t += rng.exponential(2.0 / (k * (k - 1)))
        i, j = rng.choice(k, 2, replace=False)
        new = GNode(t, children=[L[i], L[j]])
        L = [x for m, x in enumerate(L) if m not in (i, j)] + [new]
    return L[0]


def to_node(g: GNode) -> Node:
    def rec(x, parent_age):
        nd = Node()
        nd.name = x.name
        nd.length = parent_age - x.age if parent_age is not None else None
        for c in x.children:
            ch = rec(c, x.age); ch.parent = nd; nd.children.append(ch)
        return nd
    r = rec(g, None)
    renumber(r)
    return r


def root_on_outgroup(g: GNode, outgroup: str) -> bool:
    return any(c.name == outgroup for c in g.children)


def prior_centres(tree: Node, root_age: float):
    depth = {tree: 0}
    for nd in tree.preorder():
        for c in nd.children:
            depth[c] = depth[nd] + 1
    m = {}
    for nd in tree.postorder():
        m[nd] = 0 if nd.is_leaf() else 1 + max(m[c] for c in nd.children)
    internals = sorted([n for n in tree.preorder() if not n.is_leaf()], key=lambda n: n.idx)
    return np.array([root_age * m[n] / (depth[n] + m[n]) for n in internals]), internals


def ordered_samples(centres, internals, sd, nsamp, rng, root_sd=None):
    k = len(centres)
    x = centres[None, :] * np.exp(rng.normal(0, sd, (nsamp, k))) if sd > 0 else np.repeat(centres[None, :], nsamp, 0)
    if root_sd is not None:
        x[:, 0] = centres[0] * np.exp(rng.normal(0, root_sd, nsamp))
    # enforce parent > child, from tips up
    for nd in sorted(internals, key=lambda n: -n.idx):
        for c in nd.children:
            if not c.is_leaf():
                x[:, nd.idx] = np.maximum(x[:, nd.idx], x[:, c.idx] * 1.0001)
    return x


def write_outputs(tree: Node, internals, samples, path_stem: Path):
    ns = len(tree.leaves())
    mean = samples.mean(0)
    lo, hi = np.quantile(samples, [0.025, 0.975], axis=0)
    age = {n: mean[n.idx] for n in internals}

    def rec(nd):
        if nd.is_leaf():
            s = nd.name; a = 0.0
        else:
            s = "(" + ", ".join(rec(c) for c in nd.children) + ")"
            s += f" [&95%HPD={{{lo[nd.idx]:.6f}, {hi[nd.idx]:.6f}}}]"
            a = age[nd]
        if nd.parent is not None:
            s += f": {age[nd.parent] - a:.6f}"
        return s
    nwk = rec(tree)
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    Path(str(path_stem) + ".FigTree.tre").write_text(f"#NEXUS\nBEGIN TREES;\n\n\tUTREE 1 = {nwk};\n\nEND;\n")
    cols = [f"t_n{ns + 1 + n.idx}" for n in internals]
    with open(str(path_stem) + ".mcmc.txt", "w") as f:
        f.write("Gen\t" + "\t".join(cols) + "\tmu\tsigma2\tlnL\n")
        for i, row in enumerate(samples):
            f.write(f"{(i + 1) * 10}\t" + "\t".join(f"{v:.6f}" for v in row) + "\t0.1\t0.1\t-1000\n")


def add_arguments(ap):
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--outgroup", default="Gallus_gallus")
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--two-ne-g", type=float, default=10.0, help="Myr per coalescent unit (2Ne*g)")
    ap.add_argument("--root-branch-cu", type=float, default=5.0,
                    help="CU from the ingroup root to the tree root (outgroup split)")
    ap.add_argument("--tip-cu", type=float, default=1.0, help="minimum terminal branch length (CU)")
    ap.add_argument("--time-unit-ma", type=float, default=100.0, help="MCMCTree time unit in Myr")
    ap.add_argument("--occupancy", type=float, default=1.0, help="per-taxon presence probability")
    ap.add_argument("--post-sd", type=float, default=0.0, help="log-sd of posterior samples per node")
    ap.add_argument("--gene-sd", type=float, default=0.0, help="log-sd of a gene-wide dating error")
    ap.add_argument("--prior-sd", type=float, default=0.15, help="log-sd of prior-only samples per node")
    ap.add_argument("--shrink", type=float, default=1.0,
                    help="posterior centre = shrink*true + (1-shrink)*prior centre (1 = no shrinkage)")
    ap.add_argument("--ingroup-root-age-ma", "--crown-age-ma", dest="ingroup_root_age_ma", type=float,
                    help="age (Ma) of the ingroup root, i.e. the crown node of the ingroup "
                         "(e.g. crown Lepidosauria); overrides --tip-cu (terminal branches absorb it)")
    ap.add_argument("--root-age-ma", type=float,
                    help="age (Ma) of the root of the whole tree, the outgroup-ingroup split "
                         "(e.g. Gallus-Lepidosauria = crown Sauria); needs --ingroup-root-age-ma; "
                         "overrides --root-branch-cu")
    ap.add_argument("--nsamp", type=int, default=400)
    g = ap.add_argument_group("sequence mode")
    g.add_argument("--sequences", action="store_true",
                   help="write amino-acid alignments + true topologies + truth outputs instead of fake posteriors")
    g.add_argument("--sites", type=int, default=400, help="alignment length (aa)")
    g.add_argument("--mu", type=float, default=0.002, help="mean substitution rate, subs/site/Myr")
    g.add_argument("--gene-rate-cv", type=float, default=0.4, help="CV of the per-gene mean rate (gamma)")
    g.add_argument("--branch-rate-sd", type=float, default=0.2,
                   help="log-sd of independent lognormal branch rates (clock = 2 model)")
    g.add_argument("--alpha", type=float, default=0.5, help="Gamma shape for among-site rate variation")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-dir", default="sim")


def run(a):

    rng = np.random.default_rng(a.seed)
    sp = reroot_on_leaf(read_tree_file(a.species_tree), a.outgroup)
    if a.ingroup_root_age_ma is not None:
        h0 = species_heights(sp, 0.0, 0.0, a.outgroup)
        crown = [c for c in sp.children if not (c.is_leaf() and c.name == a.outgroup)][0]
        a.tip_cu = a.ingroup_root_age_ma / a.two_ne_g - h0[crown]
        if a.tip_cu <= 0:
            raise SystemExit("--ingroup-root-age-ma too young for these CU branch lengths and --two-ne-g")
        if a.root_age_ma is not None:
            if a.root_age_ma <= a.ingroup_root_age_ma:
                raise SystemExit("--root-age-ma must be older than --ingroup-root-age-ma")
            a.root_branch_cu = (a.root_age_ma - a.ingroup_root_age_ma) / a.two_ne_g
    elif a.root_age_ma is not None:
        raise SystemExit("--root-age-ma needs --ingroup-root-age-ma")
    h = species_heights(sp, a.tip_cu, a.root_branch_cu, a.outgroup)
    taxa = sorted(sp.leaf_names())
    ing = [t for t in taxa if t != a.outgroup]
    root_ma = h[sp] * a.two_ne_g
    print(f"[INFO] tree root (outgroup split) at {root_ma:.1f} Ma; ingroup root "
          f"{(h[sp] - a.root_branch_cu) * a.two_ne_g:.1f} Ma; "
          f"root branch {a.root_branch_cu:.2f} CU")
    # species-tree node ages in Ma for the absolute test
    def nwk(nd):
        if nd.is_leaf():
            s = nd.name
        else:
            s = "(" + ",".join(nwk(c) for c in nd.children) + ")"
        if nd.parent is not None:
            s += f":{(h[nd.parent] - h[nd]) * a.two_ne_g:.6f}"
        return s
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    (out / "species_tree_true_Ma.tre").write_text(nwk(sp) + ";\n")
    model = AAModel(a.alpha) if a.sequences else None
    if a.sequences:
        for d in ("seqs", "trees"):
            (out / d).mkdir(exist_ok=True)
        ru = root_ma / a.time_unit_ma
        cal = f"B({0.95 * ru:.4f},{1.05 * ru:.4f},0.025,0.025)"
        (out / "root_calibration.txt").write_text(
            f"{cal}\n# tree root (outgroup-ingroup species split) = {root_ma:.2f} Ma (+/- 5%). Gene-tree roots are older\n"
            f"# by the coalescence time in the ancestral population (mean 1 CU = {a.two_ne_g:g} Myr),\n"
            f"# exactly as in the empirical design.\n")
        print(f"[INFO] suggested root calibration: --root-cal \"{cal}\"")
    for g in range(a.n_genes):
        gid = f"SIM{g:05d}"
        while True:
            present = {t for t in ing if rng.random() < a.occupancy} | {a.outgroup}
            if len(present) >= 4:
                break
        gt = sim_gene_tree(sp, h, present, rng)
        if not root_on_outgroup(gt, a.outgroup):
            # Gallus coalesces last with probability ~1 for long root branches;
            # if not, MCMCTree would date a different root: skip like a bad locus
            continue
        tree = to_node(gt)
        internals = sorted([n for n in tree.preorder() if not n.is_leaf()], key=lambda n: n.idx)
        unit = a.two_ne_g / a.time_unit_ma
        true_age = {}
        def fill(x, nd):
            true_age[nd.idx] = x.age * unit
            for xc, nc in zip(x.children, nd.children):
                if not nc.is_leaf():
                    fill(xc, nc)
        fill(gt, tree)
        true = np.array([true_age[n.idx] for n in internals])
        if a.sequences:
            unit_myr = a.two_ne_g  # GNode ages are in CU
            def to_myr(x):
                x.age *= unit_myr
                for c in x.children:
                    to_myr(c)
            to_myr(gt)
            shape = 1 / a.gene_rate_cv ** 2 if a.gene_rate_cv > 0 else None
            mu_g = rng.gamma(shape, a.mu / shape) if shape else a.mu
            seqs = evolve(gt, model, a.sites, mu_g, a.branch_rate_sd, rng)
            with open(out / "seqs" / f"{gid}.phy", "w") as f:
                f.write(f"{len(seqs)} {a.sites}\n")
                for name in sorted(seqs):
                    f.write(f"{name}  {seqs[name]}\n")
            def topo(nd):
                return nd.name if nd.is_leaf() else "(" + ",".join(topo(c) for c in nd.children) + ")"
            (out / "trees" / f"{gid}.tre").write_text(topo(tree) + ";\n")
            write_outputs(tree, internals, np.repeat(true[None, :], 10, 0), out / "truth" / gid / gid)
            continue
        pc, _ = prior_centres(tree, root_ma / a.time_unit_ma)
        centre = a.shrink * true + (1 - a.shrink) * pc
        if a.gene_sd > 0:
            centre[1:] *= np.exp(rng.normal(0, a.gene_sd))
        ds = ordered_samples(centre, internals, a.post_sd, a.nsamp, rng)
        ps = ordered_samples(pc, internals, a.prior_sd, a.nsamp, rng, root_sd=0.02)
        write_outputs(tree, internals, ds, out / "data" / gid / gid)
        write_outputs(tree, internals, ps, out / "prior" / gid / gid)
    print(f"[DONE] wrote {out}/" + ("seqs, trees, truth" if a.sequences else "data and prior"))



def register(sub) -> None:
    p = sub.add_parser('ils-simulate', help='Simulate MSC gene trees as MCMCTree output (validation/power)', description=__doc__,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(p)
    p.set_defaults(func=run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return run(ap.parse_args(argv))
