"""Step 7 of gene_tree_dating_ILS_test.md: extract triplet ages from dated gene trees.

Reads, for every gene, the MCMCTree posterior summary (<gene>.FigTree.tre) and
samples (<gene>.mcmc.txt) of the data run and of the prior-only run, and writes

  branches.tsv               one row per species-tree branch (Section 9.1) + MSC predictions
  gene_qc.tsv                one row per gene: status, ESS, checks
  gene_nodes.tsv.gz          one row per gene x run x internal node: median, mean, var,
                             95% HPD, ESS, and (data run) 41 posterior quantiles
  gene_branch_pairs.tsv.gz   one row per gene x branch x class x (t_in node, t_out node):
                             number of branch-specific triplets reading that node pair,
                             posterior-median r = t_in/t_out (data and prior)
  triplets_<branch>.tsv.gz   optional (--triplet-branches): gene x triplet rows

Ages are kept in MCMCTree time units; `aeon ils-analyze` converts to Ma.

The node-pair table is the lossless compression of the gene x triplet table: all
branch-specific triplets of a gene that read the same two nodes collapse into one
row with a count, so class weights w(g,c), class mean ages, absolute-test
indicators and posterior probabilities can all be recomputed downstream.

MCMCTree node numbering: tips 1..s, root s+1, internal nodes numbered in preorder
(order of opening parentheses). This mapping is verified for every gene by
checking that each mcmc.txt column mean equals the FigTree.tre posterior-mean
age of the node it is mapped to; genes failing the check are flagged, not used.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import multiprocessing as mp
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .ils_common import (HPD_RE, Branch, define_branches, msc_expectations, node_ages,
                            read_clade_names, read_polytomy, read_tree_file,
                            species_node_for_branch)

QLEVELS = np.linspace(0, 1, 41)
QCOLS = [f"q{int(round(lv * 1000)):04d}" for lv in QLEVELS]
NODE_COLS = ["gene", "run", "node", "n_desc", "median", "mean", "var", "hpd_lo", "hpd_hi", "ess",
             *QCOLS, "clade"]

# --------------------------------------------------------------------------
# ESS (Geyer initial monotone sequence), vectorized over columns
# --------------------------------------------------------------------------


def ess(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    if x.ndim == 1:
        x = x[:, None]
    n, k = x.shape
    if n < 10:
        return np.full(k, np.nan)
    xc = x - x.mean(0)
    var = (xc ** 2).mean(0)
    out = np.full(k, np.nan)
    m = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(xc, n=m, axis=0)
    acov = np.fft.irfft(f * np.conj(f), n=m, axis=0)[:n] / n
    for j in range(k):
        if var[j] <= 0:
            out[j] = np.nan
            continue
        rho = acov[:, j] / acov[0, j]
        # pair sums Gamma_t = rho_2t + rho_2t+1
        npairs = (n - 1) // 2
        g = rho[0:2 * npairs:2] + rho[1:2 * npairs:2]
        pos = np.where(g <= 0)[0]
        cut = pos[0] if len(pos) else len(g)
        g = np.minimum.accumulate(g[:cut]) if cut > 0 else g[:0]
        tau = -1 + 2 * g.sum()
        out[j] = n / max(tau, 1e-12)
    return out


# --------------------------------------------------------------------------
# Reading one MCMCTree run
# --------------------------------------------------------------------------


class RunError(Exception):
    pass


def read_run(fig: Path, mcmc: Path | None, discard_frac: float, check_tol: float):
    tree = read_tree_file(fig)
    tips = tree.leaves()
    ns = len(tips)
    internals = sorted([n for n in tree.preorder() if not n.is_leaf()], key=lambda n: n.idx)
    if any(len(n.children) != 2 for n in internals):
        raise RunError("not_bifurcating")
    ages = node_ages(tree)
    fig_mean = np.array([ages[n] for n in internals])
    hpd = np.full((len(internals), 2), np.nan)
    for i, n in enumerate(internals):
        m = HPD_RE.search(n.comment or "")
        if m:
            hpd[i] = float(m.group(1)), float(m.group(2))
    samples = None
    maxdiff = np.nan
    if mcmc is not None and mcmc.exists():
        cols = [f"t_n{ns + 1 + n.idx}" for n in internals]
        hdr = pd.read_csv(mcmc, sep="\t", nrows=0).columns
        hdr = [c.strip() for c in hdr]
        missing = [c for c in cols if c not in hdr]
        if missing:
            raise RunError(f"mcmc_columns_missing:{missing[:3]}")
        df = pd.read_csv(mcmc, sep="\t", usecols=lambda c: c.strip() in set(cols))
        df.columns = [c.strip() for c in df.columns]
        samples = df[cols].to_numpy(float)
        # verify the node mapping against FigTree.tre posterior means (full chain)
        # MCMCTree's FigTree means agree with the mcmc.txt column means only
        # approximately (to ~1e-3 for data runs, ~1e-2 for short prior runs), so
        # the test is: each column is closest to the node it is mapped to, and
        # within check_tol (relative) of it.
        m = samples.mean(0)
        D = np.abs(m[:, None] - fig_mean[None, :])
        own = np.diag(D)
        maxdiff = float(own.max())
        # nodes whose means are nearly tied (very short gene-tree branches) may be
        # closer to a neighbour by rounding; allow 0.5% of the root age of slack
        slack = 0.005 * fig_mean.max()
        if np.any(own > D.min(1) + slack) or np.any(own > check_tol * np.maximum(fig_mean, 0.1)):
            raise RunError(f"mean_check_failed(maxdiff={maxdiff:.4g})")
        if discard_frac > 0:
            samples = samples[int(len(samples) * discard_frac):]
    return tree, internals, fig_mean, hpd, samples, maxdiff


# --------------------------------------------------------------------------
# Per-gene worker
# --------------------------------------------------------------------------

G = {}  # globals for workers


def _init(state):
    G.update(state)


def lca_matrices(tree, internals, tax_index, ntax):
    """A[i,j] = internal node idx of LCA(i,j); D[i,j] = its depth (edges from root)."""
    A = np.full((ntax, ntax), -1, dtype=np.int32)
    D = np.full((ntax, ntax), -1, dtype=np.int32)
    depth = {tree: 0}
    for nd in tree.preorder():
        for c in nd.children:
            depth[c] = depth[nd] + 1
    leafsets = {}
    for nd in tree.postorder():
        if nd.is_leaf():
            leafsets[nd] = [tax_index[nd.name]]
        else:
            ls = [leafsets[c] for c in nd.children]
            for x, y in itertools.combinations(range(len(ls)), 2):
                ix = np.array(ls[x]); iy = np.array(ls[y])
                A[np.ix_(ix, iy)] = nd.idx; A[np.ix_(iy, ix)] = nd.idx
                D[np.ix_(ix, iy)] = depth[nd]; D[np.ix_(iy, ix)] = depth[nd]
            leafsets[nd] = sum(ls, [])
    return A, D


def node_stats(samples, internals, fig_mean, hpd, keep_q):
    k = len(internals)
    out = {"node": np.array([n.idx for n in internals]),
           "mean": fig_mean, "hpd_lo": hpd[:, 0], "hpd_hi": hpd[:, 1]}
    if samples is not None:
        out["median"] = np.median(samples, 0)
        out["var"] = samples.var(0)
        out["ess"] = ess(samples)
        if keep_q:
            q = np.quantile(samples, QLEVELS, axis=0).T
            for j, qc in enumerate(QCOLS):
                out[qc] = q[:, j]
    else:
        out["median"] = np.full(k, np.nan)
        out["var"] = np.full(k, np.nan)
        out["ess"] = np.full(k, np.nan)
    return out


def process_gene(item):
    gene, dfig, dmc, pfig, pmc = item
    st = G
    qc = {"gene": gene, "status": "ok", "reason": "", "n_taxa": np.nan,
          "min_ess_data_used": np.nan, "min_ess_prior_used": np.nan,
          "mean_check_data": np.nan, "mean_check_prior": np.nan, "has_prior": pfig is not None}
    try:
        tree, internals, fmean, hpd, smp, md = read_run(dfig, dmc, st["discard"], st["tol"])
    except RunError as e:
        qc["status"], qc["reason"] = "excluded", f"data:{e}"
        return qc, None, None, None
    except Exception as e:  # parse errors
        qc["status"], qc["reason"] = "excluded", f"data:parse_error:{e}"
        return qc, None, None, None
    qc["mean_check_data"] = md
    tips = [lf.name for lf in tree.leaves()]
    qc["n_taxa"] = len(tips)
    unknown = set(tips) - set(st["tax_index"])
    if unknown:
        qc["status"], qc["reason"] = "excluded", f"unknown_taxa:{sorted(unknown)[:3]}"
        return qc, None, None, None
    og = st["outgroup"]
    if og not in tips:
        qc["status"], qc["reason"] = "excluded", "no_outgroup"
        return qc, None, None, None
    if not any(c.is_leaf() and c.name == og for c in tree.children):
        qc["status"], qc["reason"] = "excluded", "outgroup_not_at_root"
        return qc, None, None, None
    if smp is None and st["require_samples"]:
        qc["status"], qc["reason"] = "excluded", "data:no_mcmc_samples"
        return qc, None, None, None

    # prior run
    psmp = None; pinternals = None
    if pfig is not None:
        try:
            ptree, pinternals, pfmean, phpd, psmp, pmd = read_run(pfig, pmc, st["discard"], st["tol"])
            qc["mean_check_prior"] = pmd
        except Exception as e:
            qc["status"], qc["reason"] = "excluded", f"prior:{e}"
            return qc, None, None, None
        dcl = [frozenset(n.leaf_names()) for n in internals]
        pcl = [frozenset(n.leaf_names()) for n in pinternals]
        if dcl != pcl:
            qc["status"], qc["reason"] = "excluded", "prior_topology_mismatch"
            return qc, None, None, None
    elif st["require_prior"]:
        qc["status"], qc["reason"] = "excluded", "no_prior_run"
        return qc, None, None, None

    tax_index = st["tax_index"]; ntax = len(tax_index)
    A, D = lca_matrices(tree, internals, tax_index, ntax)
    present = np.zeros(ntax, bool)
    present[[tax_index[t] for t in tips]] = True

    pair_rows = []
    trip_rows = []
    for br in st["branches"]:
        a, b, s = st["trip_idx"][br.bid]
        v = present[a] & present[b] & present[s]
        if not v.any():
            continue
        a_, b_, s_ = a[v], b[v], s[v]
        dd = np.stack([D[a_, b_], D[a_, s_], D[b_, s_]], 1)
        aa = np.stack([A[a_, b_], A[a_, s_], A[b_, s_]], 1)
        cls = dd.argmax(1)
        # the sister pair is strictly deeper than the other two pairs
        srt = np.sort(dd, 1)
        ok = srt[:, 2] > srt[:, 1]
        if not ok.all():
            qc["reason"] += f"{br.bid}:unresolved_triplets;"
        rows = np.arange(len(cls))
        nin = aa[rows, cls]
        nout = aa[rows, (cls + 1) % 3]
        assert np.all(nout == aa[rows, (cls + 2) % 3]), "outer LCA mismatch"
        key = np.stack([cls[ok], nin[ok], nout[ok]], 1)
        uk, cnt = np.unique(key, axis=0, return_counts=True)
        n_g = int(ok.sum())
        for (c, i_, o_), k in zip(uk, cnt):
            pair_rows.append((gene, br.bid, br.classes[c], int(i_), int(o_), int(k), n_g))
        if br.bid in st["trip_branches"]:
            names = st["taxa"]
            for t_a, t_b, t_s, c, i_, o_ in zip(a_[ok], b_[ok], s_[ok], cls[ok], nin[ok], nout[ok]):
                trip_rows.append((gene, br.bid, names[t_a], names[t_b], names[t_s],
                                  br.classes[c], int(i_), int(o_)))

    if not pair_rows:
        qc["status"], qc["reason"] = "excluded", "no_branch_specific_triplets"
        return qc, None, None, None

    pairs = pd.DataFrame(pair_rows, columns=["gene", "branch", "cls", "node_in", "node_out",
                                             "n_trip", "n_gene_branch"])
    # node tables
    nd = pd.DataFrame(node_stats(smp, internals, fmean, hpd, keep_q=True))
    nd.insert(0, "run", "data"); nd.insert(0, "gene", gene)
    nd["n_desc"] = [len(n.leaf_names()) for n in internals]
    nd["clade"] = [",".join(sorted(n.leaf_names())) if st["write_clades"] else "" for n in internals]
    frames = [nd]
    if pinternals is not None:
        pn = pd.DataFrame(node_stats(psmp, pinternals, pfmean, phpd, keep_q=False))
        pn.insert(0, "run", "prior"); pn.insert(0, "gene", gene)
        pn["n_desc"] = nd["n_desc"].values; pn["clade"] = ""
        frames.append(pn)
    nodes = pd.concat(frames, ignore_index=True).reindex(columns=NODE_COLS)

    # ESS of the nodes actually used, and posterior-median r per node pair
    used = np.unique(np.concatenate([pairs.node_in.values, pairs.node_out.values]))
    if smp is not None:
        e = nd.set_index("node")["ess"]
        qc["min_ess_data_used"] = float(e.loc[used].min())
    if psmp is not None:
        e = pn.set_index("node")["ess"]
        qc["min_ess_prior_used"] = float(e.loc[used].min())
    up = pairs[["node_in", "node_out"]].drop_duplicates()
    rd, rp = {}, {}
    for i_, o_ in up.itertuples(index=False):
        if smp is not None:
            rd[(i_, o_)] = float(np.median(smp[:, i_] / smp[:, o_]))
        else:
            rd[(i_, o_)] = float(fmean[i_] / fmean[o_])
        if psmp is not None:
            rp[(i_, o_)] = float(np.median(psmp[:, i_] / psmp[:, o_]))
        elif pinternals is not None:
            rp[(i_, o_)] = float(pfmean[i_] / pfmean[o_])
        else:
            rp[(i_, o_)] = np.nan
    keys = list(zip(pairs.node_in, pairs.node_out))
    pairs["r_data"] = [rd[k] for k in keys]
    pairs["r_prior"] = [rp[k] for k in keys]
    trips = pd.DataFrame(trip_rows, columns=["gene", "branch", "tax_K1", "tax_K2", "tax_S",
                                             "cls", "node_in", "node_out"]) if trip_rows else None
    return qc, nodes, pairs, trips


# --------------------------------------------------------------------------
# File discovery
# --------------------------------------------------------------------------


def discover(root: Path, fig_suffix: str, mcmc_suffix: str, id_regex: str | None,
             skip_prior_tag: bool = True):
    """Map gene ID -> (FigTree.tre, mcmc.txt or None). Files named
    <stem>.prior.FigTree.tre (aeon --prior-only without --out-dir) are skipped
    unless fig_suffix itself is the prior suffix."""
    out = {}
    if root is None:
        return out
    rx = re.compile(id_regex) if id_regex else None
    for fig in sorted(root.rglob(f"*{fig_suffix}")):
        stem = fig.name[: -len(fig_suffix)]
        if skip_prior_tag and stem.endswith(".prior") and not fig_suffix.startswith(".prior"):
            continue
        if any(part.startswith("_aeon_work_") for part in fig.parts):
            continue
        gid = stem
        if rx:
            m = rx.search(stem)
            if not m:
                continue
            gid = m.group(1) if m.groups() else m.group(0)
        mc = fig.with_name(stem + mcmc_suffix)
        if gid in out:
            raise SystemExit(f"duplicate gene id {gid}: {out[gid][0]} and {fig}")
        out[gid] = (fig, mc if mc.exists() else None)
    return out


def add_arguments(ap):
    ap.add_argument("--species-tree", required=True,
                    help="ASTRAL species tree (CU internal branch lengths; q1/q2/q3 annotations optional)")
    ap.add_argument("--data-dir", required=True, help="root of the data-run outputs (searched recursively)")
    ap.add_argument("--prior-dir", help="root of the prior-only outputs (searched recursively)")
    ap.add_argument("--outgroup", default="Gallus_gallus")
    ap.add_argument("--gene-list", help="file with the final locus IDs (Section 4.1); others ignored")
    ap.add_argument("--gene-id-regex", help="regex applied to '<stem>' of <stem>.FigTree.tre; "
                                            "group 1 (or whole match) is the gene ID")
    ap.add_argument("--fig-suffix", default=".FigTree.tre")
    ap.add_argument("--mcmc-suffix", default=".mcmc.txt")
    ap.add_argument("--polytomy", help="TSV: group<TAB>comma-separated taxa, three rows (Section 13)")
    ap.add_argument("--polytomy-T", type=float, default=float("nan"),
                    help="CU length to report for the polytomy node (e.g. 0.016)")
    ap.add_argument("--clade-names", help="TSV: name<TAB>comma-separated taxa, for readable branch labels")
    ap.add_argument("--triplet-branches", default="",
                    help="comma-separated branch ids (e.g. poly,b05) for which to write gene x triplet rows; "
                         "run once without it to see ids in branches.tsv")
    ap.add_argument("--discard-frac", type=float, default=0.0,
                    help="extra burn-in fraction to drop from mcmc.txt (MCMCTree already drops burnin)")
    ap.add_argument("--mean-check-tol", type=float, default=0.05,
                    help="relative tolerance for the mcmc.txt vs FigTree.tre node-mapping check")
    ap.add_argument("--allow-no-prior", action="store_true",
                    help="keep genes without a prior-only run (DiD will be NaN for them)")
    ap.add_argument("--allow-no-samples", action="store_true",
                    help="keep genes without mcmc.txt (uses FigTree means; no ESS, no medians)")
    ap.add_argument("--no-clades", action="store_true", help="do not store node clades in gene_nodes")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out-dir", default="ilsdate_out")


def run(a):

    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    sp = read_tree_file(a.species_tree)
    taxa = sorted(sp.leaf_names())
    if a.outgroup not in taxa:
        sys.exit(f"outgroup {a.outgroup} not in species tree")
    ingroup = set(taxa) - {a.outgroup}
    branches, warns = define_branches(sp, a.outgroup, ingroup, read_polytomy(a.polytomy),
                                      read_clade_names(a.clade_names), a.polytomy_T)
    for w in warns:
        print("[WARN]", w)
    tax_index = {t: i for i, t in enumerate(taxa)}
    trip_idx = {}
    for br in branches:
        g = [np.array([tax_index[t] for t in grp]) for grp in br.groups]
        mesh = np.array(np.meshgrid(*g, indexing="ij")).reshape(3, -1)
        trip_idx[br.bid] = (mesh[0], mesh[1], mesh[2])

    # branches.tsv
    rows = []
    for br in branches:
        pred = msc_expectations(br.T) if br.kind == "binary" else msc_expectations(np.nan)
        tin_taxa, tout_taxa = species_node_for_branch(br)
        rows.append(dict(branch=br.bid, kind=br.kind, label=br.label(), T_CU=br.T,
                         n_triplets=br.n_triplets, classes=",".join(br.classes),
                         group1=",".join(br.groups[0]), group2=",".join(br.groups[1]),
                         group3=",".join(br.groups[2]),
                         tau_in_taxa=",".join(tin_taxa), tau_out_taxa=",".join(tout_taxa),
                         **{f"pred_{k}": v for k, v in pred.items()},
                         **{f"q_astral_{c}": br.q_astral.get(c, np.nan) for c in br.classes[:3]},
                         q_note=br.q_note))
    bdf = pd.DataFrame(rows)
    bdf.to_csv(out / "branches.tsv", sep="\t", index=False)
    print(f"[INFO] {len(branches)} testable branches; {bdf.n_triplets.sum()} branch-specific triplets "
          f"(of {len(ingroup) * (len(ingroup) - 1) * (len(ingroup) - 2) // 6} lepidosaur triplets)")

    dmap = discover(Path(a.data_dir), a.fig_suffix, a.mcmc_suffix, a.gene_id_regex)
    pmap = {}
    if a.prior_dir:
        # aeon --prior-only writes <stem>.prior.* unless given --out-dir
        pmap = discover(Path(a.prior_dir), ".prior" + a.fig_suffix, ".prior" + a.mcmc_suffix, a.gene_id_regex)
        if not pmap:
            pmap = discover(Path(a.prior_dir), a.fig_suffix, a.mcmc_suffix, a.gene_id_regex)
    genes = sorted(dmap)
    if a.gene_list:
        keep = {ln.strip() for ln in Path(a.gene_list).read_text().splitlines() if ln.strip()}
        missing = sorted(keep - set(dmap))
        if missing:
            print(f"[WARN] {len(missing)} listed genes have no data-run output (e.g. {missing[:3]})")
            (out / "genes_missing_output.txt").write_text("\n".join(missing) + "\n")
        genes = [g for g in genes if g in keep]
    print(f"[INFO] {len(genes)} genes with data-run output; {sum(g in pmap for g in genes)} with prior-only output")

    items = [(g, dmap[g][0], dmap[g][1], pmap.get(g, (None, None))[0], pmap.get(g, (None, None))[1])
             for g in genes]
    state = dict(branches=branches, trip_idx=trip_idx, tax_index=tax_index, taxa=taxa,
                 outgroup=a.outgroup, discard=a.discard_frac, tol=a.mean_check_tol,
                 require_prior=not a.allow_no_prior, require_samples=not a.allow_no_samples,
                 trip_branches=set(x for x in a.triplet_branches.split(",") if x),
                 write_clades=not a.no_clades)

    qcs = []
    node_path = out / "gene_nodes.tsv.gz"
    pair_path = out / "gene_branch_pairs.tsv.gz"
    trip_files = {}
    first_n = first_p = True
    with gzip.open(node_path, "wt") as fn, gzip.open(pair_path, "wt") as fp:
        if a.threads > 1:
            pool = mp.Pool(a.threads, initializer=_init, initargs=(state,))
            it = pool.imap_unordered(process_gene, items, chunksize=8)
        else:
            _init(state); it = map(process_gene, items)
        for i, (qc, nodes, pairs, trips) in enumerate(it, 1):
            qcs.append(qc)
            if nodes is not None:
                nodes.to_csv(fn, sep="\t", index=False, header=first_n, float_format="%.6g"); first_n = False
                pairs.to_csv(fp, sep="\t", index=False, header=first_p, float_format="%.6g"); first_p = False
            if trips is not None:
                for bid, sub in trips.groupby("branch"):
                    if bid not in trip_files:
                        trip_files[bid] = gzip.open(out / f"triplets_{bid}.tsv.gz", "wt")
                        sub.to_csv(trip_files[bid], sep="\t", index=False, header=True)
                    else:
                        sub.to_csv(trip_files[bid], sep="\t", index=False, header=False)
            if i % 500 == 0:
                print(f"[INFO] {i}/{len(items)} genes")
        if a.threads > 1:
            pool.close(); pool.join()
    for f in trip_files.values():
        f.close()
    qdf = pd.DataFrame(qcs).sort_values("gene")
    qdf.to_csv(out / "gene_qc.tsv", sep="\t", index=False)
    print(qdf.groupby(["status"]).size().to_string())
    bad = qdf[qdf.status != "ok"]
    if len(bad):
        print("[INFO] exclusion reasons:")
        print(bad.reason.str.replace(r"[(/].*", "", regex=True).value_counts().to_string())
    print(f"[DONE] wrote {out}/")



def register(sub) -> None:
    p = sub.add_parser('ils-extract', help='Extract triplet node ages from dated gene trees (ILS dating test)', description=__doc__,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(p)
    p.set_defaults(func=run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return run(ap.parse_args(argv))
