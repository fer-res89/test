"""Sections 8-13 and 16 of gene_tree_dating_ILS_test.md, from the tables written
by `aeon ils-extract`.

Outputs (in --out-dir):
  gene_accounting.tsv     genes by QC status and ESS exclusion (Table 5)
  branch_overview.tsv     per branch: T, triplets, genes, F(c), predicted Delta, information ratios
  astral_check.tsv        F(c)/sum F against ASTRAL q (Table 1b)
  branch_contrasts.tsv    per branch x contrast x quantity: data, prior, DiD with 95% CI,
                          bootstrap p, BH q (DiD and data families), weighted medians (Tables 1, 4)
  absolute_test.tsv       fraction of genes with t_in > tau_out, per tau source (Table 2)
  implied_2Neg.tsv        2Ne*g from Delta t_in, Delta t_out (data and DiD) and from BPP (Table 3)
  regression.tsv          observed Delta (Ma) on predicted Delta (CU) across branches, bootstrapped
  per_triplet_<b>.tsv     descriptive per-triplet comparisons (Section 8), if triplet tables exist
  gene_branch_class.tsv.gz  gene x branch x class weights and class-mean ages (all quantities)
  dates_long.tsv.gz       tidy gene-tree dates for plotting: gene x branch x class x run,
                          t_in_Ma, t_out_Ma, r, weight w (R-ready)
  branch_tau.tsv          species-tree tau_in / tau_out per branch and tau source
  plot_branch_distributions.R   ggplot2 script: faceted histograms from the two tables above
  branch_distributions.pdf, branch_distributions/<branch>.png
                          per-branch histograms of t_in, t_out, r by class (data + prior)
  fig1_forest.png, fig2_obs_vs_pred.png, fig4_absolute.png

Weighting (Section 9.3): w(g,c) = n(g,c)/n(g); class age = mean over the gene's
class-c triplets; mu(c) = sum_g w a / sum_g w. Pooled D sums the D1 and D2
contributions. All intervals: gene-level bootstrap (genes resampled with
replacement; data and prior statistics recomputed on the same replicate, so the
DiD is paired). Two-sided bootstrap p = 2 min(P*(<=0), P*(>=0)) with the
(k+1)/(B+1) correction, so the smallest attainable p is about 2/(B+1).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .ils_common import mrca_age, msc_expectations, node_ages, read_tree_file

QTY = ["tin", "tout", "r"]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def wquantile(x, w, q):
    x = np.asarray(x, float); w = np.asarray(w, float)
    m = np.isfinite(x) & (w > 0)
    if not m.any():
        return np.nan
    o = np.argsort(x[m]); xs = x[m][o]; cw = np.cumsum(w[m][o])
    return float(xs[np.searchsorted(cw, q * cw[-1])])


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    m = np.isfinite(p)
    if not m.any():
        return out
    pv = p[m]; n = len(pv); o = np.argsort(pv)
    adj = pv[o] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(n); res[o] = np.minimum(adj, 1)
    out[m] = res
    return out


def boot_p(rep):
    rep = rep[np.isfinite(rep)]
    B = len(rep)
    if B == 0:
        return np.nan
    lo = ((rep <= 0).sum() + 1) / (B + 1)
    hi = ((rep >= 0).sum() + 1) / (B + 1)
    return float(min(1.0, 2 * min(lo, hi)))


def parse_tau(spec: str):
    """NAME=PATH, NAME=PATH@root=280 (rescale so the tree root is 280 Ma),
    NAME=PATH@x100 (multiply ages by 100)."""
    name, rest = spec.split("=", 1)
    path, mod = (rest.split("@", 1) + [None])[:2]
    tr = read_tree_file(path)
    ages = node_ages(tr)
    if mod:
        if mod.startswith("root="):
            f = float(mod[5:]) / ages[tr]
        elif mod.startswith("x"):
            f = float(mod[1:])
        else:
            raise SystemExit(f"bad tau modifier {mod}")
        ages = {k: v * f for k, v in ages.items()}
    return name, tr, ages


# --------------------------------------------------------------------------
# Building the gene x branch x class table
# --------------------------------------------------------------------------

def load(a):
    ext = Path(a.ext_dir)
    br = pd.read_csv(ext / "branches.tsv", sep="\t")
    qc = pd.read_csv(ext / "gene_qc.tsv", sep="\t", dtype={"reason": str})
    qc["reason"] = qc.reason.fillna("")
    pairs = pd.read_csv(ext / "gene_branch_pairs.tsv.gz", sep="\t")
    use = ["gene", "run", "node", "median", "mean", "var"]
    qcols = None
    hdr = pd.read_csv(ext / "gene_nodes.tsv.gz", sep="\t", nrows=0).columns
    qcols = [c for c in hdr if c.startswith("q") and c[1:].isdigit()]
    nodes = pd.read_csv(ext / "gene_nodes.tsv.gz", sep="\t", usecols=use + qcols)
    return br, qc, pairs, nodes, qcols


def gene_filter(qc, a):
    qc = qc.copy()
    qc["use"] = qc.status.eq("ok")
    qc["ess_fail_data"] = qc.use & (qc.min_ess_data_used < a.min_ess)
    qc["ess_fail_prior"] = qc.use & (qc.min_ess_prior_used < a.min_ess_prior)
    qc.loc[qc.ess_fail_data | qc.ess_fail_prior, "use"] = False
    if a.exclude_genes:
        ex = {ln.strip() for ln in Path(a.exclude_genes).read_text().splitlines() if ln.strip()}
        qc["excluded_list"] = qc.gene.isin(ex)
        qc.loc[qc.excluded_list, "use"] = False
    return qc


def attach_ages(pairs, nodes, a, qcols):
    stat = a.age_stat
    d = nodes[nodes.run == "data"].set_index(["gene", "node"])
    p = nodes[nodes.run == "prior"].set_index(["gene", "node"])
    out = pairs.copy()
    ki = pd.MultiIndex.from_arrays([out.gene, out.node_in])
    ko = pd.MultiIndex.from_arrays([out.gene, out.node_out])
    U = a.time_unit_ma
    out["tin_d"] = d[stat].reindex(ki).to_numpy() * U
    out["tout_d"] = d[stat].reindex(ko).to_numpy() * U
    out["tin_p"] = p[stat].reindex(ki).to_numpy() * U
    out["tout_p"] = p[stat].reindex(ko).to_numpy() * U
    out["r_d"] = out.r_data
    out["r_p"] = out.r_prior
    vd_i = d["var"].reindex(ki).to_numpy(); vp_i = p["var"].reindex(ki).to_numpy()
    vd_o = d["var"].reindex(ko).to_numpy(); vp_o = p["var"].reindex(ko).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        out["info_in"] = 1 - vd_i / vp_i
        out["info_out"] = 1 - vd_o / vp_o
    qin = (d[qcols].reindex(ki).to_numpy(np.float32) * U) if qcols else None
    return out, qin


def post_prob_greater(qmat, tau):
    """P(t > tau) from a row-wise table of posterior quantiles at levels 0..1
    (linear interpolation of the empirical CDF; 0/1 outside the sampled range)."""
    lv = np.linspace(0, 1, qmat.shape[1])
    tau = np.asarray(tau, float)
    out = np.full(len(qmat), np.nan)
    ok = np.isfinite(tau) & np.isfinite(qmat).all(1)
    q = qmat[ok]; t = tau[ok]
    j = (q < t[:, None]).sum(1)                       # number of quantiles below tau
    res = np.empty(len(t))
    res[j == 0] = 1.0
    top = j == q.shape[1]
    res[top] = 0.0
    mid = ~(top | (j == 0))
    jm = j[mid]; r = np.nonzero(mid)[0]
    q0, q1 = q[r, jm - 1], q[r, jm]
    frac = np.where(q1 > q0, (t[mid] - q0) / np.where(q1 > q0, q1 - q0, 1), 0.0)
    res[mid] = 1 - (lv[jm - 1] + frac * (lv[jm] - lv[jm - 1]))
    out[ok] = res
    return out


def aggregate(pairs, quantities, br):
    """gene x branch x class: w and w*abar for each quantity; pooled D added."""
    df = pairs.copy()
    for q in quantities:
        df[f"_s_{q}"] = df[q] * df.n_trip
        df[f"_k_{q}"] = df[q].notna().astype(int)
    agg = df.groupby(["gene", "branch", "cls"], sort=False).agg(
        n_gc=("n_trip", "sum"), n_g=("n_gene_branch", "first"),
        **{f"_s_{q}": (f"_s_{q}", "sum") for q in quantities},
        **{f"_k_{q}": (f"_k_{q}", "sum") for q in quantities}).reset_index()
    for q in quantities:  # pandas sums all-NaN groups to 0; keep them missing
        agg.loc[agg[f"_k_{q}"] == 0, f"_s_{q}"] = np.nan
    agg["w"] = agg.n_gc / agg.n_g
    for q in quantities:
        agg[f"a_{q}"] = agg[f"_s_{q}"] / agg.n_gc          # class mean age of the gene
        agg[f"wa_{q}"] = agg.w * agg[f"a_{q}"]
    agg = agg.drop(columns=[f"_s_{q}" for q in quantities] + [f"_k_{q}" for q in quantities])
    binary = set(br.loc[br.kind == "binary", "branch"])
    dd = agg[agg.branch.isin(binary) & agg.cls.isin(["D1", "D2"])]
    if len(dd):
        pooled = dd.groupby(["gene", "branch"], sort=False).agg(
            n_gc=("n_gc", "sum"), n_g=("n_g", "first"), w=("w", "sum"),
            **{f"wa_{q}": (f"wa_{q}", lambda x: x.sum(min_count=1)) for q in quantities}).reset_index()
        pooled["cls"] = "D"
        for q in quantities:
            pooled[f"a_{q}"] = pooled[f"wa_{q}"] / pooled.w
        agg = pd.concat([agg, pooled], ignore_index=True)
    return agg


# --------------------------------------------------------------------------
# Statistics computed identically on the point estimate and on each replicate
# --------------------------------------------------------------------------

class StatEngine:
    def __init__(self, agg, br, quantities, min_F):
        self.genes = np.array(sorted(agg.gene.unique()))
        gi = {g: i for i, g in enumerate(self.genes)}
        self.keys = sorted(set(zip(agg.branch, agg.cls)))
        ki = {k: i for i, k in enumerate(self.keys)}
        G, K = len(self.genes), len(self.keys)
        rows = agg.gene.map(gi).to_numpy(); cols = np.array([ki[k] for k in zip(agg.branch, agg.cls)])
        self.W = np.zeros((G, K))
        self.W[rows, cols] = agg.w.to_numpy()
        self.WA = {}
        for q in quantities:
            m = np.zeros((G, K))
            m[rows, cols] = np.nan_to_num(agg[f"wa_{q}"].to_numpy())
            # weight that has a finite value (e.g. prior missing): track separately
            wv = np.zeros((G, K))
            wv[rows, cols] = np.where(np.isfinite(agg[f"wa_{q}"].to_numpy()), agg.w.to_numpy(), 0)
            self.WA[q] = (m, wv)
        self.quantities = quantities
        self.br = br.set_index("branch")
        self.min_F = min_F
        self.k = ki

    def mu(self, Mw):
        """Mw: (b, G) gene multiplicities. Returns dict q -> (b, K) means and F (b, K)."""
        F = Mw @ self.W
        out = {}
        for q, (m, wv) in self.WA.items():
            with np.errstate(invalid="ignore", divide="ignore"):
                out[q] = (Mw @ m) / (Mw @ wv)
        return out, F

    def col(self, b, c):
        return self.k.get((b, c))


def contrasts_spec(br):
    """List of (branch, contrast name, class A, class B): statistic = mu(A) - mu(B)."""
    spec = []
    for _, r in br.iterrows():
        cl = r.classes.split(",")
        if r.kind == "binary":
            spec += [(r.branch, "D-C", "D", "C"), (r.branch, "D1-D2", "D1", "D2"),
                     (r.branch, "D1-C", "D1", "C"), (r.branch, "D2-C", "D2", "C")]
        else:
            spec += [(r.branch, f"{cl[0]} - {cl[1]}", cl[0], cl[1]),
                     (r.branch, f"{cl[0]} - {cl[2]}", cl[0], cl[2]),
                     (r.branch, f"{cl[1]} - {cl[2]}", cl[1], cl[2])]
    return spec


def compute_all(eng: StatEngine, Mw, spec, abs_specs, pred, eligible):
    """Return dict name -> (b,) arrays."""
    mus, F = eng.mu(Mw)
    res = {}
    for (b, name, A, B) in spec:
        ia, ib = eng.col(b, A), eng.col(b, B)
        for q in QTY:
            for ver, suf in (("data", "d"), ("prior", "p")):
                key = f"{q}_{suf}"
                if ia is None or ib is None or key not in mus:
                    v = np.full(Mw.shape[0], np.nan)
                else:
                    v = mus[key][:, ia] - mus[key][:, ib]
                res[(b, name, q, ver)] = v
            res[(b, name, q, "did")] = res[(b, name, q, "data")] - res[(b, name, q, "prior")]
    # absolute test fractions
    for (b, c, src, method, qname) in abs_specs:
        i = eng.col(b, c)
        res[("abs", b, c, src, method)] = (mus[qname][:, i] if i is not None
                                           else np.full(Mw.shape[0], np.nan))
    # cross-branch regressions of observed Delta (Ma) on predicted Delta (CU)
    for q, pkey in (("tin", "pred_dTin_CU"), ("tout", "pred_dTout_CU")):
        bs = [b for b in eligible if np.isfinite(pred.loc[b, pkey])]
        if len(bs) < 3:
            continue
        x = pred.loc[bs, pkey].to_numpy()
        for ver in ("data", "did"):
            Y = np.stack([res[(b, "D-C", q, ver)] for b in bs], 1)
            ok = np.isfinite(Y).all(1)
            xm = x - x.mean()
            slope = np.where(ok, (Y - Y.mean(1, keepdims=True)) @ xm / (xm @ xm), np.nan)
            inter = Y.mean(1) - slope * x.mean()
            origin = Y @ x / (x @ x)
            res[("reg", q, ver, "slope")] = slope
            res[("reg", q, ver, "intercept")] = inter
            res[("reg", q, ver, "slope_origin")] = origin
    return res, mus, F


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def add_arguments(ap):
    ap.add_argument("--ext-dir", required=True, help="output directory of `aeon ils-extract`")
    ap.add_argument("--out-dir", default="ilsdate_results")
    ap.add_argument("--time-unit-ma", type=float, default=100.0, help="Myr per MCMCTree time unit")
    ap.add_argument("--age-stat", choices=["median", "mean"], default="median")
    ap.add_argument("--min-ess", type=float, default=200, help="data-run ESS threshold (Section 5.6)")
    ap.add_argument("--min-ess-prior", type=float, default=200, help="prior-run ESS threshold")
    ap.add_argument("--exclude-genes", help="file of gene IDs to exclude (e.g. failed second-chain check)")
    ap.add_argument("--min-F", type=float, default=10, help="minimum effective genes per class")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--chunk", type=int, default=100, help="replicates per matrix block")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--astral-tol", type=float, default=0.03, help="flag |F share - q| above this")
    ap.add_argument("--tau-tree", action="append", default=[],
                    help="species divergence times for the absolute test: NAME=PATH[@root=MA|@xFACTOR]; "
                         "repeatable (e.g. bpp_T1, bpp_T2, concat). Ages in Ma after the modifier.")
    ap.add_argument("--no-lower5", action="store_true", help="skip the gene-tree 5th-percentile tau source")
    ap.add_argument("--bpp-params", help="TSV with columns branch (or taxa), theta, tau_subs, tau_Ma "
                                         "for the triplet root node of each branch (Section 12)")
    ap.add_argument("--two-ne-g", type=float, help="Myr per CU, to overlay predicted Delta in Fig. 1")
    ap.add_argument("--plot-branches", default="",
                    help="comma-separated branch ids for the per-branch distribution plots (default: all)")


def run(a):
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)

    br, qc, pairs, nodes, qcols = load(a)
    qc = gene_filter(qc, a)
    acct = pd.DataFrame({
        "category": ["genes with data output", "excluded at extraction", "failed data ESS",
                     "failed prior ESS", "excluded by list", "used"],
        "n": [len(qc), int((qc.status != "ok").sum()), int(qc.ess_fail_data.sum()),
              int(qc.ess_fail_prior.sum()), int(qc.get("excluded_list", pd.Series(False)).sum()),
              int(qc.use.sum())]})
    reasons = qc.loc[qc.status != "ok", "reason"].str.replace(r"[(/].*", "", regex=True).value_counts()
    acct = pd.concat([acct, pd.DataFrame({"category": [f"  extraction: {k}" for k in reasons.index],
                                          "n": reasons.values})])
    acct.to_csv(out / "gene_accounting.tsv", sep="\t", index=False)
    print(acct.to_string(index=False))

    # ASTRAL check uses every gene with branch triplets (topology does not depend on dating)
    pairs_all = pairs.copy()
    pairs = pairs[pairs.gene.isin(set(qc.gene[qc.use]))]
    pairs, qin = attach_ages(pairs, nodes, a, qcols)

    br = br.copy()
    br_idx = br.set_index("branch")

    # ---- tau sources for the absolute test
    tau_sources = {}
    tau_in_sources = {}
    for spec in a.tau_tree:
        name, tr, ages = parse_tau(spec)
        tau_sources[name] = {r.branch: mrca_age(tr, ages, r.tau_out_taxa.split(",")) for r in br.itertuples()}
        tau_in_sources[name] = {r.branch: mrca_age(tr, ages, r.tau_in_taxa.split(",")) for r in br.itertuples()}
    quantities = ["tin_d", "tout_d", "r_d", "tin_p", "tout_p", "r_p", "info_in", "info_out"]

    # preliminary aggregate to get gene_lower5 tau (weighted 5th percentile of class t_out)
    agg0 = aggregate(pairs, ["tout_d"], br)
    if not a.no_lower5:
        sub = agg0[agg0.cls != "D"]
        tau_sources["gene_lower5"] = {b: wquantile(g.a_tout_d, g.w, 0.05) for b, g in sub.groupby("branch")}
    abs_q = []
    abs_specs = []
    for src, tmap in tau_sources.items():
        tau = pairs.branch.map(tmap).to_numpy(float)
        pairs[f"abs_pt_{src}"] = np.where(np.isfinite(tau), (pairs.tin_d > tau).astype(float), np.nan)
        if qin is not None:
            pairs[f"abs_post_{src}"] = post_prob_greater(qin, tau)
            abs_q += [f"abs_pt_{src}", f"abs_post_{src}"]
        else:
            abs_q += [f"abs_pt_{src}"]
    quantities += abs_q

    agg = aggregate(pairs, quantities, br)
    agg.to_csv(out / "gene_branch_class.tsv.gz", sep="\t", index=False, float_format="%.6g")
    write_dates_tables(agg, br, tau_sources, tau_in_sources, out)

    for r in br.itertuples():
        for c in r.classes.split(",") + (["D"] if r.kind == "binary" else []):
            for src in tau_sources:
                for meth in ("pt", "post"):
                    qn = f"abs_{meth}_{src}"
                    if qn in quantities:
                        abs_specs.append((r.branch, c, src, meth, qn))

    eng = StatEngine(agg, br, quantities, a.min_F)
    G = len(eng.genes)
    spec = contrasts_spec(br)
    pred = br_idx[["pred_dTin_CU", "pred_dTout_CU", "T_CU"]]

    # point estimates
    ones = np.ones((1, G))
    mus0, F0 = eng.mu(ones)
    Fpt = {k: F0[0, i] for k, i in eng.k.items()}
    eligible = [b for b in br.branch[br.kind == "binary"]
                if min(Fpt.get((b, "C"), 0), Fpt.get((b, "D"), 0)) >= a.min_F]
    point, _, _ = compute_all(eng, ones, spec, abs_specs, pred, eligible)

    # bootstrap
    rng = np.random.default_rng(a.seed)
    reps = {k: [] for k in point}
    done = 0
    while done < a.bootstrap:
        b = min(a.chunk, a.bootstrap - done)
        Mw = rng.multinomial(G, np.full(G, 1.0 / G), size=b).astype(float)
        r, _, _ = compute_all(eng, Mw, spec, abs_specs, pred, eligible)
        for k, v in r.items():
            reps[k].append(v)
        done += b
    reps = {k: np.concatenate(v) for k, v in reps.items()}

    def ci(k):
        v = reps[k]; v = v[np.isfinite(v)]
        if len(v) < 0.9 * a.bootstrap:
            return np.nan, np.nan
        return tuple(np.percentile(v, [2.5, 97.5]))

    # ---- branch overview
    ov = []
    for r in br.itertuples():
        row = dict(branch=r.branch, kind=r.kind, label=r.label, T_CU=r.T_CU, n_triplets=r.n_triplets,
                   n_genes=int(agg[(agg.branch == r.branch) & (agg.cls != "D")].gene.nunique()))
        cls = r.classes.split(",")
        for c in cls + (["D"] if r.kind == "binary" else []):
            row[f"F_{c}"] = Fpt.get((r.branch, c), 0.0)
        for c in cls:
            i = eng.col(r.branch, c)
            row[f"info_in_{c}"] = mus0["info_in"][0, i] if i is not None else np.nan
            row[f"info_out_{c}"] = mus0["info_out"][0, i] if i is not None else np.nan
        row.update(pred_P_C=r.pred_P_C, pred_dTin_CU=r.pred_dTin_CU, pred_dTout_CU=r.pred_dTout_CU,
                   pred_frac_C_older=r.pred_frac_C_older)
        ov.append(row)
    ov = pd.DataFrame(ov)
    ov.to_csv(out / "branch_overview.tsv", sep="\t", index=False, float_format="%.5g")

    # ---- ASTRAL check (all genes with triplets)
    ac = []
    pa = pairs_all.groupby(["gene", "branch", "cls"]).agg(n=("n_trip", "sum"), ng=("n_gene_branch", "first"))
    pa["w"] = pa.n / pa.ng
    Fall = pa.groupby(["branch", "cls"]).w.sum()
    for r in br.itertuples():
        cls = r.classes.split(",")
        tot = sum(Fall.get((r.branch, c), 0) for c in cls)
        row = dict(branch=r.branch, label=r.label, n_genes_all=tot, q_note=r.q_note)
        dmax = np.nan
        for c in cls:
            share = Fall.get((r.branch, c), 0) / tot if tot else np.nan
            row[f"share_{c}"] = share
            qcol = f"q_astral_{c}"
            q = getattr(r, qcol) if hasattr(r, qcol) else np.nan
            row[f"q_{c}"] = q
            if np.isfinite(q):
                d = 100 * (share - q); row[f"diff_pp_{c}"] = d
                dmax = np.nanmax([dmax, abs(d)]) if np.isfinite(dmax) else abs(d)
        if r.kind == "binary":
            ex = msc_expectations(r.T_CU)
            row["msc_P_C"], row["msc_P_D_each"] = ex["P_C"], ex["P_D_each"]
        row["max_abs_diff_pp"] = dmax
        row["flag"] = bool(np.isfinite(dmax) and dmax > 100 * a.astral_tol)
        ac.append(row)
    pd.DataFrame(ac).to_csv(out / "astral_check.tsv", sep="\t", index=False, float_format="%.5g")

    # ---- branch contrasts
    rows = []
    wm = {}
    for r in br.itertuples():
        sub = agg[agg.branch == r.branch]
        for c in sub.cls.unique():
            s = sub[sub.cls == c]
            for q in QTY:
                wm[(r.branch, c, q, "data")] = wquantile(s[f"a_{q}_d"], s.w, 0.5)
                wm[(r.branch, c, q, "prior")] = wquantile(s[f"a_{q}_p"], s.w, 0.5)
    for (b, name, A, B) in spec:
        rb = br_idx.loc[b]
        FA, FB = Fpt.get((b, A), 0.0), Fpt.get((b, B), 0.0)
        for q in QTY:
            row = dict(branch=b, label=rb.label, T_CU=rb.T_CU, contrast=name, quantity=q,
                       F_A=FA, F_B=FB, reportable=bool(min(FA, FB) >= a.min_F))
            if name == "D-C" and q in ("tin", "tout"):
                row["pred_CU"] = rb.pred_dTin_CU if q == "tin" else rb.pred_dTout_CU
            for ver in ("data", "prior", "did"):
                k = (b, name, q, ver)
                lo, hi = ci(k)
                row[f"{ver}"] = point[k][0]
                row[f"{ver}_lo"], row[f"{ver}_hi"] = lo, hi
                row[f"{ver}_p"] = boot_p(reps[k])
            for ver in ("data", "prior"):
                row[f"{ver}_wmedian"] = wm.get((b, A, q, ver), np.nan) - wm.get((b, B, q, ver), np.nan)
            row["did_wmedian"] = row["data_wmedian"] - row["prior_wmedian"]
            rows.append(row)
    bc = pd.DataFrame(rows)
    for ver in ("did", "data"):
        bc[f"{ver}_q_bh"] = np.nan
        for (con, q), g in bc[bc.reportable].groupby(["contrast", "quantity"]):
            if con in ("D1-C", "D2-C"):
                continue  # descriptive only; not a planned family
            bc.loc[g.index, f"{ver}_q_bh"] = bh(g[f"{ver}_p"].to_numpy())
    bc.to_csv(out / "branch_contrasts.tsv", sep="\t", index=False, float_format="%.5g")

    # ---- absolute test
    ab = []
    for (b, c, src, meth, qn) in abs_specs:
        k = ("abs", b, c, src, meth)
        lo, hi = ci(k)
        T = br_idx.loc[b, "T_CU"]
        predv = (1.0 if c in ("D", "D1", "D2") else msc_expectations(T)["frac_C_older"]) \
            if br_idx.loc[b, "kind"] == "binary" else np.nan
        ab.append(dict(branch=b, label=br_idx.loc[b, "label"], T_CU=T, cls=c, tau_source=src,
                       tau_out_Ma=tau_sources[src].get(b, np.nan), method=meth,
                       F=Fpt.get((b, c), 0.0), frac=point[k][0], lo=lo, hi=hi, predicted=predv))
    ab = pd.DataFrame(ab)
    ab.to_csv(out / "absolute_test.tsv", sep="\t", index=False, float_format="%.5g")

    # ---- implied 2Ne*g
    bpp = None
    if a.bpp_params:
        bpp = pd.read_csv(a.bpp_params, sep="\t")
        if "branch" not in bpp.columns:
            m = {frozenset(r.tau_out_taxa.split(",")): r.branch for r in br.itertuples()}
            bpp["branch"] = [m.get(frozenset(t.split(","))) for t in bpp.taxa]
        bpp["two_ne_g_bpp_Myr"] = bpp.theta * bpp.tau_Ma / (2 * bpp.tau_subs)
        bpp = bpp.set_index("branch")
    im = []
    for b in br.branch[br.kind == "binary"]:
        row = dict(branch=b, label=br_idx.loc[b, "label"], T_CU=br_idx.loc[b, "T_CU"],
                   eligible=b in eligible)
        for q, pk in (("tin", "pred_dTin_CU"), ("tout", "pred_dTout_CU")):
            pc = br_idx.loc[b, pk]
            for ver in ("data", "did"):
                k = (b, "D-C", q, ver)
                v = reps[k] / pc
                v = v[np.isfinite(v)]
                row[f"from_{q}_{ver}"] = point[k][0] / pc
                row[f"from_{q}_{ver}_lo"], row[f"from_{q}_{ver}_hi"] = (
                    np.percentile(v, [2.5, 97.5]) if len(v) else (np.nan, np.nan))
        if bpp is not None and b in bpp.index:
            row["from_bpp"] = bpp.loc[b, "two_ne_g_bpp_Myr"]
        im.append(row)
    pd.DataFrame(im).to_csv(out / "implied_2Neg.tsv", sep="\t", index=False, float_format="%.5g")
    rg = []
    for k in point:
        if k[0] == "reg":
            lo, hi = ci(k)
            rg.append(dict(quantity=k[1], version=k[2], parameter=k[3], estimate=point[k][0], lo=lo, hi=hi,
                           n_branches=len(eligible)))
    pd.DataFrame(rg).to_csv(out / "regression.tsv", sep="\t", index=False, float_format="%.5g")

    # ---- per-triplet descriptive comparisons (Section 8)
    per_triplet(Path(a.ext_dir), nodes, br_idx, a, out, set(qc.gene[qc.use]))

    # ---- figures
    try:
        figures(bc, ab, agg, br, br_idx, point, reps, eligible, a, out, tau_sources, tau_in_sources)
    except Exception as e:  # never lose the tables because of a plotting problem
        print(f"[WARN] figures failed: {e}")
    print(f"[DONE] {G} genes used; results in {out}/")


def per_triplet(ext, nodes, br_idx, a, out, use_genes):
    from scipy.stats import mannwhitneyu
    files = sorted(ext.glob("triplets_*.tsv.gz"))
    if not files:
        return
    d = nodes[nodes.run == "data"].set_index(["gene", "node"])[a.age_stat] * a.time_unit_ma
    p = nodes[nodes.run == "prior"].set_index(["gene", "node"])[a.age_stat] * a.time_unit_ma
    for f in files:
        t = pd.read_csv(f, sep="\t")
        t = t[t.gene.isin(use_genes)]
        b = t.branch.iloc[0]
        ki = pd.MultiIndex.from_arrays([t.gene, t.node_in]); ko = pd.MultiIndex.from_arrays([t.gene, t.node_out])
        t["tin_data"] = d.reindex(ki).to_numpy(); t["tout_data"] = d.reindex(ko).to_numpy()
        t["tin_prior"] = p.reindex(ki).to_numpy(); t["tout_prior"] = p.reindex(ko).to_numpy()
        t["r_data"] = t.tin_data / t.tout_data; t["r_prior"] = t.tin_prior / t.tout_prior
        classes = br_idx.loc[b, "classes"].split(",")
        if br_idx.loc[b, "kind"] == "binary":
            comps = [("D-C", ["D1", "D2"], ["C"]), ("D1-D2", ["D1"], ["D2"])]
        else:
            comps = [(f"{x} - {y}", [x], [y]) for i, x in enumerate(classes) for y in classes[i + 1:]]
        rows = []
        for trip, g in t.groupby(["tax_K1", "tax_K2", "tax_S"]):
            for name, A, B in comps:
                ga, gb = g[g.cls.isin(A)], g[g.cls.isin(B)]
                row = dict(branch=b, triplet=",".join(trip), contrast=name, nA=len(ga), nB=len(gb))
                for col in ("tin_data", "tout_data", "r_data", "tin_prior", "tout_prior", "r_prior"):
                    if len(ga) >= 10 and len(gb) >= 10:
                        row[f"{col}_dmedian"] = ga[col].median() - gb[col].median()
                        row[f"{col}_p_mwu"] = mannwhitneyu(ga[col], gb[col]).pvalue
                rows.append(row)
        pd.DataFrame(rows).to_csv(out / f"per_triplet_{b}.tsv", sep="\t", index=False, float_format="%.5g")


CLASS_COLORS = {"C": "#2166ac", "D1": "#d6604d", "D2": "#f4a582", "D": "#b2182b"}
POLY_COLORS = ["#1b9e77", "#d95f02", "#7570b3"]


def write_dates_tables(agg, br, tau_sources, tau_in_sources, out):
    """Tidy tables for plotting gene-tree dates in R (or anything else).

    dates_long.tsv.gz: one row per gene x branch x class x run (pooled D omitted;
      use cls_group to pool). Ages are the gene's class means over its
      branch-specific triplets, in Ma; w is the gene weight n(g,c)/n(g) used by
      every gene-weighted statistic, so weight histograms by w.
    branch_tau.tsv: species-tree tau_in / tau_out per branch and tau source,
      for reference lines.
    """
    bi = br.set_index("branch")
    d = agg[agg.cls != "D"]
    rows = []
    for run, suf in (("data", "d"), ("prior", "p")):
        x = d[["gene", "branch", "cls", "n_gc", "n_g", "w"]].copy()
        x["t_in_Ma"] = d[f"a_tin_{suf}"].values
        x["t_out_Ma"] = d[f"a_tout_{suf}"].values
        x["r"] = d[f"a_r_{suf}"].values
        if x[["t_in_Ma", "t_out_Ma"]].isna().all().all():
            continue
        x.insert(3, "run", run)
        rows.append(x)
    long = pd.concat(rows, ignore_index=True)
    long.insert(2, "branch_label", long.branch.map(bi.label))
    long.insert(3, "T_CU", long.branch.map(bi.T_CU))
    long.insert(4, "kind", long.branch.map(bi.kind))
    long.insert(6, "cls_group", np.where((long.kind == "binary") & long.cls.isin(["D1", "D2"]), "D", long.cls))
    long = long.rename(columns={"n_gc": "n_triplets", "n_g": "n_triplets_gene_branch"})
    order = {b: i for i, b in enumerate(br.branch)}
    long = long.sort_values(["branch", "run", "cls", "gene"], key=lambda c: c.map(order) if c.name == "branch" else c)
    long.to_csv(out / "dates_long.tsv.gz", sep="\t", index=False, float_format="%.6g")
    tr = []
    for src, m in tau_sources.items():
        for b in br.branch:
            tr.append(dict(branch=b, branch_label=bi.loc[b, "label"], tau_source=src,
                           tau_in_Ma=tau_in_sources.get(src, {}).get(b, np.nan), tau_out_Ma=m.get(b, np.nan)))
    pd.DataFrame(tr, columns=["branch", "branch_label", "tau_source", "tau_in_Ma", "tau_out_Ma"]).to_csv(
        out / "branch_tau.tsv", sep="\t", index=False, float_format="%.6g")
    (out / "plot_branch_distributions.R").write_text(R_SCRIPT)


def branch_distributions(agg, br_idx, branches, tau_sources, tau_in_sources, out):
    """One figure per branch: gene-weighted histograms of class t_in, t_out and r
    (data run filled, prior-only run as dashed outlines), class weighted means,
    and species-tree tau lines. PNG per branch + one multi-page PDF."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    odir = out / "branch_distributions"
    odir.mkdir(exist_ok=True)
    src = next((s for s in tau_sources if s != "gene_lower5"), None) or next(iter(tau_sources), None)
    with PdfPages(out / "branch_distributions.pdf") as pdf:
        for b in branches:
            if b not in br_idx.index:
                continue
            sub = agg[agg.branch == b]
            classes = br_idx.loc[b, "classes"].split(",")
            colors = (CLASS_COLORS if br_idx.loc[b, "kind"] == "binary"
                      else dict(zip(classes, POLY_COLORS)))
            fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
            for ax, q, xl in zip(axes, ("tin", "tout", "r"), ("t_in (Ma)", "t_out (Ma)", "r = t_in / t_out")):
                vals = pd.concat([sub[f"a_{q}_d"], sub[f"a_{q}_p"]]).to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if len(vals) < 5:
                    ax.set_visible(False); continue
                lo, hi = np.percentile(vals, [0.5, 99.5])
                bins = np.linspace(lo, hi, 45)
                for c in classes:
                    s = sub[sub.cls == c]
                    for ver, style in (("d", "filled"), ("p", "dashed")):
                        x = s[f"a_{q}_{ver}"].to_numpy(float); w = s.w.to_numpy(float)
                        m = np.isfinite(x)
                        if m.sum() < 5:
                            continue
                        if style == "filled":
                            ax.hist(x[m], bins=bins, weights=w[m], density=True, histtype="stepfilled",
                                    alpha=0.35, color=colors.get(c), label=f"{c}  F={w[m].sum():.0f}")
                            ax.hist(x[m], bins=bins, weights=w[m], density=True, histtype="step",
                                    lw=1.2, color=colors.get(c))
                            ax.axvline(np.average(x[m], weights=w[m]), color=colors.get(c), lw=1, ls=":")
                        else:
                            ax.hist(x[m], bins=bins, weights=w[m], density=True, histtype="step",
                                    lw=0.9, ls="--", color=colors.get(c),
                                    label=f"{c} prior" if q == "tin" else None)
                if src and q in ("tin", "tout"):
                    tau = (tau_in_sources if q == "tin" else tau_sources).get(src, {}).get(b, np.nan)
                    if np.isfinite(tau):
                        ax.axvline(tau, color="k", lw=1.2, label=f"τ ({src})" if q == "tin" else None)
                ax.set_xlabel(xl); ax.set_yticks([])
            axes[0].legend(fontsize=7, frameon=False)
            T = br_idx.loc[b, "T_CU"]
            fig.suptitle(f"{b}: {br_idx.loc[b, 'label']}   (T = {T:.3f} CU)   "
                         f"gene-weighted; dotted = weighted mean; dashed outline = prior only", fontsize=9)
            fig.tight_layout()
            fig.savefig(odir / f"{b}.png", dpi=170)
            pdf.savefig(fig)
            plt.close(fig)


R_SCRIPT = r"""# Gene-tree date distributions per species-tree branch, from ilsdate / aeon ils-analyze.
# Usage:  Rscript plot_branch_distributions.R [results_dir]   (default: this script's folder)
# Needs: ggplot2 (and readr is not required).
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args)) args[1] else "."
d    <- read.delim(gzfile(file.path(dir, "dates_long.tsv.gz")), stringsAsFactors = FALSE)
tau  <- read.delim(file.path(dir, "branch_tau.tsv"), stringsAsFactors = FALSE)

d$facet   <- paste0(d$branch, ": ", d$branch_label, "  (T=", signif(d$T_CU, 3), ")")
d$facet   <- factor(d$facet, levels = unique(d$facet))          # species-tree (preorder) order
tau$facet <- levels(d$facet)[match(tau$branch, d$branch[match(levels(d$facet), d$facet)])]
tau       <- tau[!is.na(tau$facet) & tau$tau_source != "gene_lower5", ]

plot_q <- function(var, xlab, tau_col = NULL, run = "data") {
  x <- d[d$run == run & is.finite(d[[var]]), ]
  p <- ggplot(x, aes(x = .data[[var]])) +
    geom_histogram(aes(y = after_stat(density), weight = w, fill = cls, colour = cls),
                   bins = 40, position = "identity", alpha = 0.35) +
    facet_wrap(~ facet, scales = "free", ncol = 3) +
    labs(x = xlab, y = "gene-weighted density", fill = "class", colour = "class") +
    theme_bw(base_size = 8) + theme(strip.text = element_text(size = 6))
  if (!is.null(tau_col) && nrow(tau))
    p <- p + geom_vline(data = tau, aes(xintercept = .data[[tau_col]], linetype = tau_source))
  p
}

pdf(file.path(dir, "branch_distributions_R.pdf"), width = 11, height = 8.5)
print(plot_q("t_in_Ma",  "t_in (Ma)",  "tau_in_Ma"))
print(plot_q("t_out_Ma", "t_out (Ma)", "tau_out_Ma"))
print(plot_q("r",        "r = t_in / t_out"))
if ("prior" %in% d$run) {
  print(plot_q("t_in_Ma",  "t_in (Ma), prior only",  run = "prior"))
  print(plot_q("t_out_Ma", "t_out (Ma), prior only", run = "prior"))
}
invisible(dev.off())
message("wrote ", file.path(dir, "branch_distributions_R.pdf"))

# One branch, interactively:
#   b <- subset(d, branch == "b05" & run == "data")
#   ggplot(b, aes(t_out_Ma, weight = w, fill = cls)) + geom_histogram(bins = 40, alpha = .5, position = "identity")
# Pool D1 and D2 with fill = cls_group.
"""


def figures(bc, ab, agg, br, br_idx, point, reps, eligible, a, out, tau_sources, tau_in_sources):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    binb = br[br.kind == "binary"].sort_values("T_CU")
    # Fig 1: forest of DiD (D - C), with data and prior
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * len(binb) + 1.8), sharey=True)
    for ax, q in zip(axes, ("tin", "tout")):
        s = bc[(bc.contrast == "D-C") & (bc.quantity == q)].set_index("branch").loc[binb.branch]
        y = np.arange(len(s))
        for ver, off, mk, col in (("did", 0, "o", "k"), ("data", 0.22, "^", "tab:blue"),
                                  ("prior", -0.22, "v", "tab:orange")):
            ax.errorbar(s[ver], y + off, xerr=[s[ver] - s[f"{ver}_lo"], s[f"{ver}_hi"] - s[ver]],
                        fmt=mk, color=col, ms=4, lw=1, label=ver if q == "tin" else None,
                        alpha=1 if ver == "did" else 0.6)
        if a.two_ne_g:
            ax.scatter(s.pred_CU * a.two_ne_g, y, marker="|", s=150, color="red",
                       label=f"predicted × {a.two_ne_g:g} Myr/CU" if q == "tin" else None)
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_title(f"Δ{'t_in' if q == 'tin' else 't_out'} (discordant − concordant), Ma")
        dim = ~s.reportable.to_numpy()
        for yy in y[dim]:
            ax.axhspan(yy - 0.5, yy + 0.5, color="0.92", zorder=-1)
    axes[0].set_yticks(np.arange(len(binb)))
    axes[0].set_yticklabels([f"{r.label}  (T={r.T_CU:.3f})" for r in binb.itertuples()], fontsize=7)
    axes[0].legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(out / "fig1_forest.png", dpi=200); plt.close(fig)

    # Fig 2: observed vs predicted
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ver = "did" if bc["did"].notna().any() else "data"
    for q, mk, pk in (("tin", "o", "pred_dTin_CU"), ("tout", "s", "pred_dTout_CU")):
        s = bc[(bc.contrast == "D-C") & (bc.quantity == q) & bc.branch.isin(eligible)]
        ax.errorbar(s.pred_CU, s[ver], yerr=[s[ver] - s[f"{ver}_lo"], s[f"{ver}_hi"] - s[ver]], fmt=mk,
                    ms=5, lw=0.8, label=f"{q} ({ver})")
        kS, kI = ("reg", q, ver, "slope"), ("reg", q, ver, "intercept")
        if kS in point:
            xx = np.linspace(0, max(s.pred_CU.max(), 0.01), 10)
            ax.plot(xx, point[kI][0] + point[kS][0] * xx, lw=1,
                    label=f"{q}: slope {point[kS][0]:.2f} Myr/CU")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlabel("predicted Δ (CU)")
    ax.set_ylabel("observed Δ, data − prior (Ma)" if ver == "did" else "observed Δ, data run (Ma)")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(out / "fig2_obs_vs_pred.png", dpi=200)
    plt.close(fig)

    # Fig 3: per-branch distributions of gene-tree dates (all branches unless --plot-branches)
    sel = [x for x in a.plot_branches.split(",") if x] or list(br.branch)
    branch_distributions(agg, br_idx, sel, tau_sources, tau_in_sources, out)

    # Fig 4: absolute test
    if len(ab):
        src = next(iter(tau_sources))
        s = ab[(ab.tau_source == src) & (ab.method == "pt") & ab.cls.isin(["C", "D"])
               & ab.branch.isin(br.branch[br.kind == "binary"])]
        fig, ax = plt.subplots(figsize=(5.5, 4.2))
        TT = np.linspace(0.001, max(2.0, np.nanmax(br.T_CU) * 1.05), 200)
        ax.plot(TT, [msc_expectations(t)["frac_C_older"] for t in TT], color="tab:blue", lw=1,
                label="ILS: concordant")
        ax.axhline(1, color="tab:red", lw=1, label="ILS: discordant")
        for c, col in (("C", "tab:blue"), ("D", "tab:red")):
            ss = s[s.cls == c]
            ax.errorbar(ss.T_CU, ss.frac, yerr=[ss.frac - ss.lo, ss.hi - ss.frac], fmt="o", color=col,
                        ms=4, lw=0.8, label=f"observed {c}")
        ax.set_xlabel("T (CU)"); ax.set_ylabel(f"fraction with t_in > τ_out ({src})")
        ax.set_ylim(-0.02, 1.05); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(out / "fig4_absolute.png", dpi=200); plt.close(fig)



def register(sub) -> None:
    p = sub.add_parser('ils-analyze', help='Gene-weighted ILS dating statistics, bootstrap, tables and plots', description=__doc__,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(p)
    p.set_defaults(func=run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return run(ap.parse_args(argv))
