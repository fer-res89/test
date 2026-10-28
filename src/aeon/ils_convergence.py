"""Section 5.6: compare two independent chains for the replicate subset.

For every gene with outputs in both --run1 and --run2, compare the posterior
median (from mcmc.txt; mean from FigTree.tre if no samples) of every internal
node. Writes

  convergence_nodes.tsv   gene x node: age in each chain, difference, relative difference
  convergence_genes.tsv   gene: max |relative difference|, max |difference| / mean HPD width
  convergence_fail.txt    genes exceeding the thresholds (feed to --exclude-genes)
  convergence_1to1.png    run1 vs run2 node ages (points should lie on the 1:1 line)

Node identity is checked by clade, so the comparison does not depend on node
numbering. ESS per node is computed in `aeon ils-extract` (gene_qc.tsv).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .ils_extract import RunError, discover, read_run


def add_arguments(ap):
    ap.add_argument("--run1", required=True)
    ap.add_argument("--run2", required=True)
    ap.add_argument("--fig-suffix", default=".FigTree.tre")
    ap.add_argument("--mcmc-suffix", default=".mcmc.txt")
    ap.add_argument("--gene-id-regex")
    ap.add_argument("--time-unit-ma", type=float, default=100.0)
    ap.add_argument("--max-rel-diff", type=float, default=0.05, help="fail if any node differs by more")
    ap.add_argument("--max-hpd-frac", type=float, default=0.25,
                    help="fail if any node differs by more than this fraction of its mean 95%% HPD width")
    ap.add_argument("--out-dir", default="convergence")


def run(a):
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    m1 = discover(Path(a.run1), a.fig_suffix, a.mcmc_suffix, a.gene_id_regex)
    m2 = discover(Path(a.run2), a.fig_suffix, a.mcmc_suffix, a.gene_id_regex)
    genes = sorted(set(m1) & set(m2))
    print(f"[INFO] {len(genes)} genes with two chains")
    rows, grows = [], []
    for g in genes:
        try:
            r1 = read_run(m1[g][0], m1[g][1], 0.0, 0.05)
            r2 = read_run(m2[g][0], m2[g][1], 0.0, 0.05)
        except (RunError, Exception) as e:
            grows.append(dict(gene=g, status=f"error:{e}")); continue
        t1, i1, mean1, h1, s1, _ = r1
        t2, i2, mean2, h2, s2, _ = r2
        c1 = {frozenset(n.leaf_names()): k for k, n in enumerate(i1)}
        c2 = {frozenset(n.leaf_names()): k for k, n in enumerate(i2)}
        if set(c1) != set(c2):
            grows.append(dict(gene=g, status="topology_differs")); continue
        a1 = np.median(s1, 0) if s1 is not None else mean1
        a2 = np.median(s2, 0) if s2 is not None else mean2
        mx_rel = mx_h = 0.0
        for cl, k1 in c1.items():
            k2 = c2[cl]
            d = a1[k1] - a2[k2]
            rel = abs(d) / max(a1[k1], a2[k2])
            w = np.nanmean([h1[k1, 1] - h1[k1, 0], h2[k2, 1] - h2[k2, 0]])
            hf = abs(d) / w if np.isfinite(w) and w > 0 else np.nan
            mx_rel = max(mx_rel, rel); mx_h = np.nanmax([mx_h, hf])
            rows.append(dict(gene=g, n_desc=len(cl), age1=a1[k1] * a.time_unit_ma,
                             age2=a2[k2] * a.time_unit_ma, diff_Ma=d * a.time_unit_ma, rel_diff=rel, hpd_frac=hf))
        fail = mx_rel > a.max_rel_diff or (np.isfinite(mx_h) and mx_h > a.max_hpd_frac)
        grows.append(dict(gene=g, status="fail" if fail else "ok", max_rel_diff=mx_rel, max_hpd_frac=mx_h))
    nd = pd.DataFrame(rows); gd = pd.DataFrame(grows)
    nd.to_csv(out / "convergence_nodes.tsv", sep="\t", index=False, float_format="%.5g")
    gd.to_csv(out / "convergence_genes.tsv", sep="\t", index=False, float_format="%.5g")
    bad = gd[gd.status != "ok"].gene
    (out / "convergence_fail.txt").write_text("".join(f"{x}\n" for x in bad))
    if len(nd):
        print(f"[INFO] node ages: r = {np.corrcoef(nd.age1, nd.age2)[0, 1]:.5f}; "
              f"median |rel diff| = {nd.rel_diff.median():.4f}; 99th pct = {nd.rel_diff.quantile(0.99):.4f}")
    print(gd.status.value_counts().to_string())
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        ax.scatter(nd.age1, nd.age2, s=3, alpha=0.4)
        lim = [0, max(nd.age1.max(), nd.age2.max()) * 1.02]
        ax.plot(lim, lim, color="k", lw=0.8)
        ax.set_xlabel("chain 1 node age (Ma)"); ax.set_ylabel("chain 2 node age (Ma)")
        fig.tight_layout(); fig.savefig(out / "convergence_1to1.png", dpi=200)
    except Exception as e:
        print(f"[WARN] plot failed: {e}")



def register(sub) -> None:
    p = sub.add_parser('ils-convergence', help='Compare two MCMCTree chains per gene', description=__doc__,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(p)
    p.set_defaults(func=run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return run(ap.parse_args(argv))
