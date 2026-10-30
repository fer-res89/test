"""Gene-tree topology classes and branch lengths per species-tree branch, and a
test of whether discordance depends on a gene's substitution rate.

Under ILS a gene's topology does not depend on how fast it evolves, so the
discordant share should be flat across rate classes. Excess discordance in
slow (low-signal) genes is the signature of gene-tree estimation error (GTEE).

Input: gene trees with branch lengths in substitutions/site (IQ-TREE
.treefile or any Newick), one file per gene; the species tree used for the
ILS analyses. Gene trees are rooted on --outgroup (genes where the outgroup is
missing are skipped).

For every gene and every testable species-tree branch (same branch-specific
triplets as `aeon ils-extract`), each triplet is classified as C, D1 or D2.
Per gene x branch the script records the triplet shares, the gene's majority
class, the mean internal branch length (path between the triplet's inner and
outer nodes, in subs/site), that length relative to the gene's total tree
length, and the total tree length (rate proxy; alignments of similar length).

With --true-tree-dir (simulations) the true class of every triplet is also
recorded, so estimation error is measured directly: the share of inferred
discordant triplets that are truly concordant, by rate class.

Outputs (--out-dir):
  gene_branch_lengths.tsv   gene x branch rows (the amniote all_branch_lengths.tsv analogue)
  discordance_vs_rate.tsv   per branch: discordant share by rate tercile, chi2 and Spearman tests,
                            internal-branch lengths by class (and GTEE rates with true trees)
  discordance_vs_rate.png   discordant share by rate tercile for every branch
"""
from __future__ import annotations

import argparse
import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .ils_common import (OutgroupError, define_branches, msc_expectations, parse_outgroup,
                         read_clade_names, read_polytomy, read_tree_file, reroot_on_outgroup)


def _gene_id(path: Path, suffix: str, rx):
    stem = path.name[: -len(suffix)] if path.name.endswith(suffix) else path.stem
    if rx:
        m = rx.search(stem)
        if not m:
            return None
        return m.group(1) if m.groups() else m.group(0)
    return stem


def _matrices(tree, tax_index, ntax):
    """LCA node id, LCA depth (edges) and root distance (subs) for tip pairs."""
    A = np.full((ntax, ntax), -1, dtype=np.int32)
    D = np.full((ntax, ntax), -1, dtype=np.int32)
    depth, rdist, ids = {tree: 0}, {tree: 0.0}, {}
    for nd in tree.preorder():
        for c in nd.children:
            depth[c] = depth[nd] + 1
            rdist[c] = rdist[nd] + (c.length or 0.0)
    leafsets = {}
    k = 0
    for nd in tree.postorder():
        if nd.is_leaf():
            leafsets[nd] = [tax_index[nd.name]]
            continue
        ids[k] = nd
        ls = [leafsets[c] for c in nd.children]
        for x, y in itertools.combinations(range(len(ls)), 2):
            ix, iy = np.array(ls[x]), np.array(ls[y])
            A[np.ix_(ix, iy)] = k; A[np.ix_(iy, ix)] = k
            D[np.ix_(ix, iy)] = depth[nd]; D[np.ix_(iy, ix)] = depth[nd]
        leafsets[nd] = sum(ls, [])
        k += 1
    rd = np.array([rdist[ids[i]] for i in range(k)])
    return A, D, rd


def classify(tree, br_idx, tax_index, ntax, present):
    """Per branch: (cls array, inner node, outer node, ok mask, valid mask)."""
    A, D, rd = _matrices(tree, tax_index, ntax)
    out = {}
    for bid, (a, b, s) in br_idx.items():
        v = present[a] & present[b] & present[s]
        if not v.any():
            continue
        a_, b_, s_ = a[v], b[v], s[v]
        dd = np.stack([D[a_, b_], D[a_, s_], D[b_, s_]], 1)
        aa = np.stack([A[a_, b_], A[a_, s_], A[b_, s_]], 1)
        cls = dd.argmax(1)
        srt = np.sort(dd, 1)
        ok = srt[:, 2] > srt[:, 1]           # resolved triplets only
        rows = np.arange(len(cls))
        nin, nout = aa[rows, cls], aa[rows, (cls + 1) % 3]
        out[bid] = (cls, nin, nout, ok, v, rd)
    return out


def tree_length(tree):
    return sum((n.length or 0.0) for n in tree.preorder() if n.parent is not None)


def add_arguments(ap):
    ap.add_argument("--species-tree", required=True, help="species tree used for the ILS analyses")
    ap.add_argument("--gene-tree-dir", required=True, help="directory of gene trees with branch lengths")
    ap.add_argument("--suffix", default=".treefile", help="gene tree file suffix (default .treefile)")
    ap.add_argument("--true-tree-dir", help="true gene topologies (ils-simulate trees/), for measuring GTEE")
    ap.add_argument("--true-suffix", default=".tre")
    ap.add_argument("--gene-id-regex", help="regex on the file stem; group 1 (or whole match) is the gene id")
    ap.add_argument("--outgroup", default="Gallus_gallus",
                    help="outgroup taxon or comma-separated outgroup clade")
    ap.add_argument("--polytomy", help="TSV group<TAB>taxa (three rows), as in ils-extract")
    ap.add_argument("--clade-names", help="TSV name<TAB>taxa for readable branch labels")
    ap.add_argument("--bins", type=int, default=3, help="number of rate classes (default terciles)")
    ap.add_argument("--rate", choices=["total", "mean"], default="total",
                    help="rate proxy: total tree length (default) or mean branch length (if taxon "
                         "occupancy varies a lot between genes)")
    ap.add_argument("--out-dir", default="genetree_rate_check")


def run(a):
    from scipy.stats import chi2_contingency, spearmanr
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    sp = read_tree_file(a.species_tree)
    taxa = sorted(sp.leaf_names())
    og = parse_outgroup(a.outgroup)
    if not og <= set(taxa):
        sys.exit(f"outgroup taxa not in species tree: {sorted(og - set(taxa))}")
    branches, warns = define_branches(sp, og, set(taxa) - og, read_polytomy(a.polytomy),
                                      read_clade_names(a.clade_names))
    for w in warns:
        print("[WARN]", w)
    tax_index = {t: i for i, t in enumerate(taxa)}
    ntax = len(taxa)
    bmap = {b.bid: b for b in branches}
    br_idx = {}
    for br in branches:
        g = [np.array([tax_index[t] for t in grp]) for grp in br.groups]
        m = np.array(np.meshgrid(*g, indexing="ij")).reshape(3, -1)
        br_idx[br.bid] = (m[0], m[1], m[2])

    rx = re.compile(a.gene_id_regex) if a.gene_id_regex else None
    files = sorted(Path(a.gene_tree_dir).glob(f"*{a.suffix}"))
    true_files = {}
    if a.true_tree_dir:
        for f in Path(a.true_tree_dir).glob(f"*{a.true_suffix}"):
            gid = _gene_id(f, a.true_suffix, rx)
            if gid:
                true_files[gid] = f
    if not files:
        sys.exit(f"[ERROR] no files matching *{a.suffix} in {a.gene_tree_dir}")
    print(f"[INFO] {len(files)} gene-tree files (*{a.suffix}); e.g. {files[0].name}")
    if a.true_tree_dir:
        print(f"[INFO] {len(true_files)} true trees with an ID (*{a.true_suffix} in {a.true_tree_dir})")
    rows, skipped = [], []
    n_true_matched = 0
    for f in files:
        gid = _gene_id(f, a.suffix, rx)
        if gid is None:
            skipped.append((f.name, "gene_id_regex_no_match")); continue
        try:
            t = read_tree_file(f)
        except Exception as e:
            skipped.append((gid, f"parse:{e}")); continue
        tips = t.leaf_names()
        if not (og & tips):
            skipped.append((gid, "no_outgroup")); continue
        if tips - set(taxa):
            skipped.append((gid, f"taxa_not_in_species_tree:{','.join(sorted(tips - set(taxa))[:3])}")); continue
        try:
            t = reroot_on_outgroup(t, og)
        except OutgroupError as e:
            skipped.append((gid, f"outgroup:{e}")); continue
        present = np.zeros(ntax, bool); present[[tax_index[x] for x in tips]] = True
        L = tree_length(t)
        nbr = sum(1 for n in t.preorder() if n.parent is not None)
        rate = L if a.rate == "total" else L / nbr
        res = classify(t, br_idx, tax_index, ntax, present)
        tres = None
        if gid in true_files:
            n_true_matched += 1
            tt = reroot_on_outgroup(read_tree_file(true_files[gid]), og)
            tres = classify(tt, br_idx, tax_index, ntax, present)
        for bid, (cls, nin, nout, ok, v, rd) in res.items():
            br = bmap[bid]
            n = int(ok.sum())
            if n == 0:
                continue
            c, i_, o_ = cls[ok], nin[ok], nout[ok]
            internal = rd[i_] - rd[o_]                 # path length between inner and outer node
            shares = [(c == k).mean() for k in range(3)]
            maj = int(np.argmax(shares)) if sorted(shares)[-1] > sorted(shares)[-2] else -1
            row = dict(gene=gid, branch=bid, label=br.label(), kind=br.kind, T_CU=br.T, n_triplets=n,
                       share_C=shares[0], share_D1=shares[1], share_D2=shares[2],
                       topology_class=(br.classes[maj] if maj >= 0 else "tie"),
                       internal_branch_length=float(internal.mean()),
                       internal_branch_length_majority=(float(internal[c == maj].mean()) if maj >= 0 else np.nan),
                       gene_tree_total_length=L, rate_proxy=rate, n_taxa=len(tips))
            row["internal_rel"] = row["internal_branch_length_majority"] / L if L > 0 else np.nan
            if tres is not None and bid in tres:
                tc, _, _, tok, tv, _ = tres[bid]
                # align triplets: both masks come from the same presence vector
                both = ok & tok
                tcls = tc[both]; icls = cls[both]
                row["true_share_C"] = (tcls == 0).mean() if len(tcls) else np.nan
                row["frac_triplets_wrong"] = (tcls != icls).mean() if len(tcls) else np.nan
                if br.kind == "binary":
                    dmask = icls != 0
                    row["frac_inferredD_trueC"] = ((tcls[dmask] == 0).mean() if dmask.any() else np.nan)
            rows.append(row)
    if skipped:
        pd.DataFrame(skipped, columns=["gene", "reason"]).to_csv(out / "skipped_genes.tsv", sep="\t", index=False)
        reasons = pd.Series([r.split(":")[0] for _, r in skipped]).value_counts()
        print("[INFO] skipped:", ", ".join(f"{k} {v}" for k, v in reasons.items()))
        ex = [(g, r) for g, r in skipped[:3]]
        print("[INFO] e.g.", "; ".join(f"{g}: {r}" for g, r in ex))
    if a.true_tree_dir:
        print(f"[INFO] {n_true_matched} gene trees matched to a true tree")
        if n_true_matched == 0:
            print("[WARN] no true trees matched: check that --gene-id-regex extracts the same ID from both "
                  "file names, and --true-suffix")
    if not rows:
        sys.exit("[ERROR] no gene trees could be classified (see the skipped reasons above)")
    gb = pd.DataFrame(rows)
    gb.to_csv(out / "gene_branch_lengths.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"[INFO] {gb.gene.nunique()} gene trees classified; {len(skipped)} skipped")

    # ---- discordance vs rate, per branch
    labels = ["slow", "mid", "fast"] if a.bins == 3 else [f"q{i + 1}" for i in range(a.bins)]
    summ = []
    for bid, g in gb[gb.topology_class != "tie"].groupby("branch", sort=False):
        br = bmap[bid]
        g = g.copy()
        if br.kind == "binary":
            g["disc"] = (g.topology_class != "C").astype(int)
        else:
            # polytomy: share of the least frequent-overall class pair is not meaningful; use "not majority"
            top = g.topology_class.value_counts().idxmax()
            g["disc"] = (g.topology_class != top).astype(int)
        try:
            g["rbin"] = pd.qcut(g.rate_proxy, a.bins, labels=labels)
        except ValueError:
            continue
        ct = pd.crosstab(g.rbin, g.disc)
        p_chi = chi2_contingency(ct)[1] if ct.shape == (a.bins, 2) else np.nan
        rho, p_rho = spearmanr(g.rate_proxy, g.disc)
        row = dict(branch=bid, label=br.label(), kind=br.kind, T_CU=br.T, n_genes=len(g),
                   disc_share=g.disc.mean(),
                   msc_disc_share=(1 - msc_expectations(br.T)["P_C"]) if br.kind == "binary" else 2 / 3,
                   chi2_p=p_chi, spearman_rho=rho, spearman_p=p_rho)
        for lab, s in g.groupby("rbin", observed=True):
            row[f"disc_{lab}"] = s.disc.mean()
            row[f"internal_rel_D_{lab}"] = s.loc[s.disc == 1, "internal_rel"].median()
            row[f"internal_rel_C_{lab}"] = s.loc[s.disc == 0, "internal_rel"].median()
            if "frac_inferredD_trueC" in s:
                row[f"GTEE_{lab}"] = s.loc[s.disc == 1, "frac_inferredD_trueC"].mean()
        row["median_internal_C"] = g.loc[g.disc == 0, "internal_branch_length_majority"].median()
        row["median_internal_D"] = g.loc[g.disc == 1, "internal_branch_length_majority"].median()
        row["median_rate_C"] = g.loc[g.disc == 0, "rate_proxy"].median()
        row["median_rate_D"] = g.loc[g.disc == 1, "rate_proxy"].median()
        if "true_share_C" in g:
            row["true_disc_share"] = 1 - g.true_share_C.mean()
            row["GTEE_overall"] = g.loc[g.disc == 1, "frac_inferredD_trueC"].mean() if br.kind == "binary" else np.nan
        summ.append(row)
    sm = pd.DataFrame(summ)
    if len(sm):
        from .ils_analyze import bh
        sm["chi2_q_bh"] = bh(sm.chi2_p.to_numpy())
    sm.to_csv(out / "discordance_vs_rate.tsv", sep="\t", index=False, float_format="%.5g")
    cols = ["branch", "T_CU", "n_genes", "msc_disc_share", "disc_share"] + \
           [c for c in sm.columns if c.startswith("disc_") and c != "disc_share"] + ["chi2_p", "chi2_q_bh"] + \
           [c for c in ["true_disc_share", "GTEE_overall"] if c in sm.columns]
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(sm.sort_values("T_CU")[cols].round(3).to_string(index=False))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        s = sm.sort_values("T_CU")
        fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(s) + 2), 4))
        x = np.arange(len(s))
        w = 0.8 / a.bins
        for i, lab in enumerate(labels):
            ax.bar(x + (i - (a.bins - 1) / 2) * w, s[f"disc_{lab}"], w, label=f"{lab} genes")
        ax.scatter(x, s.msc_disc_share, marker="_", s=200, color="k", label="MSC expectation (T)", zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{b}\nT={t:.2f}" + ("*" if q < 0.05 else "") for b, t, q in
                            zip(s.branch, s.T_CU, s.chi2_q_bh)], fontsize=6)
        ax.set_ylabel("discordant share of genes")
        ax.set_title("Discordance by gene rate class (* BH q < 0.05: discordance depends on rate)", fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(out / "discordance_vs_rate.png", dpi=200)
    except Exception as e:
        print(f"[WARN] plot failed: {e}")
    print(f"[DONE] wrote {out}/")


def register(sub) -> None:
    p = sub.add_parser("ils-genetrees", help="Gene-tree topology classes, branch lengths and discordance-vs-rate check",
                       description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(p)
    p.set_defaults(func=run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    return run(ap.parse_args(argv))
