#!/usr/bin/env python3
"""
12_dar_fdr_23annotated.py
BH-FDR for the DAR enrichment tests, restricted to the 23 annotated cell types.

Two correction families are written to the output:
  *_23 : BH within each direction (23 tests per direction)
  *_46 : BH across both directions together (23 cell types x 2 directions = 46 tests)
         -> used for Figure 3B in the manuscript (4 significant cell types).

FIX (v3.3): the previous bh_fdr() divided the sorted p-values by n, n-1, ..., 1
instead of the rank 1, 2, ..., n, so the smallest p-values were not corrected at all.
Note: 'perm_pvalue' produced by 11_dar_enrichment_50pct.R is a tie-corrected normal
approximation of the same rank-sum statistic (not a permutation test), so
padj_perm_* is not independent of padj_wilcoxon_*.
"""
import numpy as np
import pandas as pd

INPUT = "/mnt/results/rerconverge_cont_v2/phase_f_dar/dar_enrichment_summary_v3.csv"
OUTPUT = "/mnt/results/rerconverge_cont_v2/dar_enrichment_23annotated_v3.csv"


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (same result as R p.adjust(method='BH'))."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)          # p_(i) * n / i
    adj = np.minimum.accumulate(adj[::-1])[::-1]       # enforce monotonicity from the top rank down
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def main():
    dar = pd.read_csv(INPUT)
    annotated = dar[~dar["cell_type"].str.contains("Unannotated")].copy().reset_index(drop=True)
    print(f"Annotated cell types: {annotated['cell_type'].nunique()}")

    # family of 46: both directions together
    annotated["padj_wilcoxon_46"] = bh_fdr(annotated["pvalue_wilcoxon"].values)
    annotated["padj_perm_46"] = bh_fdr(annotated["perm_pvalue"].values)

    # family of 23: within direction
    annotated["padj_wilcoxon_23"] = np.nan
    annotated["padj_perm_23"] = np.nan
    for direction in ["negative", "positive"]:
        idx = annotated.index[annotated["direction"] == direction]
        annotated.loc[idx, "padj_wilcoxon_23"] = bh_fdr(annotated.loc[idx, "pvalue_wilcoxon"].values)
        annotated.loc[idx, "padj_perm_23"] = bh_fdr(annotated.loc[idx, "perm_pvalue"].values)

    annotated.to_csv(OUTPUT, index=False)

    for direction in ["negative", "positive"]:
        sub = annotated[annotated["direction"] == direction]
        print(f"  {direction}: Wilcoxon FDR<0.05  BH-23 = {(sub['padj_wilcoxon_23'] < 0.05).sum()}, "
              f"BH-46 = {(sub['padj_wilcoxon_46'] < 0.05).sum()}")
    sig = annotated[(annotated["direction"] == "negative") & (annotated["padj_wilcoxon_46"] < 0.05)]
    for _, r in sig.sort_values("padj_wilcoxon_46").iterrows():
        print(f"    {r['cell_type']:12s} effect={r['effect_size']:.4f} FDR46={r['padj_wilcoxon_46']:.2e}")
    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
