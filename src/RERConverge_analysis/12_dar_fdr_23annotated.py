#!/usr/bin/env python3
"""
05_dar_recompute_23annotated_v3.py
Recompute BH FDR on 23 annotated cell types only for v3 DAR enrichment.
Directions: negative and positive only (omnibus 'all' removed in v3.1).
"""
import pandas as pd
import numpy as np

INPUT = "/mnt/results/rerconverge_cont_v2/phase_f_dar/dar_enrichment_summary_v3.csv"
OUTPUT = "/mnt/results/rerconverge_cont_v2/dar_enrichment_23annotated_v3.csv"

def bh_fdr(pvals, n_tests):
    pvals = np.asarray(pvals, dtype=float)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    sorted_padj = sorted_pvals * n_tests / np.arange(n_tests, 0, -1)
    for i in range(n_tests - 2, -1, -1):
        sorted_padj[i] = min(sorted_padj[i], sorted_padj[i + 1])
    sorted_padj = np.clip(sorted_padj, 0, 1.0)
    padj = np.ones(n_tests)
    for i, idx in enumerate(sorted_idx):
        padj[idx] = sorted_padj[i]
    return padj

def main():
    dar = pd.read_csv(INPUT)
    annotated_cts = [ct for ct in dar["cell_type"].unique() if "Unannotated" not in ct]
    print(f"Annotated cell types: {len(annotated_cts)}")

    rows = []
    for direction in ["negative", "positive"]:
        sub = dar[dar["direction"] == direction].copy()
        sub_ann = sub[sub["cell_type"].isin(annotated_cts)].copy().reset_index(drop=True)
        n = len(sub_ann)
        sub_ann["padj_wilcoxon_23"] = bh_fdr(sub_ann["pvalue_wilcoxon"].values, n)
        sub_ann["padj_perm_23"] = bh_fdr(sub_ann["perm_pvalue"].values, n)
        rows.append(sub_ann)

    result = pd.concat(rows, ignore_index=True)
    result.to_csv(OUTPUT, index=False)

    for direction in ["negative", "positive"]:
        sub = result[result["direction"] == direction]
        sig_w = (sub["padj_wilcoxon_23"] < 0.05).sum()
        sig_p = (sub["padj_perm_23"] < 0.05).sum()
        print(f"  {direction}: Wilcoxon sig={sig_w}, Perm sig={sig_p}")
        if direction == "negative":
            sig = sub[sub["padj_wilcoxon_23"] < 0.05].sort_values("padj_wilcoxon_23")
            for _, r in sig.iterrows():
                tag = "GREEN" if r["effect_size"] >= 0.015 else "YELLOW"
                print(f"    {tag}: {r['cell_type']:15s} effect={r['effect_size']:.4f} "
                      f"wilcoxon_fdr={r['padj_wilcoxon_23']:.6f} perm_fdr={r['padj_perm_23']:.6f}")

    print(f"\nSaved: {OUTPUT}")

if __name__ == "__main__":
    main()
