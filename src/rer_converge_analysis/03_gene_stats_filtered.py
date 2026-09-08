#!/usr/bin/env python3
"""
13_gene_stats_filtered.py
Recompute gene-level stats with filtered attribution:
  Gene body (Intron, 5'UTR, 3'UTR) + Promoter + Downstream only.
  All Distal Intergenic excluded.
Also add "attributed" column to robust CNEEs.
"""
import pandas as pd
import numpy as np
from scipy.stats import binom

BASE = "/mnt/results/rerconverge_cont_v2"

ATTRIBUTED_CATEGORIES = {"Intron", "5' UTR", "3' UTR", "Promoter", "Downstream"}

def main():
    # Load annotation
    annot = pd.read_csv(f"{BASE}/phase_d_gsea/cnee_annotation_cont_noexon.csv")
    annot["category"] = annot["annotation"].str.split("(").str[0].str.strip()
    print(f"Total annotated CNEEs: {len(annot):,}")

    # Filter to attributed only
    annot_filt = annot[annot["category"].isin(ATTRIBUTED_CATEGORIES)].copy()
    print(f"Attributed CNEEs (gene body + promoter + downstream): {len(annot_filt):,}")
    print(f"Unique genes: {annot_filt['geneId'].nunique():,}")
    print()

    # Load robust CNEEs
    robust = pd.read_csv(f"{BASE}/phase_b_permutations/robust_cnee_sets.csv")
    robust_names = set(robust["CNEE"])
    robust_dir = dict(zip(robust["CNEE"], robust["direction"]))

    # Mark robust in filtered annotation
    annot_filt["is_robust"] = annot_filt["name"].isin(robust_names)
    annot_filt["robust_dir"] = annot_filt["name"].map(robust_dir)

    n_robust_attr = annot_filt["is_robust"].sum()
    n_total_attr = len(annot_filt)
    bg_rate = n_robust_attr / n_total_attr
    print(f"Robust attributed: {n_robust_attr:,} / {n_total_attr:,} = {bg_rate:.4f}")
    print()

    # Gene-level stats
    gene_stats = []
    for gene_id, df in annot_filt.groupby("geneId"):
        n_total = len(df)
        k_pos = int(((df["is_robust"]) & (df["robust_dir"] == "positive")).sum())
        k_neg = int(((df["is_robust"]) & (df["robust_dir"] == "negative")).sum())

        # Binomial test: P(X >= k | n_total, p_bg)
        # For positive: test if k_pos is significantly high
        if k_pos > 0:
            binom_p_pos = 1 - binom.cdf(k_pos - 1, n_total, bg_rate)
        else:
            binom_p_pos = 1.0

        if k_neg > 0:
            binom_p_neg = 1 - binom.cdf(k_neg - 1, n_total, bg_rate)
        else:
            binom_p_neg = 1.0

        stat_count_pos = k_pos
        stat_count_neg = k_neg
        stat_binom_pos = -np.log10(binom_p_pos) if binom_p_pos > 0 else 300
        stat_binom_neg = -np.log10(binom_p_neg) if binom_p_neg > 0 else 300

        gene_stats.append({
            "geneId": gene_id,
            "n_total": n_total,
            "k_pos": k_pos,
            "k_neg": k_neg,
            "stat_count_pos": stat_count_pos,
            "stat_count_neg": stat_count_neg,
            "binom_p_pos": binom_p_pos,
            "binom_p_neg": binom_p_neg,
            "stat_binom_pos": stat_binom_pos,
            "stat_binom_neg": stat_binom_neg,
        })

    gene_df = pd.DataFrame(gene_stats)
    print(f"Gene-level stats: {len(gene_df)} genes")
    print(f"Genes with k_pos >= 1: {(gene_df['k_pos'] >= 1).sum()}")
    print(f"Genes with k_neg >= 1: {(gene_df['k_neg'] >= 1).sum()}")
    print(f"Max k_pos: {gene_df['k_pos'].max()}, Max k_neg: {gene_df['k_neg'].max()}")
    print(f"Median n_total: {gene_df['n_total'].median()}, Max n_total: {gene_df['n_total'].max()}")

    # Merge gene symbols from comprehensive org.Gg.eg.db map (v3.2 fix)
    sym_map = pd.read_csv(f"{BASE}/gene_symbol_map.csv")
    gene_df["geneId"] = gene_df["geneId"].astype(int)
    sym_map["geneId"] = sym_map["geneId"].astype(int)
    gene_df = gene_df.merge(sym_map[["geneId", "symbol"]], on="geneId", how="left")
    gene_df = gene_df.rename(columns={"symbol": "geneSymbol"})
    n_mapped = gene_df["geneSymbol"].notna().sum()
    print(f"Genes with symbols: {n_mapped}/{len(gene_df)}")

    # Sort by stat_count_pos descending
    gene_df = gene_df.sort_values("stat_count_pos", ascending=False).reset_index(drop=True)
    gene_df.to_csv(f"{BASE}/gene_level_stats_v3.csv", index=False)
    print(f"\nSaved: {BASE}/gene_level_stats_v3.csv")

    # Print top genes
    print("\nTop 10 by stat_count_pos:")
    print(gene_df.head(10)[["geneId", "geneSymbol", "n_total", "k_pos", "k_neg", "stat_count_pos"]].to_string())
    print("\nTop 10 by stat_count_neg:")
    gene_neg = gene_df.sort_values("stat_count_neg", ascending=False)
    print(gene_neg.head(10)[["geneId", "geneSymbol", "n_total", "k_pos", "k_neg", "stat_count_neg"]].to_string())

    # --- Add "attributed" column to robust CNEEs (preserve existing enriched v3 file) ---
    v3_path = f"{BASE}/robust_cnees_annotated_v3.csv"
    robust_annot = pd.read_csv(v3_path)
    if "attributed" not in robust_annot.columns:
        robust_annot["attributed"] = robust_annot["category"].isin(ATTRIBUTED_CATEGORIES).fillna(False)
        robust_annot.to_csv(v3_path, index=False)
        print(f"\nAdded 'attributed' column to: {v3_path}")
    else:
        print(f"\n'attributed' column already present in: {v3_path}")
    print(f"  Attributed: {robust_annot['attributed'].sum()}")
    print(f"  Not attributed (distal): {(~robust_annot['attributed']).sum()}")

if __name__ == "__main__":
    main()
