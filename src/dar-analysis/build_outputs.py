"""
Build final output files from DESeq2 DA results:
  1. Binary accessibility matrix (peaks x cell_types) -> CSV
  2. Per-cell-type BED files -> bed_files/
  3. Summary table -> CSV
  4. Full DESeq2 results -> TSV
"""
import pandas as pd
import numpy as np
import os

WORKSPACE = os.environ.get("WORKSPACE", "/workspace")
RESULTS = os.environ.get("RESULTS", "/mnt/results")
bed_dir = os.path.join(RESULTS, "bed_files")
os.makedirs(bed_dir, exist_ok=True)

# --- Load data ---
results = pd.read_csv(os.path.join(WORKSPACE, "da_results_full.tsv"), sep="\t")
ct_desc = pd.read_csv(os.path.join(WORKSPACE, "celltype_descriptions.tsv"), sep="\t")
pb_meta = pd.read_csv(os.path.join(WORKSPACE, "pseudobulk_metadata.tsv"), sep="\t")

desc_map = dict(zip(ct_desc["celltype"], ct_desc["description_celltype"]))
cell_counts = pb_meta.groupby("cell_type")["n_cells"].sum().to_dict()

# --- Parse peak coordinates from peak_id (chr:start-end) ---
peaks = results["peak_id"].unique()
parts = pd.DataFrame({"peak_id": peaks})["peak_id"].str.split(r"[:\-]", expand=True)
peak_coords = pd.DataFrame({
    "peak_id": peaks,
    "chr": parts[0],
    "start": parts[1].astype(int),
    "end": parts[2].astype(int)
})

# --- 1. Binary accessibility matrix ---
sig_results = results[results["significant"] == True].copy()
binary_matrix = sig_results.pivot_table(
    index="peak_id", columns="cell_type", values="significant",
    aggfunc="max", fill_value=False
).astype(int)

all_tested_peaks = results["peak_id"].unique()
missing_peaks = set(all_tested_peaks) - set(binary_matrix.index)
if missing_peaks:
    missing_df = pd.DataFrame(0, index=list(missing_peaks), columns=binary_matrix.columns)
    binary_matrix = pd.concat([binary_matrix, missing_df])

binary_matrix = binary_matrix.sort_index()
binary_matrix = binary_matrix.reindex(sorted(binary_matrix.columns), axis=1)
binary_matrix.index.name = "peak_id"

binary_matrix.to_csv(os.path.join(RESULTS, "accessible_regions_binary_matrix.csv"))
print(f"Binary matrix: {binary_matrix.shape[0]:,} peaks x {binary_matrix.shape[1]} cell types")
print(f"Total DARs: {int(binary_matrix.values.sum()):,}")

# --- 2. Per-cell-type BED files ---
for ct in sorted(sig_results["cell_type"].unique()):
    ct_results = sig_results[sig_results["cell_type"] == ct].merge(peak_coords, on="peak_id")
    ct_results = ct_results.sort_values(["padj", "log2FC"], ascending=[True, False])

    bed = ct_results[["chr", "start", "end", "peak_id", "log2FC", "padj"]].copy()
    bed["neg_log10_padj"] = -np.log10(bed["padj"].clip(lower=1e-300))

    safe_ct = ct.replace(" ", "_").replace("/", "_")
    bed_path = os.path.join(bed_dir, f"{safe_ct}_accessible_regions.bed")
    with open(bed_path, "w") as f:
        f.write(f"# Cell type: {ct}\n")
        f.write(f"# Description: {desc_map.get(ct, ct)}\n")
        f.write(f"# DA peaks: {len(bed)} (padj<0.05, log2FC>1)\n")
        f.write("# chr\tstart\tend\tpeak_id\tlog2FC\tneg_log10_padj\n")
        bed[["chr", "start", "end", "peak_id", "log2FC", "neg_log10_padj"]].to_csv(
            f, sep="\t", index=False, header=False
        )

# --- 3. Summary table ---
summary_rows = []
for ct in sorted(results["cell_type"].unique()):
    ct_results = results[results["cell_type"] == ct]
    n_da = int(ct_results["significant"].sum())
    n_tested = len(ct_results)
    summary_rows.append({
        "cell_type": ct,
        "description": desc_map.get(ct, ct),
        "n_cells": cell_counts.get(ct, 0),
        "n_peaks_tested": n_tested,
        "n_da_peaks": n_da,
        "pct_da": round(100 * n_da / n_tested, 1) if n_tested > 0 else 0,
    })

summary_df = pd.DataFrame(summary_rows).sort_values("n_da_peaks", ascending=False)
summary_df.to_csv(os.path.join(RESULTS, "da_summary.csv"), index=False)

# --- 4. Copy full results ---
import shutil
shutil.copy(os.path.join(WORKSPACE, "da_results_full.tsv"),
            os.path.join(RESULTS, "da_results_full.tsv"))

print(f"\nSummary ({len(summary_df)} cell types):")
for _, row in summary_df.iterrows():
    print(f"  {row['cell_type']:<30s} {row['n_cells']:>6d} cells  {row['n_da_peaks']:>7d} DARs  {row['pct_da']:>5.1f}%")
