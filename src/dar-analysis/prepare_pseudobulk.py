"""
Prepare pseudobulk ATAC count matrix from scATAC-seq data for DESeq2 DA analysis.

Input:
  - atac_counts_matrix.mtx.gz  (sparse matrix, features x cells)
  - atac_features.tsv.gz       (peak coordinates)
  - atac_barcodes.tsv.gz        (cell barcodes)
  - meta_data.tsv.gz            (cell type and sample annotations)

Output:
  - pseudobulk_counts.tsv       (peaks x pseudobulk groups)
  - pseudobulk_metadata.tsv     (group annotations)
  - peak_coordinates.tsv        (peak identifiers)
  - celltype_descriptions.tsv   (cell type labels)
"""
import scipy.io
import scipy.sparse
import numpy as np
import pandas as pd
import os

WORKSPACE = os.environ.get("WORKSPACE", "/workspace")

# --- Load data ---
features = pd.read_csv(
    os.path.join(WORKSPACE, "atac_features.tsv.gz"),
    sep="\t", header=None, names=["peak_id", "peak_id2", "type"]
)
barcodes = pd.read_csv(
    os.path.join(WORKSPACE, "atac_barcodes.tsv.gz"),
    sep="\t", header=None, names=["cell_id"]
)
meta = pd.read_csv(os.path.join(WORKSPACE, "meta_data.tsv.gz"), sep="\t")
mtx = scipy.io.mmread(os.path.join(WORKSPACE, "atac_counts_matrix.mtx.gz")).tocsc()
mtx.data = mtx.data.astype(np.float64)

print(f"Features: {len(features):,} peaks")
print(f"Barcodes: {len(barcodes):,} cells")
print(f"Metadata: {len(meta):,} cells, {meta['celltype'].nunique()} cell types")
print(f"Matrix: {mtx.shape} (features x cells), NNZ: {mtx.nnz:,}")

# --- Match barcodes to metadata ---
meta_indexed = meta.set_index("cell_id")
barcode_list = barcodes["cell_id"].tolist()

cell_types, samples, keep_idx = [], [], []
for i, bc in enumerate(barcode_list):
    if bc in meta_indexed.index:
        cell_types.append(meta_indexed.loc[bc, "celltype"])
        samples.append(meta_indexed.loc[bc, "sample_name"])
        keep_idx.append(i)

cell_types = np.array(cell_types)
samples = np.array(samples)
keep_idx = np.array(keep_idx)
mtx_annot = mtx[:, keep_idx]

print(f"Annotated: {len(keep_idx):,}, Unannotated: {len(barcode_list) - len(keep_idx):,}")
print(f"Cell types: {len(np.unique(cell_types))}, Samples: {len(np.unique(samples))}")

# --- Pre-filter: remove peaks with < 10 total counts ---
peak_totals = np.array(mtx_annot.sum(axis=1)).flatten()
keep_peaks = peak_totals >= 10
mtx_filtered = mtx_annot[keep_peaks, :]
peak_ids = features[keep_peaks]["peak_id"].values

print(f"Peaks after total-count filter: {keep_peaks.sum():,} (removed {(~keep_peaks).sum():,})")

# --- Pseudobulk aggregation per (cell_type, sample) ---
group_labels = np.array([f"{ct}__{s}" for ct, s in zip(cell_types, samples)])
group_counts = pd.Series(group_labels).value_counts()
valid_groups = group_counts[group_counts >= 30].index.tolist()

print(f"Pseudobulk groups: {len(valid_groups)} (>= 30 cells each)")
print(f"Group size range: {group_counts[valid_groups].min()} - {group_counts[valid_groups].max()}")

pseudobulk_data = {}
group_metadata = []
for grp in valid_groups:
    mask = group_labels == grp
    pseudobulk_data[grp] = np.array(mtx_filtered[:, mask].sum(axis=1)).flatten()
    ct, s = grp.rsplit("__", 1)
    group_metadata.append({"pseudobulk_id": grp, "cell_type": ct, "sample": s, "n_cells": int(mask.sum())})

pseudobulk_df = pd.DataFrame(pseudobulk_data, index=peak_ids)

# --- Second pre-filter: peaks non-zero in < 3 pseudobulks ---
nonzero_per_peak = (pseudobulk_df > 0).sum(axis=1)
keep_peaks2 = nonzero_per_peak >= 3
pseudobulk_df = pseudobulk_df[keep_peaks2]

print(f"Peaks after prevalence filter: {pseudobulk_df.shape[0]:,}")
print(f"Final pseudobulk matrix: {pseudobulk_df.shape}")

# --- Save outputs ---
pseudobulk_df.to_csv(os.path.join(WORKSPACE, "pseudobulk_counts.tsv"), sep="\t")
pd.DataFrame(group_metadata).to_csv(
    os.path.join(WORKSPACE, "pseudobulk_metadata.tsv"), sep="\t", index=False
)
features[keep_peaks][["peak_id"]].reset_index(drop=True).to_csv(
    os.path.join(WORKSPACE, "peak_coordinates.tsv"), sep="\t", index=False
)
meta[["celltype", "description_celltype"]].drop_duplicates().sort_values("celltype").to_csv(
    os.path.join(WORKSPACE, "celltype_descriptions.tsv"), sep="\t", index=False
)
