# Cell-type-specific DAR identification from chicken telencephalon scATAC-seq

Reproducible pipeline for identifying differentially accessible regions (DARs) from
single-cell ATAC-seq data using pseudobulk aggregation and DESeq2.

## Overview

This pipeline processes scATAC-seq data from the chicken telencephalon atlas
(GEO: GSE262321) to produce cell-type-specific accessible chromatin regions. The
workflow consists of three stages:

1. **Pseudobulk preparation** (`prepare_pseudobulk.py`) — loads the sparse count
   matrix, matches barcodes to metadata, filters peaks, and aggregates counts per
   (cell type × sample) group.
2. **Differential accessibility testing** (`run_deseq2_da.R`) — runs DESeq2 with a
   one-vs-all design using sum-to-zero contrasts on 60 pseudobulk samples.
3. **Output generation** (`build_outputs.py`) — builds a binary accessibility
   matrix, per-cell-type BED files, and a summary table.

## Data requirements

Download the following files from GEO accession GSE262321 and place them in a
working directory (default: `/workspace`):

| File | Description |
|---|---|
| `atac_counts_matrix.mtx.gz` | Sparse count matrix (features × cells, Matrix Market) |
| `atac_features.tsv.gz` | Peak coordinates (3 columns: peak\_id, peak\_id2, type) |
| `atac_barcodes.tsv.gz` | Cell barcodes (1 column: cell\_id) |
| `meta_data.tsv.gz` | Cell metadata (cell\_id, celltype, sample\_name, description\_celltype) |

Reference genome: GRCg6a (GalGal6, UCSC assembly).

## Dependencies

### Python (≥ 3.11)

```
scipy >= 1.11
pandas >= 2.0
numpy
```

### R (≥ 4.3)

```
DESeq2 >= 1.42
BiocParallel >= 1.36
data.table
```

### System

```
bedtools >= 2.31
```

## Usage

Set the working directory and number of R workers via environment variables:

```bash
export WORKSPACE=/path/to/data
export RESULTS=/path/to/output
export N_WORKERS=8
```

Run the three scripts in order:

```bash
python prepare_pseudobulk.py
Rscript run_deseq2_da.R
python build_outputs.py
```

## Output files

| File | Format | Description |
|---|---|---|
| `accessible_regions_binary_matrix.csv` | CSV | Binary matrix (182,550 peaks × 30 cell types) |
| `da_results_full.tsv` | TSV | Full DESeq2 results for all 5,476,500 tests |
| `da_summary.csv` | CSV | Per-cell-type summary (cell count, peaks tested, DAR count) |
| `bed_files/` | BED | 30 BED files (one per cell type), sorted by significance |

## Repository structure

```
dar-analysis/
├── README.md
├── .gitignore
├── code/
│   ├── prepare_pseudobulk.py    # Stage 1: pseudobulk aggregation
│   ├── run_deseq2_da.R          # Stage 2: DESeq2 DA analysis
│   └── build_outputs.py         # Stage 3: output generation
└── methods/
    └── dar_methods.md           # Nature-style methods description
```

## Key parameters

| Parameter | Value |
|---|---|
| Initial peaks | 487,770 |
| Peaks after filtering | 182,550 |
| Annotated cells | 23,271 |
| Cell types | 30 (23 annotated + 7 unannotated) |
| Biological replicates | 2 |
| Pseudobulk groups | 60 (30 × 2) |
| Minimum group size | 30 cells |
| Size factor method | poscounts |
| Design formula | `~ sample + cell_type` |
| Contrast type | contr.sum (sum-to-zero, one-vs-all) |
| FDR method | Benjamini-Hochberg |
| Significance threshold | padj < 0.05 AND log2FC > 1 |
| Total DARs | 801,499 |

## License

MIT
