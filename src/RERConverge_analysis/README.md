# RERconverge Pipeline: Evolutionary Rates of CNEEs vs Technological Innovation in Birds

Complete reproducible pipeline: relative evolutionary rates (RER) of 272,811 conserved
non-coding elements (CNEEs) correlated with a continuous phenotype
(TECHINNO2025_ResEff) across 92 bird species.

## Layout

```
pipeline/
├── input/            Raw data (only external inputs, nothing derived)
├── scripts/          All code, numbered in execution order (00–20)
├── output/           Result files, grouped by the script that produced them
├── intermediates/    Large regenerable intermediates + environment snapshot
└── README.md         This file
```

## Inputs (`input/`)

| File | Size | Description |
|------|------|-------------|
| `rer_trees_birds.txt` | 906 MB | Per-CNEE phylogenetic trees (92 species) with branch lengths estimated from nucleotide alignments. Read by `readTrees()` in script 01. |
| `df_merged_for_rer_no_na.csv` | 7.5 KB | Phenotype table: 92 species (rows = GCF/GCA accessions), columns include `ScientificName` and `TECHINNO2025_ResEff` (technological innovation / research efficiency residual). |
| `lifted_annotations_genomic.bed` | 8.9 MB | Genomic coordinates of all CNEEs (galGal6): chr, start, end, name. |
| `dar_bed/*_accessible_regions.bed` | 62 MB | 30 DAR (differentially accessible region) BED files, one per brain cell type. |

## Environment

- R 4.4.3 (`/opt/conda/bin/R`), `R_LIBS_USER=/workspace/.Rlib`
- Exact package snapshot: `intermediates/Rlib.tar.gz`
- Key packages: RERconverge, ape, data.table, fgsea, ChIPseeker, TxDb.Ggallus.UCSC.galGal6.refGene,
  org.Gg.eg.db, org.Hs.eg.db, homologene, ggplot2, ggtree, cowplot; Python: pandas, numpy, openpyxl
- Scripts expect helper/intermediate files in `/workspace/` (see path map below)

## Scripts (`scripts/`, execution order)

| # | Script | What it does | Key outputs |
|---|--------|--------------|-------------|
| 00 | `00_parallel_core.R` | Parallelized `coreGetResiduals`/`getAllResiduals` helpers (sourced by 01) | — |
| 01 | `01_transform_selection.R` | Reads all 272,811 CNEE trees; computes RER matrices for 4 path transforms (sqrt, log, asinh, none) + parametric correlations | `treesObj_cont.rds`, `RERmat_*.rds`, `res_*.rds` (in `/workspace/`) |
| 02 | `02_transform_calibration.R` | Calibration on top-5,000 CNEEs (100 perms; λ computed on the top-ranked CNEEs only, so not a genome-wide inflation factor; `spearman_cal` was NA because all permpval = 1/101) → **log** has the lowest λ. `02b_genomewide_lambda.R` computes the genome-wide λ. | `output/01_transform_selection/` |
| 03 | `03_permutations_fdr.R` | Main run: correlations + **2,000 permutations** per CNEE, BH + global empirical FDR | `output/03_permutations_fdr/` (full results, 272,811), `res_final_cont_v2.rds` |
| 04 | `04_exon_cleanup_robust.R` | Removes 81 exonic CNEEs (liftOver artifacts); selects robust set: permpval < 0.05 **and** sign stability ≥ 80% over 92 leave-out iterations | `output/04_exon_cleanup_robust/` (272,730 + 7,350) |
| 05 | `05_leave_out_stability.R` | 92 leave-one-species-out reruns; Spearman rho + top-N overlap vs full run | `output/05_leave_out_stability/`, `leave_out_results.rds` |
| 06 | `06_cnee_annotation_gsea.R` | ChIPseeker annotation of ALL CNEEs (galGal6 refGene) + gene-level aggregation + intermediate GSEA (superseded by 07) | `output/06_cnee_annotation_gsea/` |
| 07 | `07_gsea_by_ontology.R` | **Final GSEA**: fgsea, rankings from count/binomial × pos/neg, GO split into BP/CC/MF, chicken + human orthologs | `output/07_gsea_by_ontology/` |
| 08 | `08_element_enrichment.R` | Fisher's exact test: 6 element categories × 3 directions among robust CNEEs | `output/08_element_enrichment/` |
| 09 | `09_dar_enrichment_fgsea.R` | CNEE–DAR overlap matrix (any overlap) + fgsea-based DAR enrichment | `output/09_dar_enrichment_fgsea/` |
| 10 | `10_dar_overlap_50pct.R` | CNEE–DAR overlap with strict >50% CNEE-length criterion | `output/10_12_dar_50pct/dar_overlap_matrix_50pct.csv` |
| 11 | `11_dar_enrichment_50pct.R` | Wilcoxon rank-sum: −log10(permpval) inside vs outside DARs, per cell type | `output/10_12_dar_50pct/` |
| 12 | `12_dar_fdr_23annotated.py` | BH FDR on the 23 annotated cell types: `*_23` within direction, `*_46` across both directions (used for Fig. 3B). BH function fixed in v3.3 | `output/10_12_dar_50pct/dar_enrichment_23annotated_v3.csv` |
| 13 | `13_gene_stats_filtered.py` | Gene-level stats with v3 attribution (gene body + promoter ≤3 kb + downstream ≤300 bp; distal intergenic excluded) | `output/13_gene_stats/gene_level_stats_v3.csv` |
| 14 | `14_gene_symbol_map.R` | Gene symbol map via org.Gg.eg.db (4,440/4,443 geneIds) | `output/14_symbol_map/gene_symbol_map.csv` |
| 15 | `15_export_tree.R` | Exports rooted species tree (outgroup: *Apteryx mantelli*) + innovation table | `output/15_export_tree/` |
| 16 | `16_generate_excel.py` | Assembles the 12-sheet Excel summary | `../rerconverge_log_noexon_summary_v3.xlsx` (project root) |
| 17–20 | `figures/` | Figure scripts: promoter RER extraction + faceted scatter, main figures (tree, top-CNEEs, DAR), cowplot assembly | `output/figures_data/`, `../figures/` |

### `output/figures_data/` — promoter figure raw data

| File | Rows | Content |
|------|------|---------|
| `promoter_rer_data_40.csv` | 2,441 | Long-format data for **all 40** selected promoter CNEEs (top 20 negative + top 20 positive by permpval): one row per CNEE × species with valid values. |
| `promoter_rer_data_20.csv` | 1,113 | Subset for the figure: the **20** CNEEs actually plotted (10 per direction by \|Rho\|, ≥5 per direction within 1 kb of TSS), plus `facet_order` (panel position: rows 1–2 negative, rows 3–4 positive). |
| `promoter_cnee_selection.csv` | 40 | Selection key: one row per CNEE with `rank_by_permpval`, `in_plot` (in the 20), `within_1kb` (≤1,000 bp of TSS). |

Columns (both data files): `CNEE` (element ID), `species_accession` (GCF/GCA genome assembly), `scientific_name` (Latin binomial of that species, e.g. *Passer domesticus*), `RER` (species-level relative evolutionary rate of the CNEE, mean over that species' branches), `TECHINNO2025_ResEff` (phenotype), `Rho`/`P`/`permpval` (CNEE-level statistics), `direction`, `geneId`/`geneSymbol`, `distanceToTSS`.

### What is NOT in `output/`

`output/` keeps only the final chain of results — files consumed by a downstream
script or reported in the write-ups. Superseded intermediates (v1/v2 GSEA variants,
earlier DAR criteria, v1 element enrichment, comparison tables) remain in
`../archive/rerconverge_cont_v2_intermediate/`.

## Path map (scripts use original absolute paths)

| Path in scripts | Meaning |
|-----------------|---------|
| `/workspace/` | Machine-local workdir for large intermediates (`treesObj_cont.rds`, `RERmat_log.rds`, `res_final_cont_v2.rds`, `leave_out_results.rds`) |
| `/workspace/parallel_core.R` | Copy of script 00 |
| `/mnt/user-uploads/` | Original location of `input/` files |
| `/mnt/results/rerconverge_cont_v2/phase_*/` | Original output locations (= `output/` here) |

To re-run: place inputs per the path map (or adjust the path constants at the top of
each script), restore `intermediates/Rlib.tar.gz` to `R_LIBS_USER`, then run scripts
in numeric order. Scripts 01–03 require ~64 GB RAM (full 272,811-tree objects);
scripts 04–20 run comfortably in 16–32 GB.

## Verification status

- Script 04 was **reconstructed** after the original was lost with a deleted machine.
  It was verified to reproduce all four surviving outputs exactly
  (`rerconverge_results_cont_v2_noexon.csv`, `robust_cnee_sets.csv`,
  `cnee_annotation_cont_noexon.csv`, `robust_cnees_annotated_v3.csv`; robust files
  identical as sets, sorted by permpval). Run with `--use-existing-stability` when
  `leave_out_results.rds` is unavailable.
- All other scripts are byte-identical to the versions that produced the results.

## Key results (see project root for deliverables)

- 0 CNEEs genome-wide significant (FDR < 0.05, three correction approaches); minimum
  BH-adjusted parametric p.adj = 0.083, and the top permpval hit mCE1022938
  (Rho = 0.570, near *BRMS1L*) reaches p.adj = 0.092 — borderline, but not significant
- 7,350 nominally associated CNEEs (permpval < 0.05 and sign stability >= 80%; 3,804 positive / 3,546 negative Rho), all 100% sign-stable; none is significant after multiple-testing correction
- GO/GSEA (script 07): one nominally significant term (human CC acrosomal vesicle, driven by a single gene, KIT); not used in the manuscript
- 4 significant cell types in DAR enrichment after BH across 46 tests (negative Rho: OL, AST, OPC, EXC_GLU-8); 6 with BH within direction (+ INH_VIP, IOL). Earlier versions reported 7 because of a bug in the BH function (fixed in v3.3). Effect sizes are small (rank-biserial 0.03-0.04)
- Full sorted per-CNEE table: `../all_cnee_results_v3.csv`
- Methods write-up: `../report_log_noexon_methods_v3.md`
