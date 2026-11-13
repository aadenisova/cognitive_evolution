# cognitive_evolution

Code for the analyses of the manuscript on avian technical innovation: phylogenetic analyses, conserved non-coding
elements (CNEEs) with RERconverge, branch-site selection tests (HyPhy), and cell-type-specific differentially
accessible regions (DARs) from chicken telencephalon scATAC-seq.

- `src/phylogenetic_analysis/` - phylogenetic signal, PGLS, RRphylo
- `src/RERConverge_analysis/` - RER of CNEEs vs technical innovation, permutations, DAR enrichment (see its README)
- `src/dar-analysis/` - pseudobulk DESeq2 pipeline for DARs (see its README)
- `src/branch_site_selections/` - alignment filtering and HyPhy tests
- `src/conservative_regions/` - earlier CNEE / liftover scripts
- `src/data_rerconverge/` - summary tables and figure scripts (Fig. 3: `src_plot/figure_article_panel.R`)

Input data (RER trees, phenotype table, DAR BED files) are not stored in this repository; scripts use the original
absolute paths documented in `src/RERConverge_analysis/README.md`.

Changes in v3.3 are listed in `CHANGES.md`.
