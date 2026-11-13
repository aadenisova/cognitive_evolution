# Changes in v3.3

- `src/RERConverge_analysis/12_dar_fdr_23annotated.py`: fixed Benjamini-Hochberg function (divisor was n..1 instead of 1..n);
  now writes BH within direction (`*_23`) and across both directions (`*_46`). DAR enrichment: 4 significant cell types
  (BH across 46 tests; OL, AST, OPC, EXC_GLU-8) instead of 7.
- `src/RERConverge_analysis/10_dar_overlap_50pct.R`: fixed `rownames(rownames(old_mat))`.
- `src/RERConverge_analysis/11_dar_enrichment_50pct.R`: documented that `perm_*` columns are a normal approximation, not a permutation test.
- `src/RERConverge_analysis/01_*.R`, `02_*.R`: lambda computed with `lower.tail = FALSE` and annotated as a top-5,000 statistic;
  new `02b_genomewide_lambda.R` computes the genome-wide lambda.
- `src/data_rerconverge/src_plot/figure_article_panel.R`: no longer needs ggtree (ape + ggplot2 only); tree with aligned tips
  (`TREE_STYLE=aligned`, old look: `depth`); PDF written with cairo_pdf (embedded fonts); paths configurable via `FIG_ROOT` / `FIG_OUT`.
- `src/data_rerconverge/rerconverge_log_noexon_summary_v3.xlsx` (Supplementary Table S2): corrected BH columns, GO sheets removed,
  sheets `Cell_Clusters` and `Transform_Calibration` added. Removed stray files (Excel lock file, `.DS_Store`, `Rplots.pdf`).
- READMEs updated.
