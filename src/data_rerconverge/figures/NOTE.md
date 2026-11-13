# Note on Figure 3 files

The files in this folder are an interim render produced with
`../src_plot/figure_article_panel.R` (ape + ggplot2, no ggtree). Panel B uses
BH FDR across all 46 tests (23 cell types x 2 directions); the statistics are
final. The tree in panel A has a different layout (tips aligned) from the
original ggtree figure in the manuscript.

`../src_plot/figure_article_panel_ggtree.R` is the original ggtree script with
the same statistical corrections (46-test FDR column, caption text). Re-running
it reproduces the original tree layout; the final journal figure should be made
with that script.
