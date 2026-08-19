library(cowplot)
library(pdftools)
library(magick)

# launch /src/Phylogenetic_tests/phylo_innovation_analysis_new.R
ROOT_vgp <-"/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests//results/TOTALINNOVATIONS2025_ResEff/"
ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/cognitive_evolution/"

ROOT_plots <- "results_500/TECHINNO2025_ResEff/figures/"

BASE <- file.path(ROOT, "results", PHENOTYPE)

OUT <- file.path(ROOT_vgp, "figures", "contMap_order_vgp.pdf")
OUT_bars <- file.path(ROOT, ROOT_plots, "phylogenetic_signal_bar.pdf")
OUT_boxs <- file.path(ROOT, ROOT_plots, "foreground_background_boxplot.pdf")
  

img <- image_read_pdf(OUT, density = 600)
tree <- ggdraw() +
  draw_image(img)
tree

img <- image_read_pdf(OUT_bars, density = 600)
signal_bars <- ggdraw() +
  draw_image(img)
signal_bars

img <- image_read_pdf(OUT_boxs, density = 600)
signal_boxs <- ggdraw() +
  draw_image(img)
signal_boxs

# launch /src/Phylogenetic_tests/phylo_innovation_analysis_new.R
#bottom <- plot_grid(
#  plots$signal_bars,
#  plots$foreground_background_box,
#  ncol = 2,
#  labels = c("B", "C")
#)

bottom <- plot_grid(
  signal_bars,
  foreground_background_box,
  ncol = 2,
  labels = c("B", "C")
)

bottom
final <- plot_grid(
  tree,
  bottom,
  ncol = 1,
  rel_heights = c(2, 1),   # верхняя и нижняя половины одинаковой высоты
  labels = c("A", "")
)

final
ggsave(
  "final_figure_1.2.pdf",
  final,
  width = 12,
  height = 10,
  units = "in"
)

ggsave(
  "/Users/aleksandradenisova/Desktop/paper_figures/final_figure_1.2.png",
  final,
  width = 12,
  height = 10,
  units = "in"
)

