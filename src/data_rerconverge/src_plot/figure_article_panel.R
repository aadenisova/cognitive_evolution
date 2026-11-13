#!/usr/bin/env Rscript
# =============================================================================
# figure_article_panel.R
# Единственный актуальный скрипт статьи (A + B + C).
#
#   A  tree + technical innovation rate (top left)
#   B  DAR enrichment, negative (top) / positive (bottom) (top right)
#   C  10 promoter correlations: 5 positive (top row), 5 negative (bottom),
#      closest to TSS (≤ TSS_MAX_BP); left = nearest
#
# Правки — блок CONFIG.
# Запуск:  Rscript src_plot/figure_article_panel.R
#
# Устаревшие 09_figures.R / 10_assemble_figures.R / 02_plot_promoter_rer_correlation.R
# только перенаправляют сюда.
# =============================================================================

is.waive <- function(x) inherits(x, "waiver")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
  library(patchwork)
  library(cowplot)
  library(readxl)
  library(scales)
})

# =============================================================================
# CONFIG
# =============================================================================

ROOT <- Sys.getenv("FIG_ROOT", "/Users/aleksandradenisova/Desktop/data_rerconverge")
TREE_STYLE   <- Sys.getenv("TREE_STYLE", "aligned")   # "aligned" (tips on one line, default) or "depth" (old ggtree-like look)
TREE_XEXPAND <- as.numeric(Sys.getenv("TREE_XEXPAND", "1.30"))  # right-hand room for tip labels (like hexpand(0.30))
OUT  <- Sys.getenv("FIG_OUT", file.path(ROOT, "figures"))

TREE_NWK     <- file.path(ROOT, "species_tree_rooted.nwk")
INNOV_CSV    <- file.path(ROOT, "innovation_table.csv")
NAMES_CSV    <- file.path(ROOT, "species_common_names.csv")
RER_CSV      <- file.path(ROOT, "promoter_rer_data_40.csv")
XLSX         <- file.path(ROOT, "rerconverge_log_noexon_summary_v3.xlsx")

# NULL = авто: 5 positive + 5 negative, ближайшие к TSS (в пределах TSS_MAX_BP)
# Или задать 10 CNEE вручную (сначала 5 pos, затем 5 neg):
# TOP_CNEES <- c(...)
TOP_CNEES  <- NULL
N_PER_DIR  <- 5
TSS_MAX_BP <- 1000   # отбор среди элементов ≤ этого расстояния до TSS

# innovation: red = low, blue = high (как в палитре Rho)
COL_POS    <- "#E94057"
COL_NEG    <- "#0279EE"
COL_SIG    <- "#3D7A22"
COL_WEAK   <- "#D48B16"
COL_NS     <- "#C8C8C8"

DAR_EFFECT_CUTOFF <- 0.015
DAR_FDR_COL <- "padj_wilcoxon_46"   # BH across all 46 tests (23 cell types x 2 directions); column added in *_corrected.xlsx

FIG_W <- 10.2
FIG_H <- 12.2
DPI   <- 600

# высоты: верх (A+B) / низ (C) — C ниже
ROW_HEIGHTS <- c(1.38, 0.88)
# ширины верха: дерево / DAR
TOP_WIDTHS  <- c(1.40, 1.00)

FONT <- "Helvetica"

# =============================================================================
# helpers
# =============================================================================

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# plotmath for facet header: gene, TSS, rho, P
facet_pm <- function(gene, dist, rho, p) {
  dist_s <- sprintf("%+d bp", as.integer(dist))
  rho_s  <- sprintf("%.2f", rho)
  if (isTRUE(p < 0.001)) {
    e <- floor(log10(p))
    m <- p / 10^e
    p_s <- paste0(sprintf("%.1f", m), " %*% 10^{", e, "}")
  } else {
    p_s <- sprintf("%.3f", p)
  }
  paste0("atop(italic('", gene, "')~'(", dist_s, ")', rho==", rho_s,
         "*','~italic(P)==", p_s, ")")
}

theme_fig <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = FONT) %+replace%
    theme(
      axis.line  = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text  = element_text(colour = "black", size = base_size - 0.5),
      axis.title = element_text(colour = "black", size = base_size, face = "plain"),
      plot.title = element_text(face = "bold", size = base_size + 0.4, hjust = 0,
                                margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.8, colour = "grey20",
                                   hjust = 0, margin = margin(b = 3)),
      plot.margin = margin(3, 5, 3, 3),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.4, face = "plain", lineheight = 1.02,
                                margin = margin(b = 1.5))
    )
}

save_both <- function(plot, name, width, height) {
  png_path <- file.path(OUT, paste0(name, ".png"))
  pdf_path <- file.path(OUT, paste0(name, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = DPI, bg = "white",
         device = ragg::agg_png, scaling = 1)
  ggsave(pdf_path, plot, width = width, height = height, bg = "white",
         device = cairo_pdf, family = FONT)  # cairo_pdf embeds fonts (journal-friendly)
  message(sprintf("  %s", png_path))
  message(sprintf("  %s", pdf_path))
}

fill_innov <- function() {
  scale_fill_gradientn(
    colours = c("#C8102E", "#E94057", "#F7C1C8", "#B7D4F5", "#0279EE", "#003B99"),
    values  = scales::rescale(c(-0.60, -0.25, 0, 0.5, 1.4, 2.60)),
    guide   = "none"
  )
}

# =============================================================================
# data
# =============================================================================

message("Loading data...")
tree  <- read.tree(TREE_NWK)
innov <- read.csv(INNOV_CSV, stringsAsFactors = FALSE)
nams  <- read.csv(NAMES_CSV, stringsAsFactors = FALSE)
rer   <- read.csv(RER_CSV, stringsAsFactors = FALSE)
rob   <- as.data.frame(read_excel(XLSX, sheet = "Robust_CNEEs"))
dar   <- as.data.frame(read_excel(XLSX, sheet = "DAR_Enrichment"))

innov$tip <- gsub(" ", "_", innov$scientific_name)
nams$tip  <- gsub(" ", "_", nams$scientific_name)
stopifnot(all(tree$tip.label %in% innov$tip))

missing_names <- setdiff(tree$tip.label, nams$tip)
if (length(missing_names)) {
  stop("No common name for: ", paste(gsub("_", " ", missing_names), collapse = ", "))
}
name_map <- setNames(nams$common_name, nams$tip)

# --- promoter CNEEs (must have per-species RER) -----------------------------
prom <- rob[rob$category == "Promoter", ]
prom <- prom[prom$CNEE %in% unique(rer$CNEE), ]
prom$abs_dist <- abs(prom$distanceToTSS)

if (is.null(TOP_CNEES)) {
  pick_dir <- function(d, dir, n) {
    x <- d[d$direction == dir, ]
    x <- x[x$abs_dist <= TSS_MAX_BP, ]
    if (nrow(x) < n) {
      stop(sprintf("Only %d %s promoter CNEEs within %d bp of TSS (need %d)",
                   nrow(x), dir, TSS_MAX_BP, n))
    }
    x <- x[order(x$abs_dist, x$P, x$permpval), ]
    head(x, n)
  }
  pos <- pick_dir(prom, "positive", N_PER_DIR)
  neg <- pick_dir(prom, "negative", N_PER_DIR)
  TOP_CNEES <- c(pos$CNEE, neg$CNEE)
}

sel <- prom[match(TOP_CNEES, prom$CNEE), ]
if (any(is.na(sel$CNEE))) {
  stop("CNEE not found: ", paste(TOP_CNEES[is.na(sel$CNEE)], collapse = ", "))
}

sel$gene_show <- ifelse(is.na(sel$geneSymbol) | sel$geneSymbol == "",
                        sel$CNEE, sel$geneSymbol)
sel$panel_lab <- mapply(facet_pm, sel$gene_show, sel$distanceToTSS, sel$Rho, sel$P,
                        USE.NAMES = FALSE)
sel$panel_lab <- factor(sel$panel_lab, levels = sel$panel_lab)

write.csv(
  sel[, c("CNEE", "geneSymbol", "direction", "Rho", "P", "permpval",
          "distanceToTSS", "abs_dist")],
  file.path(OUT, "selected_promoters.csv"),
  row.names = FALSE
)
message("Selected promoter CNEEs:")
print(sel[, c("CNEE", "geneSymbol", "direction", "Rho", "P", "permpval", "distanceToTSS")],
      row.names = FALSE)

dfb <- rer[rer$CNEE %in% TOP_CNEES, ]
dfb$panel_lab <- sel$panel_lab[match(dfb$CNEE, sel$CNEE)]
dfb$direction <- sel$direction[match(dfb$CNEE, sel$CNEE)]

# =============================================================================
# A. tree + technical innovation rate
# =============================================================================

message("Panel A: tree...")

# --- tree layout without ggtree (ape + ggplot2 only; same cladogram as ggtree(branch.length = "none")) ---
# TREE_STYLE = "depth"   : x = number of edges from the root (as before; tips ragged, dotted guides to labels)
# TREE_STYLE = "aligned" : x = distance from the tips in nodes (all tips on the right edge; no dotted guides)
tree_l <- ape::ladderize(tree, right = FALSE)
n_tip <- ape::Ntip(tree_l); n_all <- n_tip + ape::Nnode(tree_l)
ed <- tree_l$edge
kids <- split(ed[, 2], ed[, 1])
root <- n_tip + 1L
tips_in_order <- ed[ed[, 2] <= n_tip, 2]
yy <- numeric(n_all); yy[tips_in_order] <- seq_along(tips_in_order)
set_y <- function(nd) { if (nd <= n_tip) return(yy[nd]); v <- vapply(kids[[as.character(nd)]], set_y, numeric(1)); yy[nd] <<- mean(v); yy[nd] }
invisible(set_y(root))
dep <- numeric(n_all)
set_d <- function(nd, d) { dep[nd] <<- d; if (nd > n_tip) for (k in kids[[as.character(nd)]]) set_d(k, d + 1) }
set_d(root, 0)
hgt <- numeric(n_all)
set_h <- function(nd) { if (nd <= n_tip) return(0); h <- max(vapply(kids[[as.character(nd)]], set_h, numeric(1))) + 1; hgt[nd] <<- h; h }
invisible(set_h(root))
xx <- if (TREE_STYLE == "aligned") hgt[root] - hgt else dep
x_tip_max <- max(xx[1:n_tip])
seg_h <- data.frame(x0 = xx[ed[, 1]], x1 = xx[ed[, 2]], y = yy[ed[, 2]])
int_nodes <- as.integer(names(kids))
seg_v <- data.frame(x = xx[int_nodes],
                    y0 = vapply(int_nodes, function(nd) min(yy[kids[[as.character(nd)]]]), numeric(1)),
                    y1 = vapply(int_nodes, function(nd) max(yy[kids[[as.character(nd)]]]), numeric(1)))
tip_df <- data.frame(label = tree_l$tip.label, x = xx[1:n_tip], y = yy[1:n_tip],
                     display = name_map[tree_l$tip.label], stringsAsFactors = FALSE)
tip_df <- tip_df[order(tip_df$y), ]
tip_bottom_to_top <- tip_df$label

p_tree <- ggplot() +
  geom_segment(data = seg_h, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.28, colour = "grey25", lineend = "square") +
  geom_segment(data = seg_v, aes(x = x, xend = x, y = y0, yend = y1), linewidth = 0.28, colour = "grey25", lineend = "square") +
  { if (TREE_STYLE == "depth")
      geom_segment(data = tip_df, aes(x = x, xend = x_tip_max + 0.15, y = y, yend = y),
                   linewidth = 0.10, linetype = "dotted", colour = "grey40") } +
  geom_text(data = tip_df, aes(x = x_tip_max + 0.15, y = y, label = display),
            hjust = 0, size = 1.70, colour = "grey10", family = FONT) +
  scale_x_continuous(limits = c(0, x_tip_max * TREE_XEXPAND), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = 0.6)) +
  labs(x = " ") +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 8, colour = NA),
    plot.margin  = margin(4, 0, 2, 2),
    plot.background = element_rect(fill = NA, colour = NA)
  )

innov$label <- factor(innov$tip, levels = tip_bottom_to_top)
innov_bar <- ggplot(innov, aes(x = TECHINNO2025_ResEff, y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_col(aes(fill = TECHINNO2025_ResEff), width = 0.72) +
  fill_innov() +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(breaks = c(-0.5, 0, 1, 2),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Technical innovation rate") +
  theme_fig(8) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x  = element_text(size = 6.5),
    panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey90"),
    plot.margin = margin(4, 6, 2, -8)
  )

p_a <- p_tree + innov_bar + plot_layout(widths = c(3.15, 1.00))

# --- DAR --------------------------------------------------------------------
dar$fdr <- dar[[DAR_FDR_COL]]
dar$sig_cat <- ifelse(dar$fdr < 0.05 & dar$effect_size >= DAR_EFFECT_CUTOFF,
                      "FDR < 0.05",
               ifelse(dar$fdr < 0.05, "FDR < 0.05 (weak)", "n.s."))
dar$sig_cat <- factor(dar$sig_cat,
                      levels = c("FDR < 0.05", "FDR < 0.05 (weak)", "n.s."))
dar$star <- ifelse(dar$fdr < 0.001, "***",
            ifelse(dar$fdr < 0.01,  "**",
            ifelse(dar$fdr < 0.05,  "*", "")))
dar$dir_lab <- factor(dar$direction,
                      levels = c("negative", "positive"),
                      labels = c("Negative Rho", "Positive Rho"))

dar_neg_ord <- dar[dar$direction == "negative", ]
dar_neg_ord <- dar_neg_ord[order(dar_neg_ord$effect_size), ]
dar$cell_type <- factor(dar$cell_type, levels = dar_neg_ord$cell_type)

x_dar <- range(dar$effect_size, na.rm = TRUE)
xlim_dar <- c(min(0, x_dar[1]) - 0.003, max(0, x_dar[2]) + 0.010)

# =============================================================================
# C. 10 promoter correlations (top: positive, bottom: negative; each by |dist to TSS|)
# =============================================================================

message("Panel C: promoter correlations...")

p_c <- ggplot(dfb, aes(x = RER, y = TECHINNO2025_ResEff)) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x,
              colour = "grey25", fill = "grey80",
              linewidth = 0.4, alpha = 0.35) +
  geom_point(aes(colour = direction), size = 0.75, alpha = 0.72, stroke = 0) +
  facet_wrap(~ panel_lab, nrow = 2, ncol = 5, scales = "free_x",
             labeller = label_parsed) +
  scale_colour_manual(values = c(positive = COL_POS, negative = COL_NEG),
                      guide = "none") +
  labs(x = "Relative evolutionary rate (RER)",
       y = "Technical innovation rate") +
  theme_fig(8.4) +
  theme(
    panel.spacing = unit(0.22, "lines"),
    axis.text  = element_text(size = 6.6),
    axis.title = element_text(size = 8.4),
    strip.text = element_text(size = 7.6, lineheight = 1.05),
    plot.margin = margin(2, 4, 2, 2)
  )

# =============================================================================
# B. DAR enrichment (negative top, positive bottom)
# =============================================================================

message("Panel B: DAR enrichment...")

p_b <- ggplot(dar, aes(x = effect_size, y = cell_type, fill = sig_cat)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey20") +
  geom_vline(xintercept = DAR_EFFECT_CUTOFF, linetype = "dashed",
             linewidth = 0.25, colour = "grey55") +
  geom_col(width = 0.78) +
  geom_text(aes(label = star), hjust = -0.12, size = 2.2, colour = "black",
            na.rm = TRUE) +
  facet_wrap(~ dir_lab, ncol = 1) +
  scale_fill_manual(
    values = c("FDR < 0.05" = COL_SIG,
               "FDR < 0.05 (weak)" = COL_WEAK,
               "n.s." = COL_NS),
    drop = TRUE,
    name = NULL
  ) +
  coord_cartesian(xlim = xlim_dar, clip = "off") +
  labs(
    x = "Effect size",
    y = NULL,
    title = "DAR enrichment",
    subtitle = "Wilcoxon, BH FDR across 46 tests"
  ) +
  theme_fig(7.2) +
  theme(
    axis.text.y = element_text(size = 6.0),
    axis.text.x = element_text(size = 5.6),
    axis.title.x = element_text(size = 7),
    strip.text = element_text(size = 7.2, face = "bold"),
    panel.spacing.y = unit(0.12, "lines"),
    legend.position = "bottom",
    legend.margin = margin(t = 0, b = 0),
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 6.0),
    plot.margin = margin(0, 6, 0, 0),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 6.2, margin = margin(b = 2))
  )

# =============================================================================
# assemble
# =============================================================================

message("Assembling...")

p_ab <- cowplot::plot_grid(
  p_a, p_b, nrow = 1, rel_widths = TOP_WIDTHS,
  labels = c("A", "B"), label_size = 12, label_fontface = "bold",
  label_x = 0.002, label_y = 0.998
)
p_c_lab <- cowplot::plot_grid(
  p_c, labels = "C", label_size = 12, label_fontface = "bold",
  label_x = 0.002, label_y = 1.0
)
combined <- cowplot::plot_grid(
  p_ab, p_c_lab, ncol = 1, rel_heights = ROW_HEIGHTS
)

# =============================================================================
# export
# =============================================================================

message("Saving combined figure...")
save_both(combined, "figure_article_panel", FIG_W, FIG_H)

message("Saving individual panels...")
save_both(p_a, "panel_A_tree", width = 6.2, height = 6.4)
save_both(p_b, "panel_B_DAR", width = 4.2, height = 6.4)
save_both(p_c, "panel_C_promoter_corr", width = FIG_W, height = 4.4)

stale <- file.path(OUT, c(
  "panel_B_promoter_corr.png", "panel_B_promoter_corr.pdf",
  "panel_C_DAR.png", "panel_C_DAR.pdf"
))
unlink(stale[file.exists(stale)])

caption <- paste(
  "Figure 1. Technical innovation, DAR enrichment and promoter CNEEs.",
  "",
  "(A) Phylogeny of 92 bird species (cladogram; branch lengths not shown) with per-species technical innovation rate (TECHINNO2025 residual). Bars are aligned to tree tips and coloured from red (low) to blue (high). Species are labelled by English common names.",
  "",
  "(B) Enrichment of CNEEs overlapping cell-type differentially accessible regions (DARs) in 23 brain cell types. Wilcoxon rank-sum tests compare -log10 permutation P values of CNEEs overlapping a DAR with those of CNEEs that do not; effect size is the rank-biserial correlation. Tests were performed separately for CNEEs with negative ρ (top) and positive ρ (bottom). P values were adjusted by the Benjamini–Hochberg procedure across all 46 tests (23 cell types x 2 directions). Green, FDR < 0.05; grey, not significant. Asterisks: * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001. The dashed line marks an arbitrary effect-size reference of 0.015.",
  "",
  sprintf(
    "(C) Relative evolutionary rate (RER) versus technical innovation rate for ten promoter CNEEs lying within %g kb of the transcription start site (TSS): the five closest positive associations (top row) and the five closest negative associations (bottom row), each ordered left to right by increasing |distance to TSS|. Each point is one species. Red, positive ρ; blue, negative ρ (ρ and P are the CNEE-level RERconverge statistics: Pearson correlation across branches and parametric P). Lines show linear fits with 95%% confidence intervals. Gene symbol, signed distance to the TSS, ρ and parametric P are indicated above each panel.",
    TSS_MAX_BP / 1000
  ),
  sep = "\n"
)
caption_path <- file.path(OUT, "figure_article_panel_caption.txt")
writeLines(caption, caption_path)
message("  ", caption_path)

message("\nDone. Outputs in: ", OUT)
