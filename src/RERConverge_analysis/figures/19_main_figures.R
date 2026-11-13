#!/usr/bin/env Rscript
# 09_figures.R
# Generate all individual figures for v3.2.
# A: Phylogenetic tree + innovation values (rectangular)
# B1: Top 10 intron CNEEs — Rho, distance category, p-value, gene names
# B2: Promoter ≤2kb vs >2kb — |Rho| comparison
# D: DAR enrichment — effect sizes for 23 cell types (negative direction)
# S: Element enrichment distribution (supplement)
#
# All figures in R using ggplot2/ggtree. Saved as PNG + SVG.

# Workaround: ggtree 3.14.0 calls is.waive() which was removed from ggplot2 4.0.3
is.waive <- function(x) inherits(x, "waiver")

suppressPackageStartupMessages({
  library(ggtree)
  library(ggplot2)
  library(ape)
  library(treeio)
  library(tidytree)
  library(ggprism)
  library(patchwork)
})

BASE <- "/mnt/results/rerconverge_cont_v2"
OUTPUT_DIR <- "/mnt/results/figures"
DATA_DIR <- "/mnt/results/plotting/data"

# Helper: save both PNG and SVG
save_fig <- function(plot, name, width=10, height=8) {
  png_path <- file.path(OUTPUT_DIR, paste0(name, ".png"))
  svg_path <- file.path(OUTPUT_DIR, paste0(name, ".svg"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(svg_path, plot, width = width, height = height, bg = "white")
  message(sprintf("  Saved: %s + %s", png_path, svg_path))
}

# Color palette
COL_POS <- "#E94057"   # red for positive Rho
COL_NEG <- "#0279EE"   # blue for negative Rho
COL_GREEN <- "#75A025"
COL_YELLOW <- "#FF9400"
COL_GREY <- "#CCCCCC"

# ============================================================
# Figure A: Phylogenetic tree + innovation
# ============================================================
message("=== Figure A: Phylogenetic tree + innovation ===")

# Load tree
tree_rooted <- read.tree(file.path(DATA_DIR, "species_tree_rooted.nwk"))

# Load innovation table
innov <- read.csv(file.path(DATA_DIR, "innovation_table.csv"), stringsAsFactors = FALSE)

# Match innovation values to tree tip order
tip_order <- tree_rooted$tip.label
innov_matched <- innov[match(gsub(" ", "_", innov$scientific_name), tip_order), ]
innov_matched$label <- tip_order

# Tree plot
p_tree <- ggtree(tree_rooted, layout = "rectangular") +
  geom_tiplab(aes(label = label), size = 2.8, font = "italic", align = TRUE, linesize = 0.3) +
  theme_tree2() +
  xlim(0, max(node.depth.edgelength(tree_rooted)) * 1.4)

# Innovation bar chart (same y-order as tree tips)
innov_bar <- ggplot(innov_matched, aes(x = TECHINNO2025_ResEff, y = label)) +
  geom_col(aes(fill = TECHINNO2025_ResEff), width = 0.6) +
  scale_fill_viridis_c(name = "Tech.\nInnovation", option = "C", guide = "none") +
  scale_y_discrete(limits = rev(tip_order), labels = rep("", length(tip_order))) +
  labs(x = "Tech. Innovation") +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 5, 0, 0))

# Combine with patchwork
p_a <- p_tree + innov_bar +
  plot_layout(widths = c(4, 1)) &
  theme(plot.margin = margin(5, 5, 5, 5))

save_fig(p_a, "figure_A_tree_innovation", width = 12, height = 14)
message("  Figure A done")

# ============================================================
# Figure B1: Top 10 intron CNEEs
# ============================================================
message("=== Figure B1: Top 10 intron CNEEs ===")

robust <- read.csv(file.path(BASE, "robust_cnees_annotated_v3.csv"), stringsAsFactors = FALSE)
sym_map <- read.csv(file.path(BASE, "gene_symbol_map.csv"), stringsAsFactors = FALSE)

# Add gene symbols
robust$geneSymbol <- sym_map$symbol[match(as.integer(robust$geneId), sym_map$geneId)]

# Top 10 intron CNEEs by permpval
introns <- robust[robust$category == "Intron", ]
top10_intron <- introns[order(introns$permpval), ][1:10, ]

# Distance category
top10_intron$dist_abs <- abs(top10_intron$distanceToTSS)
top10_intron$dist_cat <- cut(top10_intron$dist_abs,
                              breaks = c(0, 1e4, 5e4, 1e5, Inf),
                              labels = c("≤10kb", "10-50kb", "50-100kb", ">100kb"),
                              include.lowest = TRUE)

# Label: gene symbol + CNEE ID
top10_intron$label <- paste0(top10_intron$geneSymbol, " (", top10_intron$CNEE, ")")
top10_intron$label <- factor(top10_intron$label, levels = top10_intron$label[order(top10_intron$Rho)])

# p-value annotation
top10_intron$pval_label <- sprintf("P=%.1e", top10_intron$P)

p_b1 <- ggplot(top10_intron, aes(x = Rho, y = label, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = pval_label), hjust = ifelse(top10_intron$Rho > 0, -0.1, 1.1),
            size = 3, color = "black") +
  facet_wrap(~ dist_cat, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_fill_manual(values = c("positive" = COL_POS, "negative" = COL_NEG),
                    name = "Direction") +
  labs(x = "Rho (RER-phenotype correlation)",
       y = "Gene (CNEE ID)",
       title = NULL) +
  theme_prism(base_size = 11, base_family = "sans") +
  theme(axis.text.y = element_text(size = 8, face = "italic"),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom")

save_fig(p_b1, "figure_B1_intron_top10", width = 10, height = 8)
message("  Figure B1 done")

# ============================================================
# Figure B2: Promoter ≤2kb vs >2kb comparison
# ============================================================
message("=== Figure B2: Promoter ≤2kb vs >2kb ===")

prom <- robust[robust$category == "Promoter", ]
prom$dist_group <- ifelse(abs(prom$distanceToTSS) <= 2000, "≤2kb", ">2kb")
prom$dist_group <- factor(prom$dist_group, levels = c("≤2kb", ">2kb"))
prom$abs_Rho <- abs(prom$Rho)

# Wilcoxon test
wt <- wilcox.test(prom$abs_Rho ~ prom$dist_group)
pval_str <- sprintf("Wilcoxon P = %.3f", wt$p.value)

p_b2 <- ggplot(prom, aes(x = dist_group, y = abs_Rho, fill = dist_group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, aes(color = dist_group)) +
  annotate("text", x = 1.5, y = max(prom$abs_Rho) * 1.05,
           label = pval_str, size = 3.5, hjust = 0.5) +
  scale_fill_manual(values = c("≤2kb" = "#0279EE", ">2kb" = "#FD9BED"), guide = "none") +
  scale_color_manual(values = c("≤2kb" = "#0279EE", ">2kb" = "#FD9BED"), guide = "none") +
  labs(x = "Distance from TSS",
       y = "|Rho| (absolute RER-phenotype correlation)",
       title = NULL) +
  theme_prism(base_size = 11, base_family = "sans") +
  theme(legend.position = "none")

save_fig(p_b2, "figure_B2_promoter_comparison", width = 6, height = 6)
message("  Figure B2 done")

# ============================================================
# Figure D: DAR enrichment
# ============================================================
message("=== Figure D: DAR enrichment ===")

dar <- read.csv(file.path(BASE, "dar_enrichment_23annotated_v3.csv"), stringsAsFactors = FALSE)
dar_neg <- dar[dar$direction == "negative", ]
dar_neg <- dar_neg[order(dar_neg$effect_size, decreasing = TRUE), ]

# Significance category
dar_neg$sig_cat <- ifelse(dar_neg$padj_wilcoxon_23 < 0.05 & dar_neg$effect_size >= 0.015, "Significant (FDR<0.05)",
                   ifelse(dar_neg$padj_wilcoxon_23 < 0.05, "FDR<0.05 (weak effect)", "Not significant"))
dar_neg$sig_cat <- factor(dar_neg$sig_cat, levels = c("Significant (FDR<0.05)", "FDR<0.05 (weak effect)", "Not significant"))

# Significance asterisks
dar_neg$star <- ifelse(dar_neg$padj_wilcoxon_23 < 0.001, "***",
                ifelse(dar_neg$padj_wilcoxon_23 < 0.01, "**",
                ifelse(dar_neg$padj_wilcoxon_23 < 0.05, "*", "")))

dar_neg$cell_type <- factor(dar_neg$cell_type, levels = dar_neg$cell_type[order(dar_neg$effect_size)])

p_d <- ggplot(dar_neg, aes(x = effect_size, y = cell_type, fill = sig_cat)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = star), hjust = -0.3, size = 4, color = "black") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0.015, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c("Significant (FDR<0.05)" = COL_GREEN,
                                "FDR<0.05 (weak effect)" = COL_YELLOW,
                                "Not significant" = COL_GREY),
                    name = NULL) +
  labs(x = "Effect size (rank-biserial correlation)",
       y = "Cell type",
       title = NULL) +
  theme_prism(base_size = 11, base_family = "sans") +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 9)) +
  xlim(min(dar_neg$effect_size) * 1.3, max(dar_neg$effect_size) * 1.25)

save_fig(p_d, "figure_D_dar_enrichment", width = 9, height = 8)
message("  Figure D done")

# ============================================================
# Supplement: Element enrichment distribution
# ============================================================
message("=== Supplement: Element enrichment ===")

fisher <- read.csv(file.path(BASE, "phase_e_elements/element_fisher_summary_v2.csv"), stringsAsFactors = FALSE)

# Use positive and negative directions only
fisher_sub <- fisher[fisher$direction %in% c("positive", "negative"), ]
fisher_sub$direction <- factor(fisher_sub$direction, levels = c("positive", "negative"))

# Significance marker
fisher_sub$sig <- ifelse(fisher_sub$padj < 0.05, "*", "")

# Order categories
cat_order <- c("Distal Intergenic", "Intron", "Promoter", "Downstream", "3'UTR", "5'UTR")
fisher_sub$category <- factor(fisher_sub$category, levels = cat_order)

p_s <- ggplot(fisher_sub, aes(x = category, y = enrichment_ratio, fill = direction)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(aes(label = sig, group = direction),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("positive" = COL_POS, "negative" = COL_NEG),
                    name = "Direction",
                    labels = c("Positive Rho", "Negative Rho")) +
  labs(x = "Element category",
       y = "Enrichment ratio (robust / genome proportion)",
       title = NULL) +
  theme_prism(base_size = 11, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "bottom")

save_fig(p_s, "figure_S_element_enrichment", width = 9, height = 6)
message("  Supplement done")

message("\n=== Step 09 complete ===")
