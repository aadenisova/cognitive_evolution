# phase_e_elements.R
# Phase E: Element enrichment — fgsea on genomic element categories
# Direction-specific (positive/negative Rho separately), p-value-aware ranking

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

R_LIBS_USER <- "/workspace/.Rlib"
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_e_elements"

# ============================================================
# E1: Categorize CNEEs by genomic element
# ============================================================

message("=== E1: Categorizing CNEEs by genomic element ===")

# Load annotation from Phase D
anno_df <- read.csv(file.path("/mnt/results/rerconverge_cont_v2/phase_d_gsea", "cnee_annotation_cont_full.csv"),
                     stringsAsFactors = FALSE)

# Parse annotation column into categories
# ChIPseeker annotation values:
# "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)"
# "5' UTR", "3' UTR", "1st Exon", "Other Exon"
# "1st Intron", "Other Intron", "Distal Intergenic"
# "Downstream (<=300)", "Downstream (1-2kb)", "Downstream (2-3kb)"

get_category <- function(anno) {
  if (grepl("Promoter", anno)) return("Promoter")
  if (grepl("5' UTR", anno)) return("5'UTR")
  if (grepl("3' UTR", anno)) return("3'UTR")
  if (grepl("Exon", anno)) return("Exon")
  if (grepl("Intron", anno)) return("Intron")
  if (grepl("Distal Intergenic", anno)) return("Distal Intergenic")
  if (grepl("Downstream", anno)) return("Downstream")
  return("Other")
}

# Also get sub-categories for finer resolution
get_subcategory <- function(anno) {
  if (grepl("Promoter .<=1kb", anno, fixed=TRUE)) return("Promoter (<=1kb)")
  if (grepl("Promoter .1-2kb", anno, fixed=TRUE)) return("Promoter (1-2kb)")
  if (grepl("Promoter .2-3kb", anno, fixed=TRUE)) return("Promoter (2-3kb)")
  if (grepl("5' UTR", anno)) return("5'UTR")
  if (grepl("3' UTR", anno)) return("3'UTR")
  if (grepl("1st Exon", anno)) return("1st Exon")
  if (grepl("Other Exon", anno)) return("Other Exon")
  if (grepl("1st Intron", anno)) return("1st Intron")
  if (grepl("Other Intron", anno)) return("Other Intron")
  if (grepl("Distal Intergenic", anno)) return("Distal Intergenic")
  if (grepl("Downstream", anno)) return("Downstream")
  return("Other")
}

anno_df$category <- sapply(anno_df$annotation, get_category)
anno_df$subcategory <- sapply(anno_df$annotation, get_subcategory)

message("Category distribution:")
print(table(anno_df$category))

message("\nSub-category distribution:")
print(table(anno_df$subcategory))

# Create element sets (by category)
categories <- unique(anno_df$category)
element_sets <- lapply(categories, function(cat) {
  anno_df$name[anno_df$category == cat]
})
names(element_sets) <- categories

# Create sub-category sets
subcategories <- unique(anno_df$subcategory)
element_subsets <- lapply(subcategories, function(sc) {
  anno_df$name[anno_df$subcategory == sc]
})
names(element_subsets) <- subcategories

# ============================================================
# E2: Build direction-specific rankings
# ============================================================

message("=== E2: Building direction-specific rankings ===")

# Load Phase B results
results_df <- readRDS("/workspace/res_final_cont_v2.rds")

# Positive Rho ranking (accelerated in innovative species)
pos_idx <- results_df$Rho > 0 & !is.na(results_df$permpval)
rank_pos <- -log10(results_df$permpval[pos_idx])
names(rank_pos) <- results_df$CNEE[pos_idx]
# Break ties by parametric P
tie_break <- -log10(results_df$P[pos_idx])
names(tie_break) <- results_df$CNEE[pos_idx]
rank_pos <- rank_pos + tie_break * 1e-6  # tiny tie-breaker
rank_pos <- sort(rank_pos, decreasing = TRUE)
message(sprintf("Positive Rho: %d CNEEs", length(rank_pos)))

# Negative Rho ranking (constrained in innovative species)
neg_idx <- results_df$Rho < 0 & !is.na(results_df$permpval)
rank_neg <- -log10(results_df$permpval[neg_idx])
names(rank_neg) <- results_df$CNEE[neg_idx]
tie_break_neg <- -log10(results_df$P[neg_idx])
names(tie_break_neg) <- results_df$CNEE[neg_idx]
rank_neg <- rank_neg + tie_break_neg * 1e-6
rank_neg <- sort(rank_neg, decreasing = TRUE)
message(sprintf("Negative Rho: %d CNEEs", length(rank_neg)))

# ============================================================
# E3: Run fgsea for each direction
# ============================================================

message("=== E3: Running fgsea ===")

# Positive Rho — main categories
fgsea_pos <- fgsea(element_sets, rank_pos, minSize = 10, maxSize = 50000)
fgsea_pos$direction <- "positive"
fgsea_pos <- as.data.frame(fgsea_pos)
fgsea_pos <- fgsea_pos[order(fgsea_pos$pval), ]

# Negative Rho — main categories
fgsea_neg <- fgsea(element_sets, rank_neg, minSize = 10, maxSize = 50000)
fgsea_neg$direction <- "negative"
fgsea_neg <- as.data.frame(fgsea_neg)
fgsea_neg <- fgsea_neg[order(fgsea_neg$pval), ]

message("Positive Rho — element enrichment:")
print(fgsea_pos[, c("pathway", "NES", "pval", "padj")])
message("\nNegative Rho — element enrichment:")
print(fgsea_neg[, c("pathway", "NES", "pval", "padj")])

# Sub-categories
fgsea_pos_sub <- fgsea(element_subsets, rank_pos, minSize = 10, maxSize = 50000)
fgsea_pos_sub$direction <- "positive"
fgsea_pos_sub <- as.data.frame(fgsea_pos_sub)

fgsea_neg_sub <- fgsea(element_subsets, rank_neg, minSize = 10, maxSize = 50000)
fgsea_neg_sub$direction <- "negative"
fgsea_neg_sub <- as.data.frame(fgsea_neg_sub)

# Save
write.csv(fgsea_pos, file.path(OUTPUT_DIR, "element_enrichment_pos_rho.csv"), row.names = FALSE)
write.csv(fgsea_neg, file.path(OUTPUT_DIR, "element_enrichment_neg_rho.csv"), row.names = FALSE)
write.csv(fgsea_pos_sub, file.path(OUTPUT_DIR, "element_enrichment_pos_rho_subcat.csv"), row.names = FALSE)
write.csv(fgsea_neg_sub, file.path(OUTPUT_DIR, "element_enrichment_neg_rho_subcat.csv"), row.names = FALSE)

# ============================================================
# E4: Fisher's exact test at FDR<0.05 threshold (secondary)
# ============================================================

message("=== E4: Fisher's exact test (threshold-based) ===")

# Significant CNEEs by global FDR
sig_cnees <- results_df$CNEE[results_df$sig_global_fdr]
sig_pos <- results_df$CNEE[results_df$sig_global_fdr & results_df$Rho > 0]
sig_neg <- results_df$CNEE[results_df$sig_global_fdr & results_df$Rho < 0]
all_cnees <- results_df$CNEE

message(sprintf("Significant: %d total (pos: %d, neg: %d)",
                length(sig_cnees), length(sig_pos), length(sig_neg)))

fisher_results <- do.call(rbind, lapply(categories, function(cat) {
  cat_cnees <- anno_df$name[anno_df$category == cat]
  in_cat <- all_cnees %in% cat_cnees

  for (dir in c("positive", "negative")) {
    sig <- if (dir == "positive") all_cnees %in% sig_pos else all_cnees %in% sig_neg

    # 2x2 table: sig x in_cat
    tt <- table(factor(sig, levels=c(TRUE, FALSE)), factor(in_cat, levels=c(TRUE, FALSE)))
    ft <- fisher.test(tt)

    data.frame(
      category = cat,
      direction = dir,
      n_sig_in_cat = tt[1, 1],
      n_sig_not_cat = tt[1, 2],
      n_notsig_in_cat = tt[2, 1],
      n_notsig_not_cat = tt[2, 2],
      odds_ratio = ft$estimate,
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}))

# BH correction
fisher_results$padj <- p.adjust(fisher_results$pvalue, method = "BH")
fisher_results <- fisher_results[order(fisher_results$pvalue), ]

write.csv(fisher_results, file.path(OUTPUT_DIR, "element_enrichment_fisher.csv"), row.names = FALSE)
message(sprintf("Saved: %s/element_enrichment_fisher.csv", OUTPUT_DIR))

message("\n=== Phase E complete ===")
