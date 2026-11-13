# phase_f_dar.R
# Phase F: Open chromatin (DAR) enrichment
# fgsea on 30 cell type DAR overlaps, direction-specific
# Fisher's exact as secondary

suppressPackageStartupMessages({
  library(fgsea)
  library(GenomicRanges)
  library(ggplot2)
})

R_LIBS_USER <- "/workspace/.Rlib"
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_f_dar"

# ============================================================
# F1: Load DAR regions and compute CNEE overlaps
# ============================================================

message("=== F1: Computing CNEE-DAR overlaps ===")

# Load CNEE BED
cnee_bed <- read.table("/mnt/user-uploads/lifted_annotations_genomic.bed",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(cnee_bed) <- c("chr", "start", "end", "name")
cnee_gr <- makeGRangesFromDataFrame(cnee_bed, keep.extra.columns = TRUE,
                                     starts.in.df.are.0based = TRUE)

# Get all DAR BED files
dar_files <- list.files("/mnt/user-uploads/bed_files", pattern = "_accessible_regions.bed$", full.names = TRUE)
message(sprintf("Found %d DAR files", length(dar_files)))

# Compute overlaps for each cell type
dar_sets <- list()
overlap_matrix <- matrix(0L, nrow = nrow(cnee_bed), ncol = length(dar_files),
                          dimnames = list(cnee_bed$name, basename(dar_files)))

for (i in seq_along(dar_files)) {
  f <- dar_files[i]
  cell_name <- sub("_accessible_regions.bed$", "", basename(f))

  # Read BED (skip 3 comment lines)
  dar <- read.table(f, skip = 3, header = FALSE, stringsAsFactors = FALSE)
  colnames(dar) <- c("chr", "start", "end", "peak_id", "log2FC", "neg_log10_padj")
  dar_gr <- makeGRangesFromDataFrame(dar, starts.in.df.are.0based = TRUE)

  # Find overlaps
  hits <- findOverlaps(cnee_gr, dar_gr)
  overlapping_cnees <- unique(cnee_bed$name[queryHits(hits)])

  dar_sets[[cell_name]] <- overlapping_cnees
  overlap_matrix[, i] <- as.integer(cnee_bed$name %in% overlapping_cnees)

  message(sprintf("  %s: %d DAR regions, %d CNEEs overlapping (%.1f%%)",
                  cell_name, nrow(dar), length(overlapping_cnees),
                  100 * length(overlapping_cnees) / nrow(cnee_bed)))
}

# Save overlap matrix
colnames(overlap_matrix) <- sub("_accessible_regions.bed$", "", basename(dar_files))
write.csv(overlap_matrix, file.path(OUTPUT_DIR, "dar_overlap_matrix.csv"))
message(sprintf("Saved: %s/dar_overlap_matrix.csv", OUTPUT_DIR))

# ============================================================
# F2: Run fgsea for each direction (primary, p-value-aware)
# ============================================================

message("=== F2: Running fgsea for DAR enrichment ===")

# Load Phase B results
results_df <- readRDS("/workspace/res_final_cont_v2.rds")

# Build direction-specific rankings (same as Phase E)
pos_idx <- results_df$Rho > 0 & !is.na(results_df$permpval)
rank_pos <- -log10(results_df$permpval[pos_idx])
names(rank_pos) <- results_df$CNEE[pos_idx]
tie_break <- -log10(results_df$P[pos_idx])
names(tie_break) <- results_df$CNEE[pos_idx]
rank_pos <- sort(rank_pos + tie_break * 1e-6, decreasing = TRUE)

neg_idx <- results_df$Rho < 0 & !is.na(results_df$permpval)
rank_neg <- -log10(results_df$permpval[neg_idx])
names(rank_neg) <- results_df$CNEE[neg_idx]
tie_break_neg <- -log10(results_df$P[neg_idx])
names(tie_break_neg) <- results_df$CNEE[neg_idx]
rank_neg <- sort(rank_neg + tie_break_neg * 1e-6, decreasing = TRUE)

# Run fgsea
fgsea_dar_pos <- fgsea(dar_sets, rank_pos, minSize = 10, maxSize = 50000)
fgsea_dar_pos$direction <- "positive"
fgsea_dar_pos <- as.data.frame(fgsea_dar_pos)
fgsea_dar_pos <- fgsea_dar_pos[order(fgsea_dar_pos$pval), ]

fgsea_dar_neg <- fgsea(dar_sets, rank_neg, minSize = 10, maxSize = 50000)
fgsea_dar_neg$direction <- "negative"
fgsea_dar_neg <- as.data.frame(fgsea_dar_neg)
fgsea_dar_neg <- fgsea_dar_neg[order(fgsea_dar_neg$pval), ]

message("Positive Rho — DAR enrichment (top 10):")
print(head(fgsea_dar_pos[, c("pathway", "NES", "pval", "padj")], 10))
message("\nNegative Rho — DAR enrichment (top 10):")
print(head(fgsea_dar_neg[, c("pathway", "NES", "pval", "padj")], 10))

write.csv(fgsea_dar_pos, file.path(OUTPUT_DIR, "dar_enrichment_pos_rho.csv"), row.names = FALSE)
write.csv(fgsea_dar_neg, file.path(OUTPUT_DIR, "dar_enrichment_neg_rho.csv"), row.names = FALSE)

# ============================================================
# F3: Fisher's exact test (secondary, threshold-based)
# ============================================================

message("=== F3: Fisher's exact test ===")

sig_cnees <- results_df$CNEE[results_df$sig_global_fdr]
sig_pos <- results_df$CNEE[results_df$sig_global_fdr & results_df$Rho > 0]
sig_neg <- results_df$CNEE[results_df$sig_global_fdr & results_df$Rho < 0]
all_cnees <- results_df$CNEE

message(sprintf("Significant: %d total (pos: %d, neg: %d)",
                length(sig_cnees), length(sig_pos), length(sig_neg)))

fisher_dar <- do.call(rbind, lapply(names(dar_sets), function(cell) {
  dar_cnees <- dar_sets[[cell]]
  in_dar <- all_cnees %in% dar_cnees

  for (dir in c("positive", "negative")) {
    sig <- if (dir == "positive") all_cnees %in% sig_pos else all_cnees %in% sig_neg

    tt <- table(factor(sig, levels=c(TRUE, FALSE)), factor(in_dar, levels=c(TRUE, FALSE)))
    ft <- fisher.test(tt)

    data.frame(
      cell_type = cell,
      direction = dir,
      n_sig_in_dar = tt[1, 1],
      n_sig_not_dar = tt[1, 2],
      n_notsig_in_dar = tt[2, 1],
      n_notsig_not_dar = tt[2, 2],
      odds_ratio = ft$estimate,
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}))

fisher_dar$padj <- p.adjust(fisher_dar$pvalue, method = "BH")
fisher_dar <- fisher_dar[order(fisher_dar$pvalue), ]

write.csv(fisher_dar, file.path(OUTPUT_DIR, "dar_fisher_pos_neg.csv"), row.names = FALSE)
message(sprintf("Saved: %s/dar_fisher_pos_neg.csv", OUTPUT_DIR))

# Print significant results
sig_fisher <- fisher_dar[fisher_dar$padj < 0.05, ]
if (nrow(sig_fisher) > 0) {
  message(sprintf("\nSignificant DAR enrichment (Fisher FDR<0.05): %d", nrow(sig_fisher)))
  print(sig_fisher[, c("cell_type", "direction", "odds_ratio", "pvalue", "padj")])
} else {
  message("\nNo significant DAR enrichment (Fisher FDR<0.05)")
}

message("\n=== Phase F complete ===")
