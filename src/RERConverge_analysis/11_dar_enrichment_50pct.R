#!/usr/bin/env Rscript
# 02_dar_enrichment_50pct.R
# DAR enrichment with >50% CNEE length overlap matrix.
# Wilcoxon rank-sum test on -log10(permpval) + rank-sum permutation z-score.
# Directions: negative and positive Rho only (omnibus 'all' removed in v3.1).

suppressPackageStartupMessages({
  library(GenomicRanges)
})

set.seed(42)
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_f_dar"
N_PERM <- 10000

# ============================================================
# Load data
# ============================================================
message("=== Loading data ===")

dar_mat <- read.csv(file.path(OUTPUT_DIR, "dar_overlap_matrix_50pct.csv"),
                    row.names = 1, check.names = FALSE)
message(sprintf("DAR matrix: %d CNEEs x %d cell types", nrow(dar_mat), ncol(dar_mat)))

res <- readRDS("/workspace/res_final_cont_v2.rds")
message(sprintf("Results: %d CNEEs", nrow(res)))

common <- intersect(rownames(dar_mat), res$CNEE)
dar_mat <- dar_mat[common, , drop = FALSE]
res <- res[match(common, res$CNEE), ]
message(sprintf("Merged: %d CNEEs", nrow(res)))

res$neg_log10_perm <- -log10(res$permpval)

cell_types <- colnames(dar_mat)
annotated_cts <- cell_types[!grepl("^Unannotated", cell_types)]
message(sprintf("Cell types: %d total, %d annotated", length(cell_types), length(annotated_cts)))

# ============================================================
# Wilcoxon + rank-sum z-test (NOTE: "perm_*" columns are a tie-corrected NORMAL approximation of the
# same rank-sum statistic, not a permutation test; N_PERM is unused. They are not independent of the Wilcoxon p-values.)
# ============================================================
message("=== Running Wilcoxon density test ===")

directions <- c("negative", "positive")
results_list <- list()

for (ct in cell_types) {
  dar_col <- dar_mat[, ct]

  for (direction in directions) {
    if (direction == "all") {
      idx <- seq_len(nrow(res))
    } else if (direction == "positive") {
      idx <- which(res$Rho > 0)
    } else {
      idx <- which(res$Rho < 0)
    }

    scores <- res$neg_log10_perm[idx]
    in_dar <- dar_col[idx] == 1
    n_in <- sum(in_dar)
    n_out <- sum(!in_dar)
    n_total <- n_in + n_out

    if (n_in < 10 || n_out < 10) {
      results_list[[length(results_list) + 1]] <- data.frame(
        cell_type = ct, direction = direction,
        n_in_dar = n_in, n_out_dar = n_out,
        median_in = NA_real_, median_out = NA_real_,
        W = NA_real_, pvalue_wilcoxon = NA_real_,
        effect_size = NA_real_, padj_wilcoxon = NA_real_,
        observed_ranksum = NA_real_, expected_ranksum = NA_real_,
        perm_pvalue = NA_real_, padj_perm = NA_real_,
        perm_zscore = NA_real_
      )
      next
    }

    scores_in <- scores[in_dar]
    scores_out <- scores[!in_dar]
    median_in <- median(scores_in)
    median_out <- median(scores_out)

    # Wilcoxon rank-sum test (one-sided "greater")
    wt <- wilcox.test(scores_in, scores_out, alternative = "greater")
    pval <- wt$p.value
    W_stat <- as.numeric(wt$statistic)

    # Effect size: rank-biserial correlation
    effect_size <- 2 * W_stat / (n_in * n_out) - 1

    # Rank-sum permutation test
    ranks <- rank(scores)
    observed_ranksum <- sum(ranks[in_dar])
    expected_ranksum <- n_in * (n_total + 1) / 2

    # Tie-corrected variance of rank-sum under null
    tie_tab <- table(scores)
    tie_sum <- sum(as.numeric(tie_tab)^3 - as.numeric(tie_tab))
    var_ranksum <- (n_in * n_out / 12) * (n_total + 1 - tie_sum / (n_total * (n_total - 1)))
    sd_ranksum <- sqrt(var_ranksum)

    perm_z <- (observed_ranksum - expected_ranksum) / sd_ranksum
    perm_pval <- 1 - pnorm(perm_z)

    results_list[[length(results_list) + 1]] <- data.frame(
      cell_type = ct, direction = direction,
      n_in_dar = n_in, n_out_dar = n_out,
      median_in = round(median_in, 6), median_out = round(median_out, 6),
      W = round(W_stat, 2),
      pvalue_wilcoxon = pval,
      effect_size = round(effect_size, 6),
      padj_wilcoxon = NA_real_,
      observed_ranksum = round(observed_ranksum, 2),
      expected_ranksum = round(expected_ranksum, 2),
      perm_pvalue = perm_pval,
      padj_perm = NA_real_,
      perm_zscore = round(perm_z, 4)
    )
  }

  message(sprintf("  %s done", ct))
}

df <- do.call(rbind, results_list)

# BH FDR correction within each direction
for (direction in directions) {
  mask <- df$direction == direction & !is.na(df$pvalue_wilcoxon)
  if (sum(mask) > 0) {
    df$padj_wilcoxon[mask] <- p.adjust(df$pvalue_wilcoxon[mask], method = "BH")
    df$padj_perm[mask] <- p.adjust(df$perm_pvalue[mask], method = "BH")
  }
}

df <- df[order(df$direction, df$pvalue_wilcoxon), ]

# ============================================================
# Print results
# ============================================================
message("\n=== RESULTS: Negative Rho direction ===")
neg <- df[df$direction == "negative", ]
sig_neg <- neg[!is.na(neg$padj_wilcoxon) & neg$padj_wilcoxon < 0.05, ]
message(sprintf("Wilcoxon significant (FDR<0.05): %d", nrow(sig_neg)))
sig_neg_perm <- neg[!is.na(neg$padj_perm) & neg$padj_perm < 0.05, ]
message(sprintf("Permutation significant (FDR<0.05): %d", nrow(sig_neg_perm)))
if (nrow(sig_neg) > 0) {
  print(sig_neg[, c("cell_type", "n_in_dar", "effect_size", "pvalue_wilcoxon", "padj_wilcoxon", "perm_pvalue", "padj_perm", "perm_zscore")])
}

message("\n=== RESULTS: Positive Rho direction ===")
pos <- df[df$direction == "positive", ]
sig_pos <- pos[!is.na(pos$padj_wilcoxon) & pos$padj_wilcoxon < 0.05, ]
message(sprintf("Wilcoxon significant (FDR<0.05): %d", nrow(sig_pos)))

# ============================================================
# Save
# ============================================================
out_path <- file.path(OUTPUT_DIR, "dar_enrichment_summary_v3.csv")
write.csv(df, out_path, row.names = FALSE)
message(sprintf("\nSaved: %s", out_path))

for (direction in directions) {
  sub <- df[df$direction == direction, ]
  sub_path <- file.path(OUTPUT_DIR, sprintf("dar_wilcoxon_%s_v3.csv", direction))
  write.csv(sub, sub_path, row.names = FALSE)
}

message("\n=== Step 12 complete ===")
