#!/usr/bin/env Rscript
# Differential accessibility analysis with DESeq2
# One-vs-all per cell type using sum-to-zero contrasts (contr.sum)
# Size factors: poscounts (for sparse ATAC data)

suppressMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
})

WORKSPACE <- Sys.getenv("WORKSPACE", "/workspace")
n_workers <- as.integer(Sys.getenv("N_WORKERS", "8"))
register(MulticoreParam(n_workers))

# --- Load data ---
cat("Loading pseudobulk counts...\n")
pb_counts <- fread(file.path(WORKSPACE, "pseudobulk_counts.tsv"), sep = "\t")
peak_ids <- pb_counts[[1]]
counts <- as.matrix(pb_counts[, -1])
rm(pb_counts); gc()

pb_meta <- fread(file.path(WORKSPACE, "pseudobulk_metadata.tsv"), sep = "\t")
counts <- counts[, pb_meta$pseudobulk_id]
cat("  Matrix:", nrow(counts), "peaks x", ncol(counts), "pseudobulks\n")

# --- Pre-filter: keep peaks with >= 10 counts in >= 10 samples ---
keep <- rowSums(counts >= 10) >= 10
cat("  Keeping", sum(keep), "of", length(keep), "peaks\n")
counts <- counts[keep, ]
peak_ids <- peak_ids[keep]
rm(keep); gc()

# --- Prepare colData ---
colData <- data.frame(
  row.names = pb_meta$pseudobulk_id,
  cell_type = factor(pb_meta$cell_type),
  sample = factor(pb_meta$sample)
)

# --- Sum-to-zero contrasts on cell_type ---
n_levels <- nlevels(colData$cell_type)
contrasts(colData$cell_type) <- contr.sum(n_levels)
ct_levels <- levels(colData$cell_type)

# --- Create DESeqDataSet ---
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ sample + cell_type
)

# --- Estimate size factors (poscounts for sparse data) ---
dds <- estimateSizeFactors(dds, type = "poscounts")
cat("  Size factors range:", round(range(sizeFactors(dds)), 3), "\n")

# --- Run DESeq ---
cat("Running DESeq()...\n")
dds <- DESeq(dds, parallel = TRUE)

# --- Extract one-vs-all results for each cell type ---
coef_names <- resultsNames(dds)
ct_coef_names <- coef_names[grep("^cell_type", coef_names)]
cat("  Cell type coefficients:", length(ct_coef_names), "\n")

all_results <- list()

for (i in 1:(n_levels - 1)) {
  ct_name <- ct_levels[i]
  res <- results(dds, name = ct_coef_names[i], parallel = TRUE)
  res_df <- as.data.frame(res)
  res_df$peak_id <- peak_ids
  res_df$cell_type <- ct_name
  all_results[[ct_name]] <- res_df
}

# Last level: deviation = -sum of all other cell_type coefficients
last_ct <- ct_levels[n_levels]
contrast_vec <- rep(0, length(coef_names))
names(contrast_vec) <- coef_names
contrast_vec[ct_coef_names] <- -1
res <- results(dds, contrast = contrast_vec, parallel = TRUE)
res_df <- as.data.frame(res)
res_df$peak_id <- peak_ids
res_df$cell_type <- last_ct
all_results[[last_ct]] <- res_df

# --- Combine results ---
all_results_df <- do.call(rbind, all_results)
rownames(all_results_df) <- NULL
colnames(all_results_df)[colnames(all_results_df) == "log2FoldChange"] <- "log2FC"

all_results_df$significant <- !is.na(all_results_df$padj) &
                               all_results_df$padj < 0.05 &
                               all_results_df$log2FC > 1

# --- Save ---
fwrite(all_results_df, file.path(WORKSPACE, "da_results_full.tsv"), sep = "\t")

summary_by_ct <- as.data.table(all_results_df)[significant == TRUE, .(n_da = .N), by = cell_type]
setorder(summary_by_ct, -n_da)
fwrite(summary_by_ct, file.path(WORKSPACE, "da_summary.tsv"), sep = "\t")

cat("\nTotal tests:", nrow(all_results_df), "\n")
cat("Significant (padj<0.05 & log2FC>1):", sum(all_results_df$significant, na.rm = TRUE), "\n")
