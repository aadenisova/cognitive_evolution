#!/usr/bin/env Rscript
# 11_dar_overlap_50pct.R
# Recompute DAR overlap matrix with >50% CNEE length overlap criterion.
# Currently: any 1bp overlap (findOverlaps). New: overlap_width / cnee_width > 0.5.

suppressPackageStartupMessages({
  library(GenomicRanges)
})

OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_f_dar"

# ============================================================
# Load CNEE BED
# ============================================================
message("=== Loading CNEE BED ===")
cnee_bed <- read.table("/mnt/user-uploads/lifted_annotations_genomic.bed",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(cnee_bed) <- c("chr", "start", "end", "name")
cnee_gr <- makeGRangesFromDataFrame(cnee_bed, keep.extra.columns = TRUE,
                                     starts.in.df.are.0based = TRUE)
message(sprintf("CNEEs: %d, mean width: %.1f bp", length(cnee_gr), mean(width(cnee_gr))))

# ============================================================
# Compute overlaps with >50% CNEE length criterion
# ============================================================
message("=== Computing overlaps (>50% CNEE length) ===")

dar_files <- list.files("/mnt/user-uploads/bed_files", pattern = "_accessible_regions.bed$", full.names = TRUE)
message(sprintf("DAR files: %d", length(dar_files)))

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

  # Compute overlap width for each hit
  if (length(hits) > 0) {
    overlaps <- pintersect(cnee_gr[queryHits(hits)], dar_gr[subjectHits(hits)])
    overlap_widths <- width(overlaps)
    cnee_widths <- width(cnee_gr[queryHits(hits)])

    # Keep only CNEEs where overlap > 50% of CNEE length
    frac <- overlap_widths / cnee_widths
    passing <- queryHits(hits)[frac > 0.5]
    passing_cnees <- unique(passing)

    overlap_matrix[passing_cnees, i] <- 1L
    message(sprintf("  %s: %d DARs, %d CNEEs any-overlap, %d CNEEs >50%% overlap (%.2f%%)",
                    cell_name, nrow(dar), length(unique(queryHits(hits))),
                    length(passing_cnees),
                    100 * length(passing_cnees) / nrow(cnee_bed)))
  } else {
    message(sprintf("  %s: %d DARs, 0 overlaps", cell_name, nrow(dar)))
  }
}

# ============================================================
# Save
# ============================================================
colnames(overlap_matrix) <- sub("_accessible_regions.bed$", "", basename(dar_files))
out_path <- file.path(OUTPUT_DIR, "dar_overlap_matrix_50pct.csv")
write.csv(overlap_matrix, out_path)
message(sprintf("\nSaved: %s (%d x %d)", out_path, nrow(overlap_matrix), ncol(overlap_matrix)))

# Compare with old matrix
old_path <- file.path(OUTPUT_DIR, "dar_overlap_matrix.csv")
if (file.exists(old_path)) {
  old_mat <- read.csv(old_path, row.names = 1, check.names = FALSE)
  common_rows <- intersect(rownames(overlap_matrix), rownames(rownames(old_mat)))
  message("\n=== Comparison with old (any-overlap) matrix ===")
  for (ct in colnames(overlap_matrix)) {
    if (ct %in% colnames(old_mat)) {
      old_n <- sum(old_mat[common_rows, ct] == 1)
      new_n <- sum(overlap_matrix[common_rows, ct] == 1)
      message(sprintf("  %s: old=%d, new=%d (%.1f%% of old)", ct, old_n, new_n, 100*new_n/old_n))
    }
  }
}

message("\n=== Step 11 complete ===")
