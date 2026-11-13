#!/usr/bin/env Rscript
# 04_exon_cleanup_robust.R
# Exon cleanup + robust CNEE selection (nominal significance + sign stability).
#
# NOTE (reconstruction): the original script (v3_steps_1_2.R) ran on a deleted
# machine and was lost. This version was reconstructed from the documented
# logic and VERIFIED against the surviving outputs:
#   - rerconverge_results_cont_v2_noexon.csv  (272,730 rows = 272,811 - 81 exonic)
#   - robust_cnee_sets.csv                    (7,350 rows, sign_stability = 1)
# Sign stability requires leave_out_results.rds from script 05; if that file is
# absent, pass --use-existing-stability to read the column from the surviving
# robust_cnee_sets.csv (verification mode).

suppressPackageStartupMessages({
  library(data.table)
})

BASE        <- "/mnt/results/rerconverge_cont_v2"
RESULTS_DIR <- file.path(BASE, "phase_b_permutations")
ANNO_DIR    <- file.path(BASE, "phase_d_gsea")
WORKSPACE   <- "/workspace"

USE_EXISTING_STABILITY <- "--use-existing-stability" %in% commandArgs(trailingOnly = TRUE)

# ============================================================
# Step 1: Exon cleanup
# ============================================================
message("=== Step 1: Exon cleanup ===")

res <- readRDS(file.path(WORKSPACE, "res_final_cont_v2.rds"))
res$CNEE <- rownames(res)
message(sprintf("  Full results: %d CNEEs", nrow(res)))

anno_full <- fread(file.path(ANNO_DIR, "cnee_annotation_cont_full.csv"))
exonic <- anno_full$name[grepl("^Exon", anno_full$annotation)]
message(sprintf("  Exon-annotated CNEEs (ChIPseeker): %d", length(exonic)))

res_noexon <- res[!res$CNEE %in% exonic, ]
message(sprintf("  After exon cleanup: %d CNEEs (removed %d)",
                nrow(res_noexon), nrow(res) - nrow(res_noexon)))

# Annotation without exons (input for gene-level stats)
anno_noexon <- anno_full[!anno_full$name %in% exonic, ]
fwrite(anno_noexon, file.path(ANNO_DIR, "cnee_annotation_cont_noexon.csv"))
message(sprintf("  Saved: %s", file.path(ANNO_DIR, "cnee_annotation_cont_noexon.csv")))

write.csv(res_noexon[, c("CNEE", "Rho", "N", "P", "p.adj", "permpval",
                         "permpvaladj_bh", "observed_stat", "global_fdr",
                         "sig_global_fdr")],
          file.path(RESULTS_DIR, "rerconverge_results_cont_v2_noexon.csv"),
          row.names = FALSE)

# ============================================================
# Step 2: Robust CNEE selection
# ============================================================
message("\n=== Step 2: Robust CNEE selection ===")

# 2a: nominal significance
nominal <- res_noexon[!is.na(res_noexon$permpval) & res_noexon$permpval < 0.05, ]
message(sprintf("  Nominal significance (permpval < 0.05): %d CNEEs", nrow(nominal)))

# 2b: sign stability across 92 leave-one-out iterations (script 05 output)
lo_path <- file.path(WORKSPACE, "leave_out_results.rds")
if (!USE_EXISTING_STABILITY && file.exists(lo_path)) {
  message("  Computing sign stability from leave_out_results.rds ...")
  results_lo <- readRDS(lo_path)
  sign_stability <- sapply(nominal$CNEE, function(cnee) {
    ref_sign <- sign(nominal$Rho[nominal$CNEE == cnee])
    lo_signs <- vapply(results_lo, function(x) sign(x$rho[cnee]), numeric(1))
    mean(lo_signs == ref_sign, na.rm = TRUE)
  })
} else {
  message("  !! leave_out_results.rds not available — using sign stability ",
          "from surviving robust_cnee_sets.csv (verification mode)")
  prev <- fread(file.path(RESULTS_DIR, "robust_cnee_sets.csv"),
                select = c("CNEE", "sign_stability"))
  sign_stability <- prev$sign_stability[match(nominal$CNEE, prev$CNEE)]
}

nominal$sign_stability <- as.numeric(sign_stability)
robust <- nominal[nominal$sign_stability >= 0.8, ]
robust$direction <- ifelse(robust$Rho > 0, "positive", "negative")
robust <- robust[order(robust$permpval), ]

message(sprintf("  Robust (>= 80%% sign stability): %d CNEEs", nrow(robust)))
message(sprintf("    positive Rho: %d", sum(robust$direction == "positive")))
message(sprintf("    negative Rho: %d", sum(robust$direction == "negative")))

write.csv(robust[, c("CNEE", "Rho", "permpval", "P", "observed_stat",
                     "sign_stability", "direction")],
          file.path(RESULTS_DIR, "robust_cnee_sets.csv"),
          row.names = FALSE)

# ============================================================
# Step 3: Annotate robust CNEEs (category, geneId, distanceToTSS)
# ============================================================
message("\n=== Step 3: Annotating robust CNEEs ===")

anno_noexon <- fread(file.path(ANNO_DIR, "cnee_annotation_cont_noexon.csv"))
anno_lookup <- anno_noexon[match(robust$CNEE, anno_noexon$name), ]

robust_annot <- data.frame(
  CNEE           = robust$CNEE,
  Rho            = robust$Rho,
  permpval       = robust$permpval,
  P              = robust$P,
  observed_stat  = robust$observed_stat,
  sign_stability = robust$sign_stability,
  direction      = robust$direction,
  annotation     = anno_lookup$annotation,
  geneId         = anno_lookup$geneId,
  distanceToTSS  = anno_lookup$distanceToTSS,
  stringsAsFactors = FALSE
)
# Simplified category (first token of ChIPseeker annotation)
robust_annot$category <- sub(" \\(.*\\)", "", robust_annot$annotation)
robust_annot$category[robust_annot$category == "Promoter"] <- "Promoter"

ATTRIBUTED <- c("Intron", "5' UTR", "3' UTR", "Promoter", "Downstream")
robust_annot$attributed <- robust_annot$category %in% ATTRIBUTED

fwrite(robust_annot, file.path(BASE, "robust_cnees_annotated_v3.csv"))
message(sprintf("  Saved: %s (%d rows)", file.path(BASE, "robust_cnees_annotated_v3.csv"), nrow(robust_annot)))

message("\n=== Step 04 complete ===")
