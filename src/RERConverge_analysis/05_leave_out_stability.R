# phase_c_leave_out.R
# Phase C: Leave-one-species-out robustness analysis
# 92 iterations, removing one species at a time

suppressPackageStartupMessages({
  library(RERconverge)
  library(parallel)
})

R_LIBS_USER <- "/workspace/.Rlib"
NCORES <- 8
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_c_leave_out"

# ============================================================
# C1: Load data and compute reference analysis
# ============================================================

message("=== C1: Loading data and computing reference ===")

# Load phenotype
pheno_df <- read.csv("/mnt/user-uploads/df_merged_for_rer_no_na.csv", row.names = 1)
tip.vals <- pheno_df$TECHINNO2025_ResEff
names(tip.vals) <- rownames(pheno_df)

# Load treesObj and RERmat
treesObj <- readRDS("/workspace/treesObj_cont.rds")
best_transform <- readLines(file.path("/mnt/results/rerconverge_cont_v2/phase_a_transform", "best_transform.txt"))
RERmat <- readRDS(sprintf("/workspace/RERmat_%s.rds", best_transform))

# Reference analysis (all 92 species)
charP_ref <- char2Paths(tip.vals, treesObj)
res_ref <- correlateWithContinuousPhenotype(RERmat, charP_ref,
                                             min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
rho_ref <- res_ref$Rho
p_ref <- res_ref$P
names(rho_ref) <- rownames(res_ref)
names(p_ref) <- rownames(res_ref)

message(sprintf("Reference: %d CNEEs, Rho range [%.4f, %.4f]", length(rho_ref), min(rho_ref, na.rm=TRUE), max(rho_ref, na.rm=TRUE)))

# Top CNEE sets by |Rho|
top100_ref <- names(sort(abs(rho_ref), decreasing=TRUE))[1:100]
top500_ref <- names(sort(abs(rho_ref), decreasing=TRUE))[1:500]
top1000_ref <- names(sort(abs(rho_ref), decreasing=TRUE))[1:1000]

# ============================================================
# C2: Leave-out loop (92 iterations, parallel)
# ============================================================

message(sprintf("=== C2: Leave-one-out (%d species, %d cores) ===", length(tip.vals), NCORES))

species_list <- names(tip.vals)

leave_out_one <- function(sp, tip.vals, treesObj, RERmat) {
  tip.vals_lo <- tip.vals
  tip.vals_lo[sp] <- NA

  charP_lo <- char2Paths(na.omit(tip.vals_lo), treesObj)
  res_lo <- correlateWithContinuousPhenotype(RERmat, charP_lo,
                                              min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)

  rho_lo <- res_lo$Rho
  p_lo <- res_lo$P
  names(rho_lo) <- rownames(res_lo)
  names(p_lo) <- rownames(res_lo)

  list(rho = rho_lo, p = p_lo)
}

results_lo <- mclapply(species_list, leave_out_one,
                       tip.vals = tip.vals, treesObj = treesObj, RERmat = RERmat,
                       mc.cores = NCORES, mc.preschedule = TRUE)
names(results_lo) <- species_list

message("Leave-out iterations complete")

# ============================================================
# C3: Robustness analysis
# ============================================================

message("=== C3: Computing robustness metrics ===")

robustness <- data.frame(
  species = species_list,
  phenotype_value = tip.vals[species_list],
  spearman_rho = NA_real_,
  top100_overlap = NA_real_,
  top500_overlap = NA_real_,
  top1000_overlap = NA_real_,
  n_p_lt_0.05 = NA_integer_,
  mean_abs_rho_shift = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  rho_lo <- results_lo[[sp]]$rho
  p_lo <- results_lo[[sp]]$p

  # Spearman correlation of Rho
  common <- intersect(names(rho_ref), names(rho_lo))
  robustness$spearman_rho[i] <- cor(rho_ref[common], rho_lo[common], method = "spearman", use = "complete.obs")

  # Top CNEE overlap
  top100_lo <- names(sort(abs(rho_lo), decreasing=TRUE))[1:100]
  top500_lo <- names(sort(abs(rho_lo), decreasing=TRUE))[1:500]
  top1000_lo <- names(sort(abs(rho_lo), decreasing=TRUE))[1:1000]

  robustness$top100_overlap[i] <- length(intersect(top100_ref, top100_lo)) / 100
  robustness$top500_overlap[i] <- length(intersect(top500_ref, top500_lo)) / 500
  robustness$top1000_overlap[i] <- length(intersect(top1000_ref, top1000_lo)) / 1000

  # Significant CNEE count
  robustness$n_p_lt_0.05[i] <- sum(p_lo < 0.05, na.rm = TRUE)

  # Mean absolute Rho shift
  robustness$mean_abs_rho_shift[i] <- mean(abs(rho_lo[common] - rho_ref[common]), na.rm = TRUE)
}

# Sort by influence (lowest Spearman = most influential)
robustness <- robustness[order(robustness$spearman_rho), ]

print(head(robustness, 10))
print(tail(robustness, 10))

# ============================================================
# C4: Save results
# ============================================================

message("=== C4: Saving results ===")

write.csv(robustness, file.path(OUTPUT_DIR, "leave_out_species_results.csv"), row.names = FALSE)
message(sprintf("Saved: %s/leave_out_species_results.csv", OUTPUT_DIR))

# Save leave-out Rho vectors as RDS for downstream use
saveRDS(results_lo, "/workspace/leave_out_results.rds")

message("\n=== Phase C complete ===")
