# phase_b_permutations.R
# Phase B: Full permulations on ALL 272,669 CNEEs + global empirical FDR
# 2,000 simperm permulations, streaming approach with two accumulators

suppressPackageStartupMessages({
  library(RERconverge)
  library(ape)
  library(parallel)
})

R_LIBS_USER <- "/workspace/.Rlib"
NCORES <- 8
BATCH.SIZE <- 200
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_b_permutations"
NPERMS <- 2000
PERM_BATCH <- 80  # perms per mclapply batch

# ============================================================
# B1: Load best transform's RERmat and treesObj
# ============================================================

message("=== B1: Loading data ===")

# Read best transform
best_transform <- readLines(file.path("/mnt/results/rerconverge_cont_v2/phase_a_transform", "best_transform.txt"))
message(sprintf("Best transform: %s", best_transform))

# Load phenotype
pheno_df <- read.csv("/mnt/user-uploads/df_merged_for_rer_no_na.csv", row.names = 1)
tip.vals <- pheno_df$TECHINNO2025_ResEff
names(tip.vals) <- rownames(pheno_df)

# Load treesObj
treesObj <- readRDS("/workspace/treesObj_cont.rds")
useSpecies <- names(tip.vals)
masterTree <- treesObj$masterTree
# Root tree for simpermvec (char2Paths works with unrooted, simpermvec needs rooted)
masterTree_rooted <- root(masterTree, outgroup = "GCF_028389875", resolve.root = TRUE)
message(sprintf("Tree rooted: %s", is.rooted(masterTree_rooted)))

# Load RERmat for best transform
RERmat <- readRDS(sprintf("/workspace/RERmat_%s.rds", best_transform))
message(sprintf("RERmat: %d CNEEs x %d branches", nrow(RERmat), ncol(RERmat)))

# ============================================================
# B2: Compute observed statistics
# ============================================================

message("=== B2: Computing observed statistics ===")

charP <- char2Paths(tip.vals, treesObj)
res <- correlateWithContinuousPhenotype(RERmat, charP,
                                         min.sp = 10,
                                         winsorizeRER = 3,
                                         winsorizetrait = 3)
res$p.adj <- p.adjust(res$P, method = "BH")

# Observed signed statistic
obs_stat <- sign(res$Rho) * (-log10(res$P))
obs_stat[is.na(obs_stat)] <- 0
names(obs_stat) <- rownames(res)

# Sort observed |stat| descending — these are the thresholds for FDR curve
obs_abs <- abs(obs_stat)
obs_abs_sorted <- sort(obs_abs, decreasing = TRUE)
N <- length(obs_abs_sorted)

message(sprintf("Observed stats: N=%d, range [%.4f, %.4f]", N, min(obs_abs), max(obs_abs)))
message(sprintf("Parametric: P<0.05: %d, p.adj<0.05: %d",
                sum(res$P < 0.05, na.rm = TRUE),
                sum(res$p.adj < 0.05, na.rm = TRUE)))

# ============================================================
# B3: Streaming permulations (2,000 simperm, all CNEEs)
# ============================================================

message(sprintf("=== B3: Running %d permulations on %d CNEEs ===", NPERMS, N))

# Two accumulators
# Accumulator 1: per-CNEE counter (for permpval)
counters <- rep(0L, nrow(RERmat))
names(counters) <- rownames(RERmat)

# Accumulator 2: total N_perm at each observed threshold (for global FDR)
# Indexed by the sorted observed |stat| thresholds
total_N_perm <- rep(0, N)

# Process permulations in batches
n_batches <- ceiling(NPERMS / PERM_BATCH)
message(sprintf("Processing in %d batches of %d perms (%d cores)", n_batches, PERM_BATCH, NCORES))

perm_done <- 0
start_time <- Sys.time()

for (batch_idx in 1:n_batches) {
  batch_start <- (batch_idx - 1) * PERM_BATCH + 1
  batch_end <- min(batch_idx * PERM_BATCH, NPERMS)
  batch_perms <- batch_start:batch_end
  batch_n <- length(batch_perms)

  # Run batch in parallel
  batch_results <- mclapply(batch_perms, function(perm_idx) {
    vec <- simpermvec(tip.vals, masterTree_rooted)
    paths <- char2Paths(na.omit(vec), treesObj)
    out <- getAllCor(RERmat, paths, method = "p", min.pos = 0,
                     winsorizeRER = 3, winsorizetrait = 3)
    perm_stat <- sign(out$Rho) * (-log10(out$P))
    perm_stat[is.na(perm_stat)] <- 0
    return(perm_stat)
  }, mc.cores = NCORES, mc.preschedule = TRUE)

  # Update accumulators
  for (perm_stat in batch_results) {
    perm_abs <- abs(perm_stat)

    # Accumulator 1: per-CNEE counter
    counters <- counters + as.integer(perm_abs >= obs_abs)

    # Accumulator 2: N_perm at each observed threshold
    perm_abs_sorted <- sort(perm_abs)
    # N - findInterval(t, sorted) = count of perm values >= t
    N_perm_at_obs <- length(perm_abs) - findInterval(obs_abs_sorted, perm_abs_sorted)
    total_N_perm <- total_N_perm + N_perm_at_obs
  }

  perm_done <- perm_done + batch_n
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate <- perm_done / elapsed
  remaining <- (NPERMS - perm_done) / rate
  message(sprintf("Batch %d/%d: %d/%d perms done (%.1f perms/s, ETA %.0f sec / %.1f min)",
                  batch_idx, n_batches, perm_done, NPERMS, rate, remaining, remaining / 60))

  rm(batch_results)
  gc()
}

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
message(sprintf("Permulations complete: %d perms in %.1f sec (%.1f min)", NPERMS, total_time, total_time / 60))

# ============================================================
# B4: Compute significance metrics
# ============================================================

message("=== B4: Computing significance metrics ===")

# Per-CNEE permpval (secondary)
permpval <- (counters + 1) / (NPERMS + 1)
names(permpval) <- names(counters)

# Per-CNEE BH-adjusted permpval (with n = all valid CNEEs)
permpvaladj_bh <- p.adjust(permpval, method = "BH")

# Global empirical FDR (primary)
mean_N_perm <- total_N_perm / NPERMS
N_obs <- 1:N  # since obs_abs_sorted is sorted descending, rank k = k CNEEs with |stat| >= threshold
global_fdr <- mean_N_perm / N_obs

# Enforce monotonicity: cumulative minimum from the right (most stringent threshold)
# FDR should be non-decreasing as we include more CNEEs (less stringent threshold)
# Actually, FDR can be non-monotonic. We use the BH-style step-up procedure:
# Find the largest k where FDR(k) <= 0.05, then all CNEEs with rank <= k are significant
global_fdr_trimmed <- pmin(global_fdr, 1)

# Step-up: find the largest k where FDR(k) <= alpha
alpha <- 0.05
sig_k <- max(which(global_fdr_trimmed <= alpha), 0)
message(sprintf("Global FDR threshold: rank k=%d (FDR=%.4f)", sig_k, global_fdr_trimmed[sig_k]))

# Map FDR back to each CNEE
# obs_abs_sorted is sorted descending; the FDR at rank k corresponds to threshold obs_abs_sorted[k]
# Each CNEE's FDR = global_fdr at its rank
cnee_rank <- match(names(obs_stat), names(obs_abs_sorted))
cnee_global_fdr <- global_fdr_trimmed[cnee_rank]
names(cnee_global_fdr) <- names(obs_stat)

# Significant CNEEs
sig_global_fdr <- cnee_global_fdr <= alpha
sig_global_fdr[is.na(sig_global_fdr)] <- FALSE

message(sprintf("Significant CNEEs (global FDR < %.2f): %d", alpha, sum(sig_global_fdr)))
message(sprintf("Per-CNEE permpval: min=%.6f, <0.05: %d, BH<0.05: %d",
                min(permpval), sum(permpval < 0.05), sum(permpvaladj_bh < 0.05)))

# ============================================================
# B5: Save results
# ============================================================

message("=== B5: Saving results ===")

# Build results data frame
results_df <- data.frame(
  CNEE = rownames(res),
  Rho = res$Rho,
  N = res$N,
  P = res$P,
  p.adj = res$p.adj,
  permpval = permpval[rownames(res)],
  permpvaladj_bh = permpvaladj_bh[rownames(res)],
  observed_stat = obs_stat[rownames(res)],
  global_fdr = cnee_global_fdr[rownames(res)],
  sig_global_fdr = sig_global_fdr[rownames(res)],
  stringsAsFactors = FALSE
)

write.csv(results_df, file.path(OUTPUT_DIR, "rerconverge_results_cont_v2.csv"), row.names = FALSE)
message(sprintf("Saved: %s/rerconverge_results_cont_v2.csv (%d rows)", OUTPUT_DIR, nrow(results_df)))

# Save FDR curve
fdr_curve <- data.frame(
  rank = 1:N,
  threshold = obs_abs_sorted,
  N_obs = N_obs,
  N_perm_mean = mean_N_perm,
  FDR = global_fdr_trimmed
)
write.csv(fdr_curve, file.path(OUTPUT_DIR, "global_fdr_curve.csv"), row.names = FALSE)
message(sprintf("Saved: %s/global_fdr_curve.csv (%d rows)", OUTPUT_DIR, nrow(fdr_curve)))

# Save comparison with previous analysis
prev_results <- read.csv("/mnt/results/rerconverge_results_cont.csv")
comparison <- data.frame(
  metric = c("Total valid CNEEs", "Parametric P<0.05", "Parametric p.adj<0.05",
             "Permpval<0.05", "Permpval BH<0.05 (n=all)",
             "Global FDR<0.05", "Min permpval"),
  previous = c(sum(!is.na(prev_results$P)),
               sum(prev_results$P < 0.05, na.rm = TRUE),
               sum(prev_results$p.adj < 0.05, na.rm = TRUE),
               sum(prev_results$permpval < 0.05, na.rm = TRUE),
               sum(prev_results$permpvaladj < 0.05, na.rm = TRUE),
               NA,
               min(prev_results$permpval, na.rm = TRUE)),
  current = c(sum(!is.na(results_df$P)),
              sum(results_df$P < 0.05, na.rm = TRUE),
              sum(results_df$p.adj < 0.05, na.rm = TRUE),
              sum(results_df$permpval < 0.05, na.rm = TRUE),
              sum(results_df$permpvaladj_bh < 0.05, na.rm = TRUE),
              sum(results_df$sig_global_fdr),
              min(results_df$permpval))
)
write.csv(comparison, file.path(OUTPUT_DIR, "permutation_comparison.csv"), row.names = FALSE)
message(sprintf("Saved: %s/permutation_comparison.csv", OUTPUT_DIR))

# Save key results as RDS for downstream phases
saveRDS(results_df, "/workspace/res_final_cont_v2.rds")
saveRDS(fdr_curve, "/workspace/fdr_curve_cont_v2.rds")

message("\n=== Phase B complete ===")
