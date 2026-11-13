# parallel_core.R
# Parallelized version of RERconverge's coreGetResiduals and getAllResiduals
# Replaces the single-threaded for-loop with batched mclapply

library(RERconverge)
library(parallel)

# Parallel version of coreGetResiduals
# Processes each CNEE tree independently using mclapply
parallel_coreGetResiduals <- function(treesObj, nvMod = NULL, n.pcs = 0, cutoff = NULL,
                                       useSpecies = NULL, min.sp = 10, min.valid = 20,
                                       maxT = NULL, weights = NULL, use.weights = TRUE,
                                       ncores = 16, batch.size = 200) {

    # --- Pre-compute shared data (same as original coreGetResiduals) ---
    if (is.null(cutoff)) {
        cutoff <- quantile(treesObj$paths, 0.05, na.rm = TRUE)
        message(paste("cutoff is set to", cutoff))
    }

    tPaths <- treesObj$paths

    if (is.null(nvMod)) {
        if (n.pcs > 0) {
            stop("PC normalization not supported in parallel version. Set n.pcs=0")
        } else {
            message("using average normalization vector")
            tPathsRaw <- treesObj$transformInv(tPaths)
            nvAve <- treesObj$transformFunc(apply(tPathsRaw, 2, mean, na.rm = TRUE, trim = 0.05))
            nvMod <- model.matrix(~1 + nvAve)
        }
    }

    if (any(is.na(nvMod))) {
        stop("NA values in model: something is wrong")
    }

    if (is.null(useSpecies)) {
        useSpecies <- treesObj$masterTree$tip.label
    } else {
        useSpecies <- intersect(treesObj$masterTree$tip.label, useSpecies)
    }

    if (is.null(maxT)) {
        maxT <- treesObj$numTrees
    }

    if (!is.null(weights)) {
        message("using external weights")
    } else if (is.null(weights) && !is.null(treesObj$weights)) {
        weights <- treesObj$weights
        message("using included weights")
    }

    if (is.null(weights) | !use.weights) {
        message("No weights found. Setting weights to 1.")
        weights <- RERconverge:::copyMat(tPaths)
        weights[] <- 1
    }

    cm <- intersect(treesObj$masterTree$tip.label, useSpecies)
    sp.miss <- setdiff(treesObj$masterTree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
        message(paste0("Species from master tree not present in useSpecies: ",
                       paste(sp.miss, collapse = ",")))
    }

    maxSpecies <- rowSums(treesObj$report[, cm])
    preds <- RERconverge:::copyMat(tPaths)
    preds[] <- NA

    # --- Split CNEE indices into batches ---
    allIndices <- 1:maxT
    batches <- split(allIndices, ceiling(seq_along(allIndices) / batch.size))
    message(sprintf("Processing %d CNEEs in %d batches (%d cores, batch size %d)",
                    maxT, length(batches), ncores, batch.size))

    # --- Process batches in parallel ---
    # Each batch processes its CNEEs sequentially and returns preds rows
    processBatch <- function(batchIdx) {
        batchPreds <- matrix(NA_real_, nrow = length(batchIdx), ncol = ncol(tPaths))
        batchModelIndex <- vector("list", length(batchIdx))
        rownames(batchPreds) <- rownames(tPaths)[batchIdx]

        for (k in seq_along(batchIdx)) {
            i <- batchIdx[k]
            tree1 <- treesObj$trees[[i]]
            both <- intersect(tree1$tip.label, cm)

            if (length(both) < min.sp) {
                next
            }

            tree1 <- unroot(pruneTree(tree1, both))
            ee <- RERconverge:::edgeIndexRelativeMasterTT(tree1, treesObj$masterTree)
            iiPaths <- treesObj$matIndex[ee[, c(2, 1)]]

            allbranch <- tPaths[i, iiPaths, drop = FALSE]
            allbranchw <- weights[i, iiPaths, drop = FALSE]

            iibad <- which(allbranch < cutoff)
            if (length(iibad) > 0) {
                allbranch[iibad] <- NA
            }

            sumCanUse <- sum(!is.na(allbranch))
            if (sumCanUse < min.valid) {
                next
            }

            iiallnotna <- which(!is.na(tPaths[i, ]))
            modelIdx <- iiPaths[which(!is.na(allbranch))]
            batchModelIndex[[k]] <- modelIdx

            if (length(modelIdx) > 0) {
                batchPreds[k, iiallnotna] <- RERconverge:::fastLmResidMatWeightedPredict(
                    matrix(allbranch, nrow = 1),
                    nvMod[iiPaths, , drop = FALSE],
                    matrix(allbranchw, nrow = 1),
                    nvMod[iiallnotna, , drop = FALSE]
                )
            }
        }

        list(indices = batchIdx, preds = batchPreds, modelIndex = batchModelIndex)
    }

    results <- mclapply(batches, processBatch, mc.cores = ncores, mc.preschedule = TRUE)

    # --- Combine results ---
    modelIndexList <- vector("list", nrow(preds))
    for (r in results) {
        preds[r$indices, ] <- r$preds
        for (k in seq_along(r$indices)) {
            if (length(r$modelIndex[[k]]) > 0) {
                modelIndexList[[r$indices[k]]] <- r$modelIndex[[k]]
            }
        }
    }

    # --- Compute residuals (same as original) ---
    allresiduals <- (tPaths - preds)
    rownames(allresiduals) <- names(treesObj$trees)

    if (!is.null(weights)) {
        message("studentizing residuals with weights")
        allresiduals <- allresiduals * sqrt(weights)
    }

    return(list(allresiduals = allresiduals, index = modelIndexList,
                preds = preds, nvMod = nvMod))
}

# Parallel version of getAllResiduals
parallel_getAllResiduals <- function(treesObj, transform = "sqrt", impute = TRUE,
                                      nvMod = NULL, n.pcs = 0, cutoff = NULL,
                                      useSpecies = NULL, min.sp = 10, min.valid = 20,
                                      maxT = NULL, weights = NULL, use.weights = TRUE,
                                      norm = "scale", all = FALSE, use.rows = NULL,
                                      ncores = 16, batch.size = 200) {

    message(sprintf("parallel_getAllResiduals: transform=%s, ncores=%d, batch.size=%d",
                    transform, ncores, batch.size))

    # transformPaths modifies treesObj in place — work on a copy
    tree2 <- transformPaths(treesObj, transform = transform, impute = impute)

    resids <- parallel_coreGetResiduals(tree2, nvMod = nvMod, n.pcs = n.pcs,
                                         cutoff = cutoff, useSpecies = useSpecies,
                                         min.sp = min.sp, min.valid = min.valid,
                                         maxT = maxT, weights = weights,
                                         use.weights = use.weights,
                                         ncores = ncores, batch.size = batch.size)

    r <- RERconverge:::getRMat(resids, all = all, use.rows = use.rows, norm = norm)
    dev.off()
    return(r)
}
