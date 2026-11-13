#!/usr/bin/env Rscript
# 04_gsea_v3.R
# GSEA with filtered attribution gene-level stats.
# Split by GO ontology (BP, CC, MF) into separate output files.
# Add gene symbols to leadingEdge column.
# 4 rankings x 2 organisms x 3 ontologies.

suppressPackageStartupMessages({
  library(fgsea)
  library(org.Gg.eg.db)
  library(org.Hs.eg.db)
  library(GO.db)
  library(homologene)
  library(AnnotationDbi)
})

set.seed(42)
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_d_gsea"

# ============================================================
# Load gene-level stats
# ============================================================
message("=== Loading gene-level stats v3 ===")
gene_stats <- read.csv(file.path(OUTPUT_DIR, "..", "gene_level_stats_v3.csv"), stringsAsFactors = FALSE)
message(sprintf("Genes: %d", nrow(gene_stats)))
message(sprintf("With k_pos >= 1: %d, k_neg >= 1: %d",
                sum(gene_stats$k_pos >= 1), sum(gene_stats$k_neg >= 1)))

# ============================================================
# Build GO gene sets split by ontology
# ============================================================
message("=== Building GO gene sets (split by ontology) ===")

build_go_by_ontology <- function(org_db, org_name) {
  go_all <- AnnotationDbi::select(org_db, keys=keys(org_db),
                                  columns=c("GO","ONTOLOGY","SYMBOL"),
                                  keytype="ENTREZID")
  go_all <- go_all[!is.na(go_all$GO), ]

  # Split by ontology
  go_list <- list()
  for (ont in c("BP", "CC", "MF")) {
    sub <- go_all[go_all$ONTOLOGY == ont, ]
    sets <- split(sub$ENTREZID, sub$GO)
    # Filter to size 10-500
    sets <- sets[sapply(sets, function(x) length(unique(x)) >= 10 && length(unique(x)) <= 500)]
    go_list[[ont]] <- sets
    message(sprintf("  %s %s: %d gene sets (size 10-500)", org_name, ont, length(sets)))
  }

  # GO term descriptions
  all_go_ids <- unique(c(names(go_list$BP), names(go_list$CC), names(go_list$MF)))
  desc <- AnnotationDbi::select(GO.db, keys=all_go_ids,
                                columns=c("GOID","TERM","ONTOLOGY"), keytype="GOID")
  desc_map <- setNames(desc$TERM, desc$GOID)
  ont_map <- setNames(desc$ONTOLOGY, desc$GOID)

  # Entrez ID -> Symbol map
  sym_map <- setNames(go_all$SYMBOL, go_all$ENTREZID)

  list(go_list = go_list, desc_map = desc_map, ont_map = ont_map, sym_map = sym_map)
}

chicken_go <- build_go_by_ontology(org.Gg.eg.db, "Chicken")

# Human ortholog mapping
message("Mapping chicken -> human orthologs...")
human_map <- tryCatch({
  homologene::homologene(unique(gene_stats$geneId), inTax=9031, outTax=9606)
}, error=function(e) { message(sprintf("homologene error: %s", e$message)); NULL })

if (!is.null(human_map) && nrow(human_map) > 0) {
  message(sprintf("Mapped %d chicken -> %d human genes",
                  length(unique(human_map[["9031_ID"]])), length(unique(human_map[["9606_ID"]]))))
}

human_go <- build_go_by_ontology(org.Hs.eg.db, "Human")

# ============================================================
# Run GSEA
# ============================================================
message("=== Running GSEA ===")

run_gsea <- function(stats_vec, pathways, desc_map, ont_map, sym_map,
                     organism, ontology, method, direction) {
  stats_vec <- sort(stats_vec[!is.na(stats_vec)], decreasing=TRUE)
  if (length(stats_vec) < 10) {
    message(sprintf("  %s %s %s %s: too few genes (%d), skipping",
                    organism, ontology, method, direction, length(stats_vec)))
    return(NULL)
  }
  message(sprintf("  %s %s %s %s: %d genes, %d non-zero, range [%.3f, %.3f]",
                  organism, ontology, method, direction, length(stats_vec),
                  sum(stats_vec > 0), min(stats_vec), max(stats_vec)))

  fg <- fgsea(pathways, stats_vec, minSize=10, maxSize=500, scoreType="pos")
  fg <- as.data.frame(fg)

  if (nrow(fg) == 0) return(NULL)

  # Convert leadingEdge list to comma-separated Entrez IDs
  fg$leadingEdge <- sapply(fg$leadingEdge, function(x) paste(as.character(x), collapse=","))

  # Add gene symbols for leadingEdge
  fg$leadingEdgeSymbols <- sapply(fg$leadingEdge, function(le_str) {
    ids <- strsplit(le_str, ",")[[1]]
    syms <- sapply(ids, function(id) {
      s <- sym_map[id]
      ifelse(is.na(s) || is.null(s), id, s)
    })
    paste(syms, collapse=",")
  })

  # Add description and ontology
  fg$Description <- sapply(fg$pathway, function(go_id) {
    d <- desc_map[go_id]
    ifelse(is.na(d), go_id, d)
  })
  fg$ontology <- ontology
  fg$method <- method
  fg$direction <- direction

  fg <- fg[order(fg$padj), ]
  return(fg)
}

# Collect results per organism x ontology
for (organism in c("chicken", "human")) {
  for (ontology in c("BP", "CC", "MF")) {
    all_results <- list()

    for (method in c("count", "binom")) {
      for (direction in c("pos", "neg")) {
        stat_col <- sprintf("stat_%s_%s", method, direction)
        stats_vec <- gene_stats[[stat_col]]
        names(stats_vec) <- as.character(gene_stats$geneId)

        if (organism == "chicken") {
          pathways <- chicken_go$go_list[[ontology]]
          desc_map <- chicken_go$desc_map
          ont_map <- chicken_go$ont_map
          sym_map <- chicken_go$sym_map
        } else {
          # Map chicken gene IDs to human
          if (is.null(human_map) || nrow(human_map) == 0) next
          chicken_ids <- names(stats_vec)
          merged <- merge(data.frame(chicken_id=chicken_ids, stat=stats_vec),
                          data.frame(chicken_id=as.character(human_map[["9031_ID"]]),
                                    human_id=as.character(human_map[["9606_ID"]]),
                                    stringsAsFactors=FALSE),
                          by="chicken_id", all.x=FALSE)
          merged <- merged[!duplicated(merged$human_id), ]
          stats_vec <- merged$stat
          names(stats_vec) <- merged$human_id

          pathways <- human_go$go_list[[ontology]]
          desc_map <- human_go$desc_map
          ont_map <- human_go$ont_map
          sym_map <- human_go$sym_map
        }

        fg <- run_gsea(stats_vec, pathways, desc_map, ont_map, sym_map,
                       organism, ontology, method, direction)
        if (!is.null(fg)) {
          all_results[[length(all_results) + 1]] <- fg
        }
      }
    }

    if (length(all_results) > 0) {
      combined <- do.call(rbind, all_results)
      combined <- combined[order(combined$padj), ]
      out_path <- file.path(OUTPUT_DIR, sprintf("gsea_v3_%s_%s.csv", organism, ontology))
      write.csv(combined, out_path, row.names=FALSE)
      sig <- sum(combined$padj < 0.05, na.rm=TRUE)
      message(sprintf("  -> %s %s: %d terms, %d significant (padj<0.05)",
                      organism, ontology, nrow(combined), sig))
      if (sig > 0) {
        print(head(combined[combined$padj < 0.05,
                           c("pathway","Description","NES","padj","leadingEdgeSymbols")], 5))
      }
    }
  }
}

message("\n=== Step 04 complete ===")
