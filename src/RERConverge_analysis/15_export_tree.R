#!/usr/bin/env Rscript
# 08_export_tree.R
# Export RER consensus phylogenetic tree (masterTree) rooted at Apteryx mantelli.
# Map tip labels from GCF/GCA accessions to Latin scientific names.
# Save Newick tree + innovation table.

suppressPackageStartupMessages({
  library(ape)
})

OUTPUT_DIR <- "/mnt/results"

# ============================================================
# Load tree
# ============================================================
message("=== Loading masterTree ===")
trees_obj <- readRDS("/workspace/treesObj_cont.rds")
tree <- trees_obj$masterTree
message(sprintf("Tree: %d tips, rooted=%s", Ntip(tree), is.rooted(tree)))

# ============================================================
# Load phenotype data for name mapping
# ============================================================
message("=== Loading phenotype data ===")
pheno <- read.csv("/mnt/user-uploads/df_merged_for_rer_no_na.csv",
                  row.names = 1, stringsAsFactors = FALSE)
message(sprintf("Phenotype: %d species", nrow(pheno)))

# Verify all tips match
stopifnot(all(tree$tip.label %in% rownames(pheno)))

# ============================================================
# Root tree at Apteryx mantelli (GCF_036417845)
# ============================================================
outgroup <- "GCF_036417845"
message(sprintf("Rooting at %s (%s)", outgroup, pheno[outgroup, "ScientificName"]))
tree_rooted <- root(tree, outgroup = outgroup, resolve.root = TRUE)
message(sprintf("Rooted: %s", is.rooted(tree_rooted)))

# ============================================================
# Map tip labels to Latin names
# ============================================================
message("=== Mapping tip labels to Latin names ===")
latin_names <- pheno[tree_rooted$tip.label, "ScientificName"]
# Replace underscores/spaces for Newick compatibility
latin_names_newick <- gsub(" ", "_", latin_names)
tree_rooted$tip.label <- latin_names_newick

# Save Newick tree
tree_path <- file.path(OUTPUT_DIR, "species_tree_rooted.nwk")
write.tree(tree_rooted, file = tree_path)
message(sprintf("Saved tree: %s", tree_path))

# ============================================================
# Save innovation table
# ============================================================
message("=== Saving innovation table ===")
innovation <- data.frame(
  accession = rownames(pheno),
  scientific_name = pheno$ScientificName,
  TECHINNO2025_ResEff = pheno$TECHINNO2025_ResEff,
  TOTALINNOVATIONS2025_ResEff = pheno$TOTALINNOVATIONS2025_ResEff,
  FOODINNO2025_ResEff = pheno$FOODINNO2025_ResEff,
  Relative_brain_size = pheno$Relative_brain_size,
  Wing_length = pheno$Wing_length,
  stringsAsFactors = FALSE
)
innovation <- innovation[order(-innovation$TECHINNO2025_ResEff), ]

innov_path <- file.path(OUTPUT_DIR, "innovation_table.csv")
write.csv(innovation, innov_path, row.names = FALSE)
message(sprintf("Saved innovation table: %s (%d species)", innov_path, nrow(innovation)))

# Print top/bottom 5
message("\nTop 5 most innovative:")
print(head(innovation[, 1:3], 5))
message("\nBottom 5 least innovative:")
print(tail(innovation[, 1:3], 5))

message("\n=== Step 08 complete ===")
