#!/usr/bin/env Rscript
# SUPERSEDED. Publication figure is built by figure_article_panel.R
# (tree + DAR + 10 promoter correlations nearest to TSS).
# This file remains so old run commands still produce the current figure.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
this_dir  <- if (length(this_file)) dirname(normalizePath(this_file)) else getwd()
source(file.path(this_dir, "figure_article_panel.R"), echo = FALSE)
