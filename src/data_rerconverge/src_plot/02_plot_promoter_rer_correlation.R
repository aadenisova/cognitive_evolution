#!/usr/bin/env Rscript
# SUPERSEDED. Promoter RER–innovation scatters are panel C of
# figure_article_panel.R (5 positive + 5 negative, nearest to TSS, ≤ 1 kb).
# This file remains so old run commands still produce the current figure.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
this_dir  <- if (length(this_file)) dirname(normalizePath(this_file)) else getwd()
source(file.path(this_dir, "figure_article_panel.R"), echo = FALSE)
