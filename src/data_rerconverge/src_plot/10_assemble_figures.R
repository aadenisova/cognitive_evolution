#!/usr/bin/env Rscript
# SUPERSEDED. Combined article figure is assembled inside figure_article_panel.R.
# This file remains so old run commands still produce the current figure.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
this_dir  <- if (length(this_file)) dirname(normalizePath(this_file)) else getwd()
source(file.path(this_dir, "figure_article_panel.R"), echo = FALSE)
