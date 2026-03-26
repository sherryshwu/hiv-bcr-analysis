#!/usr/bin/env Rscript
suppressMessages({
  library(tidyverse)
  library(readr)
  library(shazam)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--heavy_data"), type = "character",
              default = "data/processed/03_combined_datasets/all_combined_heavy_chain_data.tsv",
              help = "Input heavy chain TSV [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/clustering_threshold",
              help = "Output directory [default= %default]"),
  make_option(c("-n", "--nproc"), type = "integer", default = 4,
              help = "Number of processors [default= %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Clustering Threshold Determination ===\n")

# Load heavy chain data for threshold calculation
cat("Loading heavy chain data for threshold calculation...\n")
all_heavy_data <- read_tsv(opt$heavy_data, show_col_types = FALSE)

# Remove cells with multiple heavy chains
heavy_data <- all_heavy_data %>%
  group_by(cell_id) %>%
  mutate(n_heavy = sum(locus == "IGH")) %>%
  filter(n_heavy == 1)

cat("Heavy chain sequences for threshold calculation:", nrow(heavy_data), "after removing", nrow(all_heavy_data) - nrow(heavy_data), "cells with multiple heavy chains.\n")

# Find threshold using SHazaM (used by both methods)
cat("Computing distances with SHazaM...\n")
nn <- shazam::distToNearest(
  heavy_data,
  sequenceColumn = "junction",
  vCallColumn = "v_call",
  jCallColumn = "j_call",
  cellIdColumn = "cell_id",
  model = "ham",
  normalize = "len",
  nproc = opt$nproc
)
distances <- nn$dist_nearest[is.finite(nn$dist_nearest)]
threshold <- shazam::findThreshold(distances, method = "density")@threshold

cat("Computed threshold:", round(threshold, 4), "\n")

# Save threshold plot
pdf(file.path(opt$outdir, "threshold_density_plot.pdf"))
plot(density(distances), main = "NN distance density (SHazaM)")
abline(v = threshold, lty = 2, col = "red")
legend("topright", legend = paste("Threshold =", round(threshold, 4)),
       col = "red", lty = 2)
dev.off()

# Save threshold and Change-O input files
write.csv(data.frame(method = "shazam_density", threshold = threshold),
          file.path(opt$outdir, "threshold.csv"), row.names = FALSE)

cat("SHazaM analysis complete! Results saved to:", opt$outdir, "\n")