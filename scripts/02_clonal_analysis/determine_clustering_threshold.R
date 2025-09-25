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

cat("\n=== THRESHOLD DETERMINATION (SHazaM) ===\n")
# Load heavy chain data for threshold calculation
cat("Loading heavy chain data for threshold calculation...\n")
heavy_data <- read_tsv(opt$heavy_data, show_col_types = FALSE)

cat("Heavy chain sequences for threshold calculation:", nrow(heavy_data), "\n")

# Find threshold using SHazaM (used by both methods)
cat("Computing distances with SHazaM...\n")
nn <- shazam::distToNearest(
  heavy_data,
  sequenceColumn = "junction",
  vCallColumn = "v_call",
  jCallColumn = "j_call",
  model = "ham",
  normalize = "len",
  nproc = opt$nproc
)

d <- nn$dist_nearest[is.finite(nn$dist_nearest)]
threshold <- shazam::findThreshold(d, method = "density")@threshold

cat("Computed threshold:", round(threshold, 4), "\n")

# Save threshold plot
pdf(file.path(opt$outdir, "threshold_density_plot.pdf"))
plot(density(d), main = "NN distance density (SHazaM)")
abline(v = threshold, lty = 2, col = "red")
legend("topright", legend = paste("Threshold =", round(threshold, 4)),
       col = "red", lty = 2)
dev.off()

# Save threshold and Change-O input files
write.csv(data.frame(method = "shazam", threshold = threshold),
          file.path(opt$outdir, "threshold.csv"), row.names = FALSE)

cat("SHazaM analysis complete. Results in:", opt$outdir, "\n")