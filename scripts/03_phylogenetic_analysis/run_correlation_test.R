#!/usr/bin/env Rscript
library(tidyverse)
library(dowser)

# Define paths
tree_dir <- "data/processed/05_trees"
out_dir <- "results/correlation_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

partitions <- c("heavy_only", "light_only", "both")
all_results <- list()

for (part in partitions) {
  file_path <- file.path(tree_dir, paste0("trees_", part, ".rds"))

  if (file.exists(file_path)) {
    cat("Processing correlation for:", part, "\n")
    # Get the specific tree object
    current_trees <- readRDS(file_path)

    # Standardize time to weeks (Years * 52)
    current_trees$data <- lapply(current_trees$data, function(x) {
      x@data$time <- as.numeric(gsub("Y", "", x@data$time)) * 52
      return(x)
    })

    res <- correlationTest(current_trees, permutations = 10000)

    # Add the 'partition' column
    res$partition <- part
    all_results[[part]] <- res

    # Write individual CSV files
    write_csv(res, file.path(out_dir, paste0("correlation_test_results_", part, ".csv")))
  } else {
    warning(paste("File missing for partition:", part))
  }
}

# Combine all results into one master dataframe
if (length(all_results) > 0) {
  final_df <- bind_rows(all_results)
  write_csv(final_df, file.path(out_dir, "correlation_test_results_combined.csv"))
  cat("Correlation analysis complete.\n")
}