suppressMessages({
  library(tidyverse)
  library(readr)
  library(shazam)
  library(scoper)
  library(dowser)
  library(optparse)
})

option_list <- list(
  make_option(c("-a", "--all_data"), type = "character",
              default = "data/processed/03_combined_datasets/all_combined_data.csv",
              help = "Input combined BCR data CSV [default= %default]"),
  make_option(c("-t", "--threshold_file"), type = "character",
              default = "results/clustering_threshold/threshold.csv",
              help = "Input threshold CSV from SHazaM [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/clustering_scoper",
              help = "Output directory [default= %default]"),
  make_option(c("-n", "--nproc"), type = "integer", default = 4,
              help = "Number of processors [default= %default]"),
  make_option(c("-g", "--germline_dir"), type = "character",
              default = "/dartfs/rc/home/5/f0070d5/share/germlines/imgt/human/vdj",
              help = "IMGT germline reference directory [default= %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Load the pre-computed threshold
threshold_df <- tryCatch({
  read_csv(opt$threshold_file, show_col_types = FALSE)
}, error = function(e) {
  stop("Failed to load threshold file: ", opt$threshold_file, "\n", e$message)
})

threshold <- threshold_df$threshold[1]
cat("\n=== SCOPER APPROACH ===\n")
cat("Using threshold:", round(threshold, 4), "for Scoper analysis\n")

# Load all BCR data for Scoper
cat("Loading all BCR data for Scoper...\n")
all_data <- tryCatch({
  read_csv(opt$all_data, show_col_types = FALSE)
}, error = function(e) {
  stop("Failed to load data file: ", opt$all_data, "\n", e$message)
})

# Prepare data for clustering
scoper_data <- all_data %>%
  filter(!is.na(cell_id)) %>%
  filter(productive == TRUE) %>%
  filter(!is.na(v_call), !is.na(j_call), !is.na(junction)) %>%
  filter(locus %in% c("IGH", "IGK", "IGL")) %>%
  distinct() %>%
  # Remove cells with multiple heavy chains
  group_by(cell_id) %>%
  mutate(n_heavy = sum(locus == "IGH")) %>%
  filter(n_heavy == 1) %>%
  mutate(clonal_family_original = clonal_family)

# Count and report non-unique cell_ids
n_non_unique_cells <- scoper_data %>%
  group_by(cell_id) %>%
  summarise(n_sequences = n(), .groups = "drop") %>%
  filter(n_sequences > 1) %>%
  nrow()

# Check the results
final_summary <- scoper_data %>%
  group_by(cell_id) %>%
  summarise(
    n_heavy = sum(locus == "IGH"),
    n_light = sum(locus %in% c("IGK", "IGL")),
    .groups = "drop"
  )

cat("Processed", nrow(scoper_data), "sequences for Scoper analysis\n")
cat("Found", n_non_unique_cells, "cell_ids with multiple sequences\n")
cat("Total cells:", nrow(final_summary), "\n")
cat("Heavy chains:", sum(final_summary$n_heavy), "\n")
cat("Light chains:", sum(final_summary$n_light), "\n")

# Load IMGT references
references <- tryCatch({
  readIMGT(dir = opt$germline_dir)
}, error = function(e) {
  stop("Failed to load IMGT references: ", e$message)
})

# Run Scoper with both split_light settings
for (split_light in c(FALSE, TRUE)) {
  suffix <- ifelse(split_light, "split_light", "no_split_light")
  cat("\n--- Running hierarchicalClones with split_light =", split_light, "---\n")

  start_time <- Sys.time()

  # Print the timestamp
  cat("hierarchicalClones (split_light =", split_light, ") started at:", as.character(start_time), "\n")

  # Run hierarchical clones
  clone_results <- hierarchicalClones(
    scoper_data,
    threshold = threshold,
    split_light = split_light,
    cell_id = "cell_id",
    summarize_clones = FALSE,
    nproc = opt$nproc
  )

  end_time <- Sys.time()

  # Print the timestamps and runtime
  cat("hierarchicalClones (split_light =", split_light, ") finished at:", as.character(end_time), "\n")
  cat("hierarchicalClones (split_light =", split_light, ") runtime:", round(end_time - start_time, 2), "hours\n")

  # Save results with appropriate suffix
  write_csv(clone_results, file.path(opt$outdir, paste0("scoper_clones_", suffix, ".csv")))
# clone_results_split=read_csv(file.path(opt$outdir, paste0("scoper_clones_", suffix, ".csv")))
  # Basic clone statistics
  clone_stats <- clone_results %>%
    group_by(clone_id) %>%
    summarise(
      clone_size = n(),
      datasets = paste(unique(source), collapse = ","),
      n_datasets = length(unique(source)),
      heavy_chains = sum(locus == "IGH"),
      light_chains = sum(locus %in% c("IGK", "IGL")),
      .groups = "drop"
    )

  write_csv(clone_stats, file.path(opt$outdir, paste0("scoper_clone_stats_", suffix, ".csv")))
  cat("Scoper (split_light =", split_light, ") found", length(unique(clone_results$clone_id)), "clones\n")

  # Create germlines with Dowser
  cat("Creating germlines with Dowser (split_light =", split_light, ")...\n")
  clone_id_value <- ifelse(split_light == TRUE, "clone_subgroup_id", "clone_id")

  tryCatch({
    # Clean data for createGermlines
    clone_results_for_germlines <- clone_results %>%
      # Remove NA for clone_id
      filter(!is.na(clone_id)) %>%
      group_by(clone_id) %>%
      # Add locus and row number to sequence_id to make unique
      mutate(
        sequence_id_original = sequence_id,
        sequence_id = paste0(clone_id, "_", locus, "_", row_number())
      ) %>%
      ungroup()

    # Verify uniqueness
    if (any(duplicated(clone_results_for_germlines$sequence_id))) {
      stop("Sequence IDs are not unique after transformation!")
    }

    # Resolve light chains only when split_light = TRUE
    if (split_light == TRUE) {
      cat("Resolving light chains with resolveLightChains() started at:", as.character(Sys.time()), "\n")
      clone_results_for_germlines <- resolveLightChains(clone_results_for_germlines, text_fields = NULL, num_fields = NULL, seq_fields = NULL)
      cat("Resolving light chains with resolveLightChains() ended at:", as.character(Sys.time()), "\n")
      write_csv(clone_results_for_germlines, file.path(opt$outdir, paste0("resolved_scoper_germlines_", suffix, ".csv")))
    }

    germlines <- createGermlines(
      clone_results_for_germlines,
      references = references,
      nproc = opt$nproc,
      clone = clone_id_value,
      trim_lengths = TRUE
    )
    write_csv(germlines, file.path(opt$outdir, paste0("scoper_germlines_", suffix, ".csv")))
  }, error = function(e) {
    cat("Warning: Germline creation failed for split_light =", split_light, ":", e$message, "\n")
  })
}

cat("\nAnalysis complete! Results in:", opt$outdir, "\n")

# See if there's split in clones in split_light=false
# how many light chains / how many cells have a light chain -> should we build tree with light chain
# build tree (without and with light chain - igphyml, light chain partition / dowser vignette) label with time points and data source -> whether sorted cells are grouped under same clonal lineage