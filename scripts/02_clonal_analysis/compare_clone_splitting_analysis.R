#!/usr/bin/env Rscript
suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/clustering_scoper",
              help = "Scoper results directory [default= %default]"),
  make_option(c("--analysis_outdir"), type = "character",
              default = "results/comparative_analysis",
              help = "Output directory for analysis results [default= %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$analysis_outdir, recursive = TRUE, showWarnings = FALSE)

#--------------- Anlayze light chain distribution ---------------#
# Load both results with and without split_light
no_split_data <- read_csv(file.path(opt$outdir, "scoper_clones_no_split_light.csv"))
split_data <- read_csv(file.path(opt$outdir, "scoper_clones_split_light.csv"))

# Load clone statistics
no_split_stats <- read_csv(file.path(opt$outdir, "scoper_clone_stats_no_split_light.csv"))
split_stats <- read_csv(file.path(opt$outdir, "scoper_clone_stats_split_light.csv"))

# Basic comparison
cat("No split_light: ", length(unique(no_split_data$clone_id)), "clones\n")
cat("With split_light: ", length(unique(split_data$clone_id)), "clones\n")

# Analyze light chain distribution
light_chain_analysis <- function(df, suffix) {
  stats <- df %>%
    summarise(
      total_cells = n(),
      cells_with_heavy = sum(locus == "IGH"),
      cells_with_light = sum(locus %in% c("IGK", "IGL")),
      cells_with_kappa = sum(locus == "IGK"),
      cells_with_lambda = sum(locus == "IGL"),
      unique_cells = length(unique(cell_id_unique))
    ) %>%
    mutate(
      pct_cells_with_light = round(cells_with_light / total_cells * 100, 2),
      light_chain_ratio = round(cells_with_kappa / cells_with_lambda, 2),
    )

  cat("\n --- Light Chain Analysis (", suffix, ") ---\n", sep = "")
  print(stats, width = Inf)
  return(stats)
}

no_split_lc <- light_chain_analysis(no_split_data, "No split_light")
split_lc <- light_chain_analysis(split_data, "split_light")

#--------------- Compare clonal assignment ---------------#
# Get cell-level assignments
no_split_assignments <- no_split_data %>%
  filter(locus == "IGH") %>%
  select(cell_id_unique, clone_id, data_source, time) %>%
  mutate(clone_id = as.character(clone_id)) %>%
  distinct() %>%
  rename(clone_id_no_split = clone_id)

split_assignments <- split_data %>%
  filter(locus == "IGH") %>%
  select(cell_id_unique, clone_id, data_source, time) %>%
  mutate(clone_id = as.character(clone_id)) %>%
  distinct() %>%
  rename(clone_id_split = clone_id)

# Merge assignments
clone_comparison <- full_join(no_split_assignments, split_assignments, by = c("cell_id_unique"))
nrow(clone_comparison)

#--------------- Analyze Clone Splitting Patterns ---------------#
# What percentage of original clones (split_light = FALSE) got split?
# Analyze how original clones got split up
clone_splitting_analysis <- clone_comparison %>%
  group_by(clone_id_no_split) %>%
  summarise(
    original_clone_size = n(),
    n_split_clones = length(unique(clone_id_split)),
    split_clone_ids = paste(unique(clone_id_split), collapse = ","),
    data_sources = paste(unique(data_source.x), collapse = ","),
    timepoints = paste(unique(time.x), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    data_sources = factor(data_sources),
    timepoints = factor(timepoints),
    was_split = n_split_clones > 1
  ) %>%
  arrange(desc(n_split_clones))
table(clone_splitting_analysis$n_split_clones)
print(head(clone_splitting_analysis, 10), width = Inf)

# Overall splitting statistics
print(summary(clone_splitting_analysis, width = Inf))
splitting_summary <- clone_splitting_analysis %>%
  summarise(
    total_original_clones = n(),
    clones_that_split = sum(was_split),
    pct_clones_split = round(100 * sum(was_split) / n(), 2),
    avg_splits_per_clone = round(mean(n_split_clones), 2),
    max_splits = max(n_split_clones)
  )

cat("\n=== CLONE SPLITTING SUMMARY ===\n")
print(splitting_summary, width = Inf)

#--------------- Clonal assignments for sorted cells ---------------#
# Which clones are sorted cells assigned to?
cat("\n=== SORTED CELL CLONE ASSIGNMENTS ===\n")
# Extract sorted cell clone IDs
sorted_clone_ids_split <- split_data %>%
  filter(Source == "PT30_sorted cell") %>%
  pull(clone_id) %>%
  unique() %>%
  sort
sorted_clone_ids_no_split <- no_split_data %>%
  filter(Source == "PT30_sorted cell") %>%
  pull(clone_id) %>%
  unique() %>%
  sort

cat("Sorted cells (C02 and G11) clone assignments (with split_light):\n")
print(sorted_clone_ids_split)
cat("Sorted cells (C02 and G11) clone assignments (without split_light):\n")
print(sorted_clone_ids_no_split)

cat("Number of clones with sorted cells (split_light):", length(sorted_clone_ids_split), "\n")
cat("Number of clones with sorted cells (no split_light):", length(sorted_clone_ids_no_split), "\n")

g11_clone_id_comparison <- clone_comparison %>%
  filter(data_source.x == "g11_sorted_cultured") %>%
  select(clone_id_no_split, clone_id_split) %>%
  distinct
c02_clone_id_comparison <- clone_comparison %>%
  filter(data_source.x == "c02_sorted_cultured") %>%
  select(clone_id_no_split, clone_id_split) %>%
  distinct

cat("\nSorted and cultured cells (G11) clone assignments: \n")
print(g11_clone_id_comparison)
cat("\nSorted and cultured cells (C02) clone assignments: \n")
print(c02_clone_id_comparison)

# ============ SAVE SORTED CLONE IDs FOR TREE BUILDING ============
cat("\n=== SAVING CLONE IDs FOR TREE BUILDING ===\n")

# Save as simple text files
writeLines(as.character(sorted_clone_ids_no_split),
           file.path(opt$analysis_outdir, "sorted_clone_ids_no_split.txt"))
writeLines(as.character(sorted_clone_ids_split),
           file.path(opt$analysis_outdir, "sorted_clone_ids_split.txt"))

# Save as CSV with metadata
sorted_clones_no_split_info <- no_split_data %>%
  filter(clone_id %in% sorted_clone_ids_no_split) %>%
  group_by(clone_id) %>%
  summarise(
    clone_size = n(),
    n_heavy = sum(locus == "IGH"),
    n_light = sum(locus %in% c("IGK", "IGL")),
    data_sources = paste(unique(data_source), collapse = ","),
    timepoints = paste(unique(time), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(clone_size))

sorted_clones_split_info <- split_data %>%
  filter(clone_id %in% sorted_clone_ids_split) %>%
  group_by(clone_id) %>%
  summarise(
    clone_size = n(),
    n_heavy = sum(locus == "IGH"),
    n_light = sum(locus %in% c("IGK", "IGL")),
    data_sources = paste(unique(data_source), collapse = ","),
    timepoints = paste(unique(time), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(clone_size))

write_csv(sorted_clones_no_split_info,
          file.path(opt$analysis_outdir, "sorted_clones_no_split_info.csv"))
write_csv(sorted_clones_split_info,
          file.path(opt$analysis_outdir, "sorted_clones_split_info.csv"))

cat("Saved clone IDs to:\n")
cat("  -", file.path(opt$analysis_outdir, "sorted_clone_ids_no_split.txt"), "\n")
cat("  -", file.path(opt$analysis_outdir, "sorted_clone_ids_split.txt"), "\n")
cat("  -", file.path(opt$analysis_outdir, "sorted_clones_no_split_info.csv"), "\n")
cat("  -", file.path(opt$analysis_outdir, "sorted_clones_split_info.csv"), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")