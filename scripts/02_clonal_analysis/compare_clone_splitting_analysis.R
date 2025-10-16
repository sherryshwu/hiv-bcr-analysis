#!/usr/bin/env Rscript
suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--indir"), type = "character",
              default = "results/clustering_scoper",
              help = "Scoper results directory [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/clone_analysis")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

#--------------- Anlayze light chain distribution ---------------#
# Load both results with and without split_light
no_split_data <- read_csv(file.path(opt$indir, "scoper_germlines_no_split_light.csv"))
split_data <- read_csv(file.path(opt$indir, "scoper_germlines_split_light.csv"))

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
      unique_cells = length(unique(cell_id))
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
# Get sequence-level assignments
no_split_assignments <- no_split_data %>%
  select(cell_id, locus, clone_id, source, time) %>%
  mutate(clone_id = as.character(clone_id)) %>%
  distinct() %>%
  rename(clone_id_no_split = clone_id)

split_assignments <- split_data %>%
  select(cell_id, locus, clone_id, source, time) %>%
  mutate(clone_id = as.character(clone_id)) %>%
  distinct() %>%
  rename(clone_id_split = clone_id)

# Merge assignments
clone_comparison <- full_join(no_split_assignments, split_assignments, by = c("cell_id", "locus"))

#--------------- Analyze Clone Splitting Patterns ---------------#
# What percentage of original clones (split_light = FALSE) got split?
# Analyze how original clones got split up
clone_splitting_analysis <- clone_comparison %>%
  group_by(clone_id_no_split) %>%
  summarise(
    original_clone_size = n(),
    n_split_clones = length(unique(clone_id_split)),
    split_clone_ids = paste(unique(clone_id_split), collapse = ","),
    sources = paste(unique(source.x), collapse = ","),
    timepoints = paste(unique(time.x), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    sources = factor(sources),
    timepoints = factor(timepoints),
    was_split = n_split_clones > 1
  ) %>%
  arrange(desc(n_split_clones))

cat("\n--- Number of clones split into ---\n")
table(clone_splitting_analysis$n_split_clones)
print(head(clone_splitting_analysis, 10), width = Inf)

# Overall splitting statistics
cat("\n=== CLONE SPLITTING SUMMARY ===\n")
print(summary(clone_splitting_analysis, width = Inf))
splitting_summary <- clone_splitting_analysis %>%
  summarise(
    total_original_clones = n(),
    clones_that_split = sum(was_split),
    pct_clones_split = round(100 * sum(was_split) / n(), 2),
    avg_splits_per_clone = round(mean(n_split_clones), 2),
    max_splits = max(n_split_clones)
  )
print(splitting_summary, width = Inf)

#--------------- Clonal assignments for sorted cells ---------------#
# Which clones are sorted cells assigned to?
cat("\n=== SORTED CELL CLONE ASSIGNMENTS ===\n")
extract_sorted_clones <- function(df, mode_label, outdir = "results/clone_analysis") {
  # Find heavy chain cells in cultured and sorted samples
  heavy_keep <- df %>%
    filter(source %in% c("cultured", "sorted"), locus == "IGH")

  # Extract cultured and sorted cell IDs
  targeted_cell_ids <- heavy_keep %>% distinct(cell_id)

  # Get BOTH heavy + light for those cells
  hl_for_kept_cells <- df %>%
    semi_join(targeted_cell_ids, by = "cell_id") %>%
    filter(source %in% c("cultured", "sorted"),
           locus %in% c("IGH", "IGK", "IGL"))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_csv(hl_for_kept_cells,
            file.path(outdir, paste0("hl_for_kept_cells_", mode_label, ".csv")))

  cat("Saved: hl_for_kept_cells_", mode_label)
  invisible(list(hl_for_kept_cells = hl_for_kept_cells,
                 targeted_cell_ids = targeted_cell_ids))
}

# Extract HL inputs (no_split & split_light)
res_no_split <- extract_sorted_clones(df = no_split_data, mode_label = "no_split_light")
res_split <- extract_sorted_clones(df = split_data, mode_label = "split_light")

# Extract cultured and sorted cell IDs
targeted_cell_ids_split <- res_split$targeted_cell_ids %>% pull(cell_id)
targeted_cell_ids_no_split <- res_no_split$targeted_cell_ids %>% pull(cell_id)

cat("Sorted cells (C02 and G11) clone assignments (without split_light):\n")
print(targeted_cell_ids_no_split)
cat("Sorted cells (C02 and G11) clone assignments (with split_light):\n")
print(targeted_cell_ids_split)

cat("Number of sorted cells (split_light):", length(targeted_cell_ids_split), "\n")
cat("Number of sorted cells (no split_light):", length(targeted_cell_ids_no_split), "\n")

# ============ SAVE SORTED CELL IDs FOR TREE BUILDING ============
cat("\n=== SAVING CELL IDs FOR TREE BUILDING ===\n")

# Save as simple text files
writeLines(as.character(targeted_cell_ids_no_split),
           file.path(opt$outdir, "targeted_cell_ids_no_split_light.txt"))
writeLines(as.character(targeted_cell_ids_split),
           file.path(opt$outdir, "targeted_cell_ids_split_light.txt"))

# Save as CSV with metadata
targeted_cells_no_split_info <- no_split_data %>%
  filter(cell_id %in% targeted_cell_ids_no_split) %>%
  group_by(clone_id) %>%
  summarise(
    clone_size = n(),
    n_heavy = sum(locus == "IGH"),
    n_light = sum(locus %in% c("IGK", "IGL")),
    sources = paste(unique(source), collapse = ","),
    timepoints = paste(unique(time), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(clone_size))

targeted_cells_split_info <- split_data %>%
  filter(cell_id %in% targeted_cell_ids_split) %>%
  group_by(clone_id) %>%
  summarise(
    clone_size = n(),
    n_heavy = sum(locus == "IGH"),
    n_light = sum(locus %in% c("IGK", "IGL")),
    sources = paste(unique(source), collapse = ","),
    timepoints = paste(unique(time), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(clone_size))

write_csv(targeted_cells_no_split_info,
          file.path(opt$outdir, "targeted_cells_no_split_info.csv"))
write_csv(targeted_cells_split_info,
          file.path(opt$outdir, "targeted_cells_split_info.csv"))

cat("Saved clone IDs to:\n")
cat("  -", file.path(opt$outdir, "targeted_cell_ids_no_split.txt"), "\n")
cat("  -", file.path(opt$outdir, "targeted_cell_ids_split.txt"), "\n")
cat("  -", file.path(opt$outdir, "targeted_cells_no_split_info.csv"), "\n")
cat("  -", file.path(opt$outdir, "targeted_cells_split_info.csv"), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")