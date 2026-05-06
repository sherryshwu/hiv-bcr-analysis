#!/usr/bin/env Rscript
suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(optparse)
  library(airr)
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

#--------------- Clonal assignments for cultured and sorted cells ---------------#
# Which clones are cultured and sorted cells assigned to?
split_data <- read_csv(file.path(opt$indir, "scoper_germlines_split_light.tsv"), show_col_types = FALSE)
cat("\n=== CULTURED AND SORTED CELL CLONE ASSIGNMENTS ===\n")
extract_targeted_clones <- function(df, mode_label, outdir) {
  # Identify clone IDs containing cultured and sorted samples
  targeted_clone_ids <- df %>%
    filter(source %in% c("cultured", "sorted"), locus == "IGH") %>%
    distinct(clone_id)

  # Extract all sequences (heavy + light) belonging to those clones
  targeted_clones <- df %>%
    semi_join(targeted_clone_ids, by = "clone_id") %>%
    filter(locus %in% c("IGH", "IGK", "IGL"))

  # Save full sequence data
  output_file <- file.path(outdir, paste0("targeted_clones_", mode_label, ".csv"))
  write_csv(targeted_clones, output_file)

  cat("Saved full sequences to:", output_file, "\n")
  invisible(list(data = targeted_clones, clone_ids = targeted_clone_ids))
}

# Extract inputs
res <- extract_targeted_clones(df = split_data, mode_label = "split_light", outdir = opt$outdir)
targeted_clone_ids <- res$clone_ids %>% pull(clone_id)

cat("Targeted clones (C02 and G11) clone assignments (with split_light):\n")
print(targeted_clone_ids)

cat("Number of targeted clones (with split_light):", length(targeted_clone_ids), "\n")

# ============ SAVE TARGETED CLONE IDs FOR TREE BUILDING ============
cat("\n=== SAVING CLONE IDs FOR TREE BUILDING ===\n")

# Save as simple text files
writeLines(as.character(targeted_clone_ids), file.path(opt$outdir, "targeted_clone_ids_split_light.txt"))

# Save as CSV with metadata
summary_info <- split_data %>%
  filter(clone_id %in% targeted_clone_ids) %>%
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

write_csv(summary_info, file.path(opt$outdir, "targeted_clones_summary_metadata.csv"))

cat("Saved clone IDs to:\n")
cat("  -", file.path(opt$outdir, "targeted_clone_ids_split_light.txt"), "\n")

# Quick reference for V(D)J calls
gene_summary <- res$data %>%
  select(clone_id, any_of(c("clone_subgroup_id", "locus", "v_call", "d_call", "j_call", "junction_length"))) %>% 
  distinct()

write_csv(gene_summary, file.path(opt$outdir, "targeted_clones_gene_calls.csv"))

cat("\n=== ANALYSIS COMPLETE ===\n")