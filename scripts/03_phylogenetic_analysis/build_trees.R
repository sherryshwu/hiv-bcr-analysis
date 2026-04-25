#!/usr/bin/env Rscript
# Dowser phylogenetic tree building from scoper clones
# ============ SETTINGS ============
suppressMessages({
  library(tidyverse)
  library(dowser)
  library(ggtree)
  library(optparse)
  library(readxl)
})

# ============ COMMAND LINE OPTIONS ============
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "results/clone_analysis/targeted_clones_split_light.csv",
              help = "CSV containing sequences to build trees from (e.g., targeted_clones_* .csv)"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/phylogenetic_trees",
              help = "Output directory [default= %default]"),
  make_option(c("--min_size"), type = "integer", default = 3,
              help = "Minimum clone size [default= %default]"),
  make_option(c("--partition"), type = "character", default = "both",
              help = "Partition type: 'heavy_only', 'light_only', or 'both' [default= %default]"),
  make_option(c("-n", "--nproc"), type = "integer", default = 1,
              help = "Number of processors [default= %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ============ LOAD DATA ============
cat("\n=== BUILDING PHYLOGENETIC TREES FOR SORTED CLONES ===\n")
cat("Input:", opt$input, "\n")
cat("Output:", opt$outdir, "\n")
cat("Min clone size:", opt$min_size, "\n")
cat("Include light chains:", opt$partition, "\n\n")

# Load clone data
cat("Loading clone data...\n")
dat <- read_csv(opt$input, show_col_types = FALSE) %>%
  filter(productive == TRUE,
         !is.na(sequence_alignment),
         !is.na(germline_alignment))

# ============ FILTER CLONES ============
cat("Filtering data by locus...\n")
# Filter by locus based on partition type
if (opt$partition == "heavy_only") {
  dat <- dat %>% filter(locus == "IGH")
  chain_value <- "H"
} else if (opt$partition == "light_only") {
  dat <- dat %>% filter(locus %in% c("IGK", "IGL"))
  chain_value <- "L"
} else if (opt$partition == "both") {
  dat <- dat %>% filter(locus %in% c("IGH", "IGK", "IGL"))
  chain_value <- "HL"
} else {
  stop("Partition must be one of: heavy_only, light_only, both")
}

cat("Filtered to", nrow(dat), "sequences\n")

# Verify all cell IDs from C02 & G11 data are present in the final clone
c02_g11_cell_ids <- read_csv(file.path("data/processed/03_combined_datasets/c02_g11_cell_ids.csv")) %>% pull(cell_id)
table(!is.na(c02_g11_cell_ids %in% dat$cell_id))

# ============ BUILD TREES ============
cat("=== BUILDING TREES ===\n")
trees <- list()
# Format clones
clones <- formatClones(dat, traits = c("time", "source", "cell_id"), split_light = FALSE, chain = chain_value, minseq = opt$min_size)
clones_dir <- file.path("data/processed/04_target_clones")
dir.create(clones_dir, showWarnings = FALSE)
saveRDS(clones, file.path(clones_dir, paste0("target_clones_", opt$partition, ".rds")))

cat("Running IgPhyML via getTrees...\n")
# Define the path to igphyml executable
igphyml_path <- "/dartfs/rc/lab/H/HoehnK/Sherry/igphyml/src/igphyml"

# Print the timestamp
cat("getTrees (partition =", opt$partition, ") started at:", as.character(Sys.time()), "\n")

if (opt$partition == "both") {
  trees <- getTrees(clones, build = "igphyml", nproc = opt$nproc, partition = "hl", optimize = "tlr", exec = igphyml_path)
} else {
  trees <- getTrees(clones, build = "igphyml", nproc = opt$nproc, optimize = "tlr", exec = igphyml_path)
}
cat("getTrees (partition =", opt$partition, ") finished at:", as.character(Sys.time()), "\n")
trees_dir <- file.path("data/processed/05_trees")
dir.create(trees_dir, showWarnings = FALSE)
saveRDS(trees, file.path(trees_dir, paste0("trees", "_", opt$partition, ".rds")))
cat("Saved trees RDS to:", file.path(trees_dir, paste0("trees", "_", opt$partition, ".rds")))

# ============ PLOT ============
cat("\nCreating plots...\n")
# Define custom shapes
my_shapes <- c("10x" = 21, # Filled circle
               "bulk" = 22, # Filled square
               "sorted" = 25, # Open inverted triangle
               "cultured" = 24) # Open triangle

targets_df <- read_excel("data/raw/Sequences for Trees 012726.xlsx")
target_patterns <- targets_df$`Original Sequence ID` %>%
  na.omit() %>%
  trimws() %>%
  unique()
fuzzy_patterns <- gsub("_", ".*", target_patterns)
search_regex <- paste(fuzzy_patterns, collapse = "|")

p <- list()
for (i in seq_len(nrow(trees))) {
  clone_id <- trees$clone_id[i]
  tree <- trees$trees[[i]]

  # Prepare data for plotting
  plot_data <- trees$data[[i]]@data %>%
    select(sequence_id, time, source)

  # Plot by timepoint
  p[[i]] <- ggtree(tree) %<+% plot_data +
    geom_treescale(x = 0, y = -8, width = 0.05, fontsize = 3, linesize = 0.5) +
    geom_tippoint(aes(fill = time, shape = source), size = 1.5) +
    geom_tiplab(aes(label = label),
                data = td_filter(grepl(search_regex, label)),
                size = 3,
                hjust = -0.1,
                offset = 0.001) +
    scale_color_discrete(na.translate = FALSE) +
    scale_fill_discrete(na.translate = FALSE) +
    scale_shape_manual(values = my_shapes, na.translate = FALSE) +
    ggtitle(paste("Clone", "-", trees$clone_id[i])) +
    theme(legend.text = element_text(size = 10)) +
    guides(shape = guide_legend(title = "source"),
           fill = guide_legend(title = "sample time", override.aes = list(shape = 21, color = "transparent")))
}
pdf(file.path(opt$outdir, paste0("genetic_distance_trees_with_label_", opt$partition, ".pdf")),
    width = 8, height = 10)
for (i in seq_along(p)) {
  print(p[[i]])
}
dev.off()
cat("\nDone! Check", opt$outdir, "for results\n")