#!/usr/bin/env Rscript
# Dowser phylogenetic tree building from scoper clones
# ============ SETTINGS ============
suppressMessages({
  library(tidyverse)
  library(dowser)
  library(ggtree)
  library(optparse)
})

# ============ COMMAND LINE OPTIONS ============
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "results/clone_analysis/hl_for_kept_cells_split_light.csv",
              help = "CSV containing sequences to build trees from (e.g., hl_for_kept_cells_* .csv)"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/phylogenetic_trees",
              help = "Output directory [default= %default]"),
  make_option(c("--min_size"), type = "integer", default = 3,
              help = "Minimum clone size [default= %default]"),
  make_option(c("--partition"), type = "character", default = "heavy_only",
              help = "Partition type: 'heavy_only', 'light_only', or 'both' [default= %default]"),
  make_option(c("-n", "--nproc"), type = "integer", default = 1,
              help = "Number of processors [default= %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Decide mode-specific subfolder
mode_outdir <- file.path(opt$outdir, opt$partition)
dir.create(mode_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ============ LOAD DATA ============
cat("\n=== BUILDING PHYLOGENETIC TREES FOR SORTED CELL CLONES ===\n")
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
c02_g11_cell_ids <- file.path("data/processed/combined_datasets/c02_g11_cell_ids.csv")
table(is.na(c02_g11_cell_ids %in% dat$cell_id))

# ============ BUILD TREES ============
cat("=== BUILDING TREES ===\n")
trees <- list()
# Format clones
clones <- formatClones(dat, traits = c("time", "source"), chain = chain_value, minseq = opt$min_size)
saveRDS(clones, file.path(mode_outdir, "target_clones.rds"))

cat("Running IgPhyML via getTrees...\n")
# Define the path to igphyml executable
igphyml_path <- "/dartfs/rc/lab/H/HoehnK/Sherry/igphyml/src/igphyml"

# Print the timestamp
cat("getTrees (partition =", opt$partition, ") started at:", as.character(Sys.time()), "\n")

if (opt$partition == "both") {
  trees <- getTrees(clones, build = "igphyml", nproc = opt$nproc, partition = "hl", exec = igphyml_path)
} else {
  trees <- getTrees(clones, build = "igphyml", nproc = opt$nproc, exec = igphyml_path)
}
cat("getTrees (partition =", opt$partition, ") finished at:", as.character(Sys.time()), "\n")
saveRDS(trees, file.path(mode_outdir, "trees.rds"))
cat("Saved trees RDS to:", file.path(mode_outdir, "trees.rds"), "\n")

# ============ PLOT ============
plots_dir <- file.path(mode_outdir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

cat("\nCreating plots...\n")
p <- list()
for (i in seq_len(nrow(trees))) {
  clone_id <- trees$clone_id[i]
  tree <- trees$trees[[i]]

  # Prepare data for plotting
  plot_data <- trees$data[[i]]@data %>%
    select(sequence_id, time, source)

  # Plot by timepoint
  p[[i]] <- ggtree(tree) %<+% plot_data +
    geom_tippoint(aes(color = time, shape = source), size = 1.5) +
    scale_color_discrete(na.translate = FALSE) +
    scale_shape_discrete(na.translate = FALSE) +
    ggtitle(paste("Clone", "-", opt$partition)) +
    theme(legend.text = element_text(size = 10)) +
    guides(shape = guide_legend(title = "Sequence Source"),
           color = guide_legend(title = "Sample Time"))
}
pdf(file.path(plots_dir, paste0("genetic_distance_trees_", opt$partition, ".pdf")), 
    width = 8, height = 10)
for (i in seq_along(p)) {
  print(p[[i]])
}
dev.off()
cat("\nDone! Check", mode_outdir, "for results\n")

# TODO: compare AA changes from germline to sequences in C02 and G11 by position, AA mutations at the same site among lineages
# remake SHM frequency figure
