#!/usr/bin/env Rscript

suppressMessages({
  library(tidyverse)
  library(alakazam)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "results/clustering_scoper/scoper_germlines_no_split_light.csv",
              help = "Input germline data [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "results/mutation_analysis",
              help = "Output directory [default= %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== ANALYZING C02 AND G11 AA MUTATIONS ===\n")

# ==============================================================================
# Load and Filter Data
cat("Loading data from:", opt$input, "\n")
raw_data <- read_csv(opt$input, show_col_types = FALSE)

germline_data <- raw_data %>%
  rename(lineage_id = all_of("Clonal.Family")) %>%
  filter(
    # locus == "IGH",
    source %in% c("cultured", "sorted"),
    productive == TRUE,
    !is.na(sequence_alignment),
    !is.na(germline_alignment),
    # Filter strictly for C02 and G11
    lineage_id %in% c("C02-like", "G11-like")) %>%
  mutate(lineage_id = recode(lineage_id,
                             "C02-like" = "C02",
                             "G11-like" = "G11")
  )

cat("Sequences retained for analysis: ", nrow(germline_data), "\n")
cat("Lineage counts:\n")
print(table(germline_data$lineage_id))

# ==============================================================================
# Translate and Map Mutations
cat("Translating sequences and mapping mutations...\n")

mutation_data <- germline_data %>%
  mutate(
    seq_aa = translateDNA(sequence_alignment, trim = TRUE),
    germ_aa = translateDNA(germline_alignment, trim = TRUE)
  ) %>%
  filter(!is.na(seq_aa), !is.na(germ_aa)) %>%
  rowwise() %>%
  mutate(
    # Get sequence length
    seq_len = min(nchar(seq_aa), nchar(germ_aa)),
    positions = list(1:seq_len),
    seq_chars = list(strsplit(seq_aa, "")[[1]][1:seq_len]),
    germ_chars = list(strsplit(germ_aa, "")[[1]][1:seq_len]),
    is_mutation = list(seq_chars != germ_chars)
  ) %>%
  unnest(cols = c(positions, seq_chars, germ_chars, is_mutation)) %>%
  select(sequence_id, lineage_id, locus, position = positions,
         germline_aa = germ_chars, observed_aa = seq_chars, is_mutation)

# ==============================================================================
# Calculate Frequencies per Lineage
cat("Calculating frequencies...\n")

# Summarize by Lineage AND Position
mutation_freq <- mutation_data %>%
  group_by(lineage_id, locus, position, germline_aa) %>%
  summarise(
    n_sequences = n(),
    n_mutations = sum(is_mutation),
    mutation_freq = n_mutations / n_sequences,
    # Identify dominant mutation
    dominant_aa = {
      mutants <- observed_aa[is_mutation] # subset only mutated AAs
      if (length(mutants) == 0) {
        NA_character_ # Return NA if no mutations
      } else {
        names(sort(table(mutants), decreasing = TRUE))[1] # Return most frequent mutation
      }
    },
    .groups = "drop"
  ) %>%
  mutate(dominant_aa = ifelse(is.na(dominant_aa), germline_aa, dominant_aa))

# Save the raw frequency table
write_csv(mutation_freq, file.path(opt$outdir, "mutation_freq_by_lineage.csv"))

# ==============================================================================
# Visualization: Butterfly Plot
cat("Generating plots...\n")

# Prepare data for plotting
plot_data <- mutation_freq %>%
  group_by(locus, position) %>%
  filter(sum(n_mutations) > 0) %>%
  ungroup() %>%
  mutate(
    # Make G11 negative for the mirror effect
    plot_freq = ifelse(lineage_id == "G11", -mutation_freq, mutation_freq),
    # Create the label text (e.g., "V->A")
    label_text = ifelse(mutation_freq > 0,
                        paste0(germline_aa, "->", dominant_aa),
                        NA),
    # Calculate label Y position
    label_y = ifelse(lineage_id == "G11", plot_freq - 0.05, plot_freq + 0.05)
  )

# Create the Butterfly Plot
p <- ggplot(plot_data, aes(x = position, y = plot_freq, fill = lineage_id)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(aes(y = label_y, label = label_text),
            size = 2.5,
            angle = 90,
            na.rm = TRUE) +
  # Formatting axes
  scale_y_continuous(labels = abs, limits = c(-1.1, 1.1)) +
  scale_x_continuous(breaks = seq(0, max(plot_data$position), by = 5)) +
  scale_fill_manual(values = c("C02" = "#E69F00", "G11" = "#56B4E9")) +
  facet_wrap(~locus, scales = "free_x", ncol = 1) +
  labs(
    title = "Amino acid mutation frequency by position across C02 vs G11 lineages",
    x = "Amino acid position",
    y = "Mutation frequency\n(G11 <--- | ---> C02)",
    fill = "Lineage"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Save Plot
ggsave(file.path(opt$outdir, "mutation_butterfly_plot_split.pdf"), p, width = 12, height = 6)

cat("Done! Results saved to:", opt$outdir, "\n")