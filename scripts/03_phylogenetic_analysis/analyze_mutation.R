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

# Load and filter data
germline_data <- read_csv(opt$input, show_col_types = FALSE) %>%
  filter(source %in% c("cultured", "sorted"),
         productive == TRUE,
         !is.na(sequence_alignment),
         !is.na(germline_alignment))

cat("Analyzing", nrow(germline_data), "C02/G11 sequences\n")

# Translate sequences and calculate mutations by position
mutation_data <- germline_data %>%
  mutate(
    seq_aa = translateDNA(sequence_alignment, trim = TRUE),
    germ_aa = translateDNA(germline_alignment, trim = TRUE)
  ) %>%
  filter(!is.na(seq_aa), !is.na(germ_aa)) %>%
  rowwise() %>%
  mutate(
    # Calculate mutations at each position
    positions = list(1:min(nchar(seq_aa), nchar(germ_aa))),
    seq_chars = list(strsplit(seq_aa, "")[[1]]),
    germ_chars = list(strsplit(germ_aa, "")[[1]]),
    mutations = list(seq_chars != germ_chars)
  ) %>%
  unnest(cols = c(positions, seq_chars, germ_chars, mutations)) %>%
  select(sequence_id, clonal_family, clone_id, source,
         position = positions,
         germline_aa = germ_chars,
         observed_aa = seq_chars,
         is_mutation = mutations)

# Calculate mutation frequency by position
mutation_freq <- mutation_data %>%
  group_by(position) %>%
  summarise(
    n_sequences = n(),
    n_mutations = sum(is_mutation),
    mutation_freq = n_mutations / n_sequences,
    .groups = "drop"
  )

write_csv(mutation_freq, file.path(opt$outdir, "mutation_freq_by_position.csv"))

# Define CDR / FR boundaries (example values; adjust based on IMGT numbering)
cdr_regions <- data.frame(
  region = c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3"),
  start = c(1, 27, 38, 56, 66, 105),
  end   = c(26, 37, 55, 65, 104, 130)
) %>%
  mutate(midpoint = (start + end) / 2)

# Plot mutation frequency
p <- ggplot(mutation_freq, aes(x = position, y = mutation_freq)) +
  # shaded regions
  geom_rect(data = cdr_regions,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.08, inherit.aes = FALSE, color = NA) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 1.5, alpha = 0.6) +
  labs(
    title = "AA Mutation Frequency by Position (C02 & G11)",
    x = "Amino Acid Position",
    y = "Mutation Frequency"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank())

# Add vertical dashed lines at region boundaries
for (xpos in unique(c(cdr_regions$start, cdr_regions$end))) {
  p <- p + geom_vline(xintercept = xpos, linetype = "dashed", color = "grey60", linewidth = 0.3)
}

# Add region labels at midpoints
p <- p +
  geom_text(
    data = cdr_regions,
    aes(x = midpoint, y = 1.05 * max(mutation_freq$mutation_freq), label = region),
    color = "black", size = 4, vjust = 0
  )

ggsave(file.path(opt$outdir, "mutation_freq_by_position.pdf"),
       p, width = 10, height = 6)

# Calculate mutation frequency by clonal family
mutation_freq_family <- mutation_data %>%
  mutate(clonal_family = case_when(
    clone_id == 36344 ~ "C02",
    clone_id == 36362 ~ "G11",
    TRUE ~ "Other"
  )) %>%
  group_by(clonal_family, position) %>%
  summarise(mutation_freq = mean(is_mutation), .groups = "drop")

p2 <- ggplot(mutation_freq_family, aes(x = position, y = mutation_freq, color = factor(clonal_family))) +
  # shaded regions
  geom_rect(data = cdr_regions,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.08, inherit.aes = FALSE, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("C02" = "steelblue", "G11" = "darkorange"), name = "Clonal family") +
  labs(title = "AA Mutation Frequency by Clonal Family",
       x = "Position", y = "Mutation Frequency") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# Add vertical dashed lines at region boundaries
for (xpos in unique(c(cdr_regions$start, cdr_regions$end))) {
  p2 <- p2 + geom_vline(xintercept = xpos, linetype = "dashed", color = "grey60", linewidth = 0.3)
}

# Add region labels at midpoints
p2 <- p2 +
  geom_text(
    data = cdr_regions,
    aes(x = midpoint, y = 1.05 * max(mutation_freq$mutation_freq), label = region),
    color = "black", size = 4, vjust = 0
  )

ggsave(file.path(opt$outdir, "mutation_freq_by_family.pdf"),
       p2, width = 10, height = 6)

cat("\nResults saved to:", opt$outdir, "\n")
