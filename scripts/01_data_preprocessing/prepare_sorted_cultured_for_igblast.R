#!/usr/bin/env Rscript

suppressMessages({
  library(tidyverse)
  library(Biostrings)
})

# Load the raw C02 and G11 data
cat("Loading C02 and G11 data...\n")
c02_data <- read.csv("data/raw/c02_g11_sequences/C02_like_C&S_H&L_Ken_071725_jj.csv")
g11_data <- read.csv("data/raw/c02_g11_sequences/G11_like_C&S_H&L_Ken_071725_jj.csv")

cat("Raw data loaded:\n")
cat("  C02:", nrow(c02_data), "rows\n")
cat("  G11:", nrow(g11_data), "rows\n")

# Combine and filter out invalid data
combined <- bind_rows(
  c02_data %>% mutate(dataset = "C02"),
  g11_data %>% mutate(dataset = "G11")
)

cat("\nBefore filtering:", nrow(combined), "rows\n")

# Filter out invalid data
combined_clean <- combined %>%
  # Remove rows with missing or empty sequences
  filter(!is.na(sequence), sequence != "", nchar(sequence) > 0) %>%
  # Remove rows with missing or invalid sequence_id
  filter(
    !is.na(sequence_id),
    sequence_id != "",
    sequence_id != "sequence_id",  # Remove header rows
    sequence_id != "locus"         # Remove any other header artifacts
  ) %>%
  # Remove rows with invalid locus
  filter(
    !is.na(locus),
    locus != "",
    locus != "locus"
  )

cat("After filtering:", nrow(combined_clean), "rows\n")
cat("Removed:", nrow(combined) - nrow(combined_clean), "invalid rows\n")

combined_prep <- combined_clean %>%
  mutate(
    # Clean trailing underscores
    sequence_id_clean = str_remove(sequence_id, "_+$"),
    # Create unique ID for IgBLAST
    igblast_id = paste(dataset, sequence_id_clean, row_number(), sep = "_")
  )

cat("\nFinal sequences to annotate:", nrow(combined_prep), "\n")
cat("By dataset:\n")
print(table(combined_prep$dataset))
cat("By locus:\n")
print(table(combined_prep$locus))
cat("By source:\n")
print(table(combined_prep$Source))

# Validate sequences before writing
cat("\nValidating sequences...\n")
invalid_seqs <- combined_prep %>%
  filter(nchar(sequence) == 0 | is.na(sequence))

if (nrow(invalid_seqs) > 0) {
  cat("ERROR: Found", nrow(invalid_seqs), "invalid sequences after filtering!\n")
  print(invalid_seqs %>% select(sequence_id, sequence))
  stop("Invalid sequences detected. Please review filtering logic.")
}

# Check sequence length distribution
seq_lengths <- nchar(combined_prep$sequence)
cat("Sequence length summary:\n")
summary(seq_lengths)

# Create output directory
dir.create("data/processed/00_sorted_cultured_prep", recursive = TRUE, showWarnings = FALSE)

# Write FASTA file for IgBLAST
cat("\nWriting FASTA file...\n")
fasta_sequences <- DNAStringSet(combined_prep$sequence)
names(fasta_sequences) <- combined_prep$igblast_id

writeXStringSet(
  fasta_sequences,
  filepath = "data/processed/00_sorted_cultured_prep/sorted_cultured_sequences.fasta",
  format = "fasta"
)

# Verify FASTA was written correctly
fasta_check <- readDNAStringSet("data/processed/00_sorted_cultured_prep/sorted_cultured_sequences.fasta")
cat("FASTA verification: wrote", length(fasta_check), "sequences\n")

if (length(fasta_check) != nrow(combined_prep)) {
  stop("ERROR: FASTA file has different number of sequences than expected!")
}

# Save complete mapping file
cat("\nSaving mapping file with all original columns...\n")
mapping <- combined_prep %>%
  select(
    igblast_id,
    sequence_id_original = sequence_id,
    sequence_id_clean,
    dataset,
    locus_original = locus,
    Source,
    everything(),
    -sequence
  )

write_csv(mapping, "data/processed/00_sorted_cultured_prep/sorted_cultured_mapping.csv")

# Summary report
cat("\n=== PREPARATION COMPLETE ===\n")
cat("Input:", nrow(combined), "sequences\n")
cat("Filtered out:", nrow(combined) - nrow(combined_prep), "invalid sequences\n")
cat("Output:", nrow(combined_prep), "sequences\n")
cat("\nFiles created:\n")
cat("  - FASTA: data/processed/00_sorted_cultured_prep/sorted_cultured_sequences.fasta\n")
cat("  - Mapping: data/processed/00_sorted_cultured_prep/sorted_cultured_mapping.csv\n")
cat("\nBreakdown of removed sequences:\n")
cat("  - Empty/NA sequences:", sum(is.na(combined$sequence) | combined$sequence == ""), "\n")
cat("  - Invalid sequence_id:", sum(is.na(combined$sequence_id) | combined$sequence_id == "" | combined$sequence_id == "sequence_id"), "\n")
cat("\nReady for IgBLAST annotation!\n")