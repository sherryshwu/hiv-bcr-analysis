#!/usr/bin/env Rscript
# setwd("~/Library/CloudStorage/GoogleDrive-sherry.wu.gr@dartmouth.edu/My Drive/dartmouth/research/hiv")

suppressMessages({
  library(tidyverse)
  library(readr)
  library(purrr)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "data/processed",
              help = "Output directory [default= %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Load all datasets
cat("Loading datasets...\n")
all_clones_data <- read.csv("data/raw/All_clones_Bulkcy16cy19_10xcy16_byFelipe_20250514New Susan Moir.csv")
c02_data <- read.csv("data/raw/c02_g11_sequences/C02_like_C&S_H&L_Ken_071725_jj.csv")
g11_data <- read.csv("data/raw/c02_g11_sequences/G11_like_C&S_H&L_Ken_071725_jj.csv")
tenX_heavy <- read_tsv("data/processed/02_parsed_sequences/filtered_contig_Susan_Moir_heavy_parse-select.tsv")
tenX_light <- read_tsv("data/processed/02_parsed_sequences/filtered_contig_Susan_Moir_light_parse-select.tsv")

cat("=== BCR Dataset Integration Pipeline ===\n")
# --------------------- Helpers ---------------------
clean_sorted_df <- function(df) {
  names(df)[names(df) == "v_identity_.percentage."] <- "v_identity"
  names(df)[names(df) == "j_identity_.percentage."] <- "j_identity"

  df %>%
    # Remove the Locus column and keep the standard locus column
    dplyr::select(-any_of("Locus")) %>%
    # Remove empty or typo strings in sequence_id and locus
    filter(sequence_id != "", sequence_id != "sequence_id", locus != "locus") %>%
    # Step 1: Clean sequence_id and create prefix/suffix
    dplyr::mutate(
      sequence_id = str_remove(sequence_id, "_+$"),
      prefix      = str_extract(sequence_id, "^[^_]+"),
      suffix      = str_extract(sequence_id, "[^_]+$")
    ) %>%
    # Step 2: Use prefix/suffix to create new sequence_id and cell_id
    dplyr::mutate(
      sequence_id = paste(prefix, locus, suffix, sep = "_"),
      cell_id = paste0(prefix, "_", suffix)
    )
}

clean_data <- function(df, source_name) {
  tryCatch({
    cleaned_df <- df %>%
      # Tag source
      dplyr::mutate(data_source = source_name, .before = 1) %>%
      # Normalize factors/logicals
      dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%
      # Explicit character columns commonly seen in AIRR-like tables
      dplyr::mutate(dplyr::across(dplyr::any_of(c(
        "sequence_id", "sequence", "rev_comp", "productive",
        "v_call", "d_call", "j_call", "c_call",
        "sequence_alignment", "germline_alignment",
        "junction", "junction_aa", "cdr3_aa",
        "v_cigar", "d_cigar", "j_cigar",
        "v_quals", "d_quals", "j_quals",
        "locus", "cell_id", "clonal_family", "Source",
        "v_call_10x", "d_call_10x", "j_call_10x",
        "junction_10x", "junction_10x_aa"
      )), as.character)) %>%
      # Convert numeric-ish columns to numeric
      dplyr::mutate(dplyr::across(
        tidyselect::matches(
          "(^|_)(length|start|end|score|identity|support|count|aal)$|^mu_freq$"
        ),
        ~ suppressWarnings(as.numeric(.x))
      )) %>%
      # Coerce np/n lengths across files
      dplyr::mutate(dplyr::across(
        tidyselect::any_of(c("np1_length", "np2_length", "n1_length", "n2_length")),
        ~ suppressWarnings(as.numeric(.x))
      )) %>%
      # Convert remaining logicals to character
      dplyr::mutate(dplyr::across(where(is.logical), as.character))

    cat("✓ Cleaned", source_name, ":", nrow(cleaned_df), "sequences\n")
    return(cleaned_df)
  }, error = function(e) {
    cat("✗ Error cleaning", source_name, ":", e$message, "\n")
    stop(e)
  })
}

# Specific cleaning for C02/G11
c02_data <- clean_sorted_df(c02_data)
g11_data <- clean_sorted_df(g11_data)
all_clones_bulk_data <- all_clones_data %>% filter(type == "Bulk BCR ")

# Prepare bulk data
cat("Processing bulk data...\n")
bulk_std <- all_clones_bulk_data %>%
  clean_data("bulk_bcr") %>%
  mutate(
    cell_id = NA_character_,
    consensus_count = NA_real_,
    umi_count = NA_real_,
    Source = "bulk",
    clonal_family = NA_character_,
    c_call = NA_character_
  )

# Prepare C02 and G11 data
cat("Processing C02 and G11 data...\n")
c02_std <- c02_data %>%
  clean_data("c02_sorted_cultured") %>%
  rename(clonal_family = Clonal.Family) %>%
  mutate(
    time = "Y18",
    consensus_count = NA_real_,
    umi_count = NA_real_,
    v_call_10x = NA_character_,
    j_call_10x = NA_character_,
    junction_10x = NA_character_,
    junction_10x_aa = NA_character_,
    c_call = NA_character_
  )

g11_std <- g11_data %>%
  clean_data("g11_sorted_cultured") %>%
  rename(clonal_family = Clonal.Family) %>%
  mutate(
    time = "Y18",
    consensus_count = NA_real_,
    umi_count = NA_real_,
    v_call_10x = NA_character_,
    j_call_10x = NA_character_,
    junction_10x = NA_character_,
    junction_10x_aa = NA_character_,
    c_call = NA_character_
  )

# Prepare 10X data
cat("Processing 10X heavy chain data...\n")
tenx_heavy_std <- tenX_heavy %>%
  clean_data("10x_heavy") %>%
  mutate(
    time = "Y16",
    Source = "10X_heavy",
    clonal_family = NA_character_,
    clone_id = NA_integer_
  )

cat("Processing 10X light chain data...\n")
tenx_light_std <- tenX_light %>%
  clean_data("10x_light") %>%
  mutate(
    time = "Y16",
    Source = "10X_light",
    clonal_family = NA_character_,
    clone_id = NA_integer_
  )

# Get all unique column names
all_datasets <- list(bulk_std, c02_std, g11_std, tenx_heavy_std, tenx_light_std)
all_columns <- unique(unlist(map(all_datasets, names)))

cat("Total unique columns across datasets:", length(all_columns), "\n")

# Check dimensions of each dataset
cat("Dataset dimensions:\n")
cat("all_clones_bulk_data:", nrow(all_clones_bulk_data), "x", ncol(all_clones_bulk_data), "\n")
cat("c02_data_clean:", nrow(c02_std), "x", ncol(c02_std), "\n")
cat("g11_data_clean:", nrow(g11_std), "x", ncol(g11_std), "\n")
cat("tenX_heavy:", nrow(tenX_heavy), "x", ncol(tenX_heavy), "\n")
cat("tenX_light:", nrow(tenX_light), "x", ncol(tenX_light), "\n")

# Combine all datasets
cat("Combining datasets...\n")
combined_bcr_data <- dplyr::bind_rows(all_datasets) %>% distinct()

# Clean up the combined dataset
cat("Cleaning combined dataset...\n")
combined_clean <- combined_bcr_data %>%
  # Remove rows with no sequence_id
  filter(!is.na(sequence_id) & sequence_id != "") %>%
  # Standardize locus names
  mutate(locus = case_when(
    locus %in% c("IGH", "IgH") ~ "IGH",
    locus %in% c("IGK", "IgK") ~ "IGK",
    locus %in% c("IGL", "IgL") ~ "IGL",
    TRUE ~ locus
  )) %>%
  # Assign cell ID based on data source
  mutate(cell_id = case_when(
    data_source == "bulk_bcr" ~ sequence_id,  # Use sequence_id for bulk data
    data_source %in% c("c02_sorted_cultured", "g11_sorted_cultured") ~ cell_id,  # Keep existing prefix_suffix pattern
    data_source %in% c("10x_heavy", "10x_light") ~ cell_id,  # Keep existing 10X cell_id
    TRUE ~ cell_id
  )) %>%
  # Create unified sample identifier
  mutate(sample_id = case_when(
    data_source == "bulk_bcr" ~ paste0("bulk", "_", coalesce(time, "unknown")),
    data_source == "c02_sorted_cultured" ~ "PT30_sorted_cultured",
    data_source == "g11_sorted_cultured" ~ "PT30_sorted_cultured",
    data_source == "10x_heavy" ~ "10X_heavy",
    data_source == "10x_light" ~ "10X_light",
    TRUE ~ data_source
  )) %>%
  mutate(productive = case_when(
    productive %in% c("TRUE", "True", "true", "T") ~ TRUE,
    productive %in% c("FALSE", "False", "false", "F") ~ FALSE,
    TRUE ~ NA
  )) %>%
  mutate(junction = str_replace_all(junction, "[^ACGTMRWSYKVHDBN.?-]", "N"))

# Print summary statistics
cat("\n=== DATASET SUMMARY ===\n")
cat("Combined dataset dimensions:", nrow(combined_clean), "rows x", ncol(combined_clean), "columns\n")
cat("\nData sources:\n")
print(table(combined_clean$data_source))
cat("\nSample IDs:\n")
print(table(combined_clean$sample_id))
cat("\nLocus distribution:\n")
print(table(combined_clean$locus, useNA = "ifany"))
cat("\nProductivity:\n")
print(table(combined_clean$productive, useNA = "ifany"))

# Save the combined dataset
cat("\nSaving combined dataset...\n")
write_csv(combined_clean, file.path(opt$outdir, "all_combined_bcr_data.csv"))

# Create heavy chain only dataset for clonal clustering
heavy_only <- combined_clean %>%
  filter(locus == "IGH" & productive == TRUE) %>%
  # Remove rows with missing critical fields for clonal clustering
  filter(!is.na(v_call) & !is.na(j_call) & !is.na(junction))

cat("\nHeavy chain productive sequences:", nrow(heavy_only), "\n")
cat("Heavy chain data sources:\n")
print(table(heavy_only$data_source))

# Save heavy chain dataset for clonal clustering
write_tsv(heavy_only, file.path(opt$outdir, "all_combined_heavy_chain_data.tsv"))

cat("\nFiles created:")
cat("\n- all_combined_bcr_data.csv (all data)")
cat("\n- all_combined_heavy_chain_data.tsv (IGH productive for clonal clustering)")
cat("\nData integration complete!\n")