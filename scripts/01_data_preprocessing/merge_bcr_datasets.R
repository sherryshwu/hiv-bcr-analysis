#!/usr/bin/env Rscript
# setwd("~/Library/CloudStorage/GoogleDrive-sherry.wu.gr@dartmouth.edu/My Drive/dartmouth/research/hiv")

suppressMessages({
  library(tidyverse)
  library(readr)
  library(purrr)
  library(alakazam)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "data/processed/combined_datasets",
              help = "Output directory [default= %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Load all datasets
cat("Loading datasets...\n")
bulk_data <- read.csv("data/raw/All_clones_Bulkcy16cy19_10xcy16_byFelipe_20250514New Susan Moir.csv")
c02_data <- read.csv("data/raw/c02_g11_sequences/C02_like_C&S_H&L_Ken_071725_jj.csv")
g11_data <- read.csv("data/raw/c02_g11_sequences/G11_like_C&S_H&L_Ken_071725_jj.csv")
tenx_heavy <- read_tsv("data/processed/02_parsed_sequences/filtered_contig_Susan_Moir_heavy_parse-select.tsv")
tenx_light <- read_tsv("data/processed/02_parsed_sequences/filtered_contig_Susan_Moir_light_parse-select.tsv")

cat("=== BCR Dataset Integration Pipeline ===\n")
# --------------------- Helpers ---------------------
clean_sorted_cultured <- function(df) {
  df %>%
    # Preserve the original Source column with sort and culture info
    rename(source = Source) %>%
    mutate(source = case_when(
      source == "PT30_culture" ~ "cultured",
      source == "PT30_sorted cell" ~ "sorted",
      TRUE ~ source
    )) %>%
    # Remove the Locus column and keep the standard locus column
    dplyr::select(-any_of("Locus")) %>%
    # Remove empty or typo strings in sequence_id and locus
    filter(sequence_id != "", sequence_id != "sequence_id", locus != "locus") %>%
    rename(clonal_family = Clonal.Family,
           v_identity = v_identity_.percentage.,
           j_identity = j_identity_.percentage.) %>%
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
    ) %>%
    # Step 3: Remove "Homsap" prefix
    dplyr::mutate(
      v_call = getGene(v_call, first = FALSE, strip_d = FALSE),
      d_call = getGene(d_call, first = FALSE, strip_d = FALSE),
      j_call = getGene(j_call, first = FALSE, strip_d = FALSE)
    ) %>%
    dplyr::mutate(
      time = "Y18",
      v_identity = as.numeric(v_identity),
      j_identity = as.numeric(j_identity),
      consensus_count = NA_real_,
      umi_count = NA_real_,
      v_call_10x = NA_character_,
      j_call_10x = NA_character_,
      junction_10x = NA_character_,
      junction_10x_aa = NA_character_,
      c_call = NA_character_
    )
}

clean_data <- function(df, source_name = NULL) {
  tryCatch({
    # Tag source if needed
    if (!"source" %in% names(df) && !is.null(source_name)) {
      df <- df %>% dplyr::mutate(source = source_name, .before = 1)
    }
    cleaned_df <- df %>%
      # Remove existing clone_id if present
      dplyr::select(-any_of("clone_id")) %>%
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

    source_label <- if (!is.null(source_name)) source_name else "data"
    cat("✓ Cleaned", source_label, ":", nrow(cleaned_df), "sequences\n")
    return(cleaned_df)
  }, error = function(e) {
    error_label <- if (!is.null(source_name)) source_name else "data"
    cat("✗ Error cleaning", error_label, ":", e$message, "\n")
    stop(e)
  })
}

# --------------------- Process datasets ---------------------
# Prepare C02 and G11 data
cat("Processing C02 and G11 data...\n")
c02_clean <- c02_data %>% clean_sorted_cultured() %>% clean_data()
g11_clean <- g11_data %>% clean_sorted_cultured() %>% clean_data()

# Prepare bulk data
cat("Processing bulk data...\n")
bulk_clean <- bulk_data %>%
  filter(type == "Bulk BCR ") %>%
  clean_data("bulk") %>%
  mutate(
    cell_id = sequence_id,
    clonal_family = NA_character_
  )

# Prepare 10X data
cat("Processing 10X heavy chain data...\n")
tenx_heavy_clean <- tenx_heavy %>%
  clean_data("10x_heavy") %>%
  mutate(
    time = "Y16",
    clonal_family = NA_character_
  )

cat("Processing 10X light chain data...\n")
tenx_light_clean <- tenx_light %>%
  clean_data("10x_light") %>%
  mutate(
    time = "Y16",
    clonal_family = NA_character_
  )

# --------------------- Combine datasets ---------------------
# Get all unique column names
all_datasets <- list(bulk_clean, c02_clean, g11_clean, tenx_heavy_clean, tenx_light_clean)

# Combine all datasets
cat("Combining datasets...\n")
combined_data <- dplyr::bind_rows(all_datasets) %>% distinct()

# Clean up the combined dataset
cat("Cleaning combined dataset...\n")
combined_clean <- combined_data %>%
  # Remove rows with no sequence_id
  filter(!is.na(sequence_id) & sequence_id != "") %>%
  # Standardize locus names
  mutate(locus = case_when(
    locus %in% c("IGH", "IgH") ~ "IGH",
    locus %in% c("IGK", "IgK") ~ "IGK",
    locus %in% c("IGL", "IgL") ~ "IGL",
    TRUE ~ locus
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
print(table(combined_clean$source))
cat("\nLocus distribution:\n")
print(table(combined_clean$locus, useNA = "ifany"))
cat("\nProductivity:\n")
print(table(combined_clean$productive, useNA = "ifany"))

# --------------------- Save files ---------------------
# Save cell IDs from C02 and G11 data for downstream checks
c02_g11_cell_ids <- unique(c(g11_clean$cell_id, c02_clean$cell_id))
write_csv(data.frame(cell_id = c02_g11_cell_ids), file.path(opt$outdir, "c02_g11_cell_ids.csv"))

# Save the combined dataset
cat("\nSaving combined dataset...\n")
write_csv(combined_clean, file.path(opt$outdir, "all_combined_data.csv"))

# Create heavy chain only dataset for clonal clustering
heavy_only <- combined_clean %>%
  filter(locus == "IGH" & productive == TRUE) %>%
  # Remove rows with missing critical fields for clonal clustering
  filter(!is.na(v_call) & !is.na(j_call) & !is.na(junction))

cat("\nHeavy chain productive sequences:", nrow(heavy_only), "\n")
cat("Heavy chain data sources:\n")
print(table(heavy_only$source))

# Save heavy chain dataset for clonal clustering
write_tsv(heavy_only, file.path(opt$outdir, "all_combined_heavy_chain_data.tsv"))

cat("\nFiles created:")
cat("\n- all_combined_data.csv (all data)")
cat("\n- all_combined_heavy_chain_data.tsv (IGH productive for clonal clustering)")
cat("\nData integration complete!\n")