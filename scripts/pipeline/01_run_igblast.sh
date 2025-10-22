#!/usr/bin/env bash
#SBATCH --job-name=run_igblast
#SBATCH --mem=8g
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --account=hoehnlab
#SBATCH --time=2:00:00
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

export PATH=$PATH:/dartfs-hpc/rc/home/c/f0070d5/.local/bin/
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
module load python/anaconda3
conda activate r_phylo

# Define run parameters and input files
PROJECT_DIR="/dartfs/rc/lab/H/HoehnK/Sherry/hiv_analysis"
RAW_DIR="$PROJECT_DIR/data/raw"
PROC_DIR="$PROJECT_DIR/data/processed"
LOG_DIR="$PROJECT_DIR/log"
PREP_DIR="$PROC_DIR/00_sorted_cultured_prep"
IGBLAST_DIR="$PROC_DIR/01_igblast_output"
PARSED_DIR="$PROC_DIR/02_parsed_sequences"

# IgBLAST parameters
IGBLAST_DB="$HOME/share/igblast/"
GERMLINE_DB="$HOME/share/germlines/imgt/human/vdj/"
IGBLAST_EXEC="$HOME/programs/ncbi-igblast-1.22.0/bin/igblastn"

# Create directories
mkdir -p "$PROC_DIR" "$LOG_DIR" "$PREP_DIR" "$IGBLAST_DIR" "$PARSED_DIR"

# Setup logging
PIPELINE_LOG="$LOG_DIR/01_run_igblast.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

cd "$PROJECT_DIR"

echo "=========================================="
echo "   IgBLAST Annotation for All BCR Data"
echo "=========================================="
echo "Started at: $(date)"
echo ""

# Start
echo "OUTPUT DIRECTORY: ${procdir}"
echo -e "START"
STEP=0

# Part 1: Sorted/Cultured Sequences (C02 & G11)
echo "=== PART 1: Sorted/Cultured Sequences ==="

# Prepare sorted/cultured sequences
printf "\n%2d: %-40s %s\n" $((++STEP)) "Preparing C02/G11 sequences" "$(date +'%H:%M %D')"
Rscript scripts/01_data_preprocessing/prepare_sorted_cultured_for_igblast.R
SORTED_FASTA="$PREP_DIR/sorted_cultured_sequences.fasta"
if [[ ! -f "$SORTED_FASTA" ]]; then
    echo "ERROR: Sorted/cultured FASTA file not created!"
    exit 1
fi

# Run IgBLAST on sorted/cultured data
printf "\n%2d: %-40s %s\n" $((++STEP)) "IgBLAST (sorted/cultured)" "$(date +'%H:%M %D')"
AssignGenes.py igblast \
    -s "$SORTED_FASTA" \
    -b "$IGBLAST_DB" \
    --outdir "$IGBLAST_DIR" \
    --outname sorted_cultured \
    --nproc 4 \
    --organism human \
    --loci ig \
    --format blast \
    --exec "$IGBLAST_EXEC"

# Create database for sorted/cultured
printf "\n%2d: %-40s %s\n" $((++STEP)) "MakeDb (sorted/cultured)" "$(date +'%H:%M %D')"
MakeDb.py igblast \
    -i "$IGBLAST_DIR/sorted_cultured_igblast.fmt7" \
    -s "$SORTED_FASTA" \
    -r "$GERMLINE_DB" \
    --outdir "$IGBLAST_DIR" \
    --outname sorted_cultured \
    --extended

# Check output
SORTED_DB="$IGBLAST_DIR/sorted_cultured_db-pass.tsv"
if [[ ! -f "$SORTED_DB" ]]; then
    echo "ERROR: Sorted/cultured database not created!"
    exit 1
fi

echo "✓ Sorted/cultured IgBLAST annotation complete"
echo "  Output: $SORTED_DB"

# Part 2: 10X Single-Cell Data
echo ""
echo "=== PART 2: 10X Single-Cell Data ==="

# Define 10X file names
CONTIG_BASE="filtered_contig_Susan_Moir"
FASTA_FILE="$RAW_DIR/10x_data/${CONTIG_BASE}.fasta"
CSV_FILE="$RAW_DIR/10x_data/filtered_contig_annotations_Susan_Moir.csv"

# Check if 10X files exist
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "ERROR: 10X FASTA file not found: $FASTA_FILE"
    exit 1
fi
if [[ ! -f "$CSV_FILE" ]]; then
    echo "ERROR: 10X CSV file not found: $CSV_FILE"
    exit 1
fi

# Run IgBLAST on 10X data
printf "\n%2d: %-40s %s\n" $((++STEP)) "IgBLAST (10X data)" "$(date +'%H:%M %D')"
AssignGenes.py igblast \
    -s "$FASTA_FILE" \
    -b "$IGBLAST_DB" \
    --outdir "$IGBLAST_DIR" \
    --outname "$CONTIG_BASE" \
    --nproc 4 \
    --organism human \
    --loci ig \
    --format blast \
    --exec "$IGBLAST_EXEC"

# Create database for 10X data
printf "\n%2d: %-40s %s\n" $((++STEP)) "MakeDb (10X data)" "$(date +'%H:%M %D')"
MakeDb.py igblast \
    -i "$IGBLAST_DIR/${CONTIG_BASE}_igblast.fmt7" \
    -s "$FASTA_FILE" \
    -r "$GERMLINE_DB" \
    --outdir "$IGBLAST_DIR" \
    --outname "$CONTIG_BASE" \
    --10x "$CSV_FILE" \
    --extended

## Check output
TENX_DB="$IGBLAST_DIR/${CONTIG_BASE}_db-pass.tsv"
if [[ ! -f "$TENX_DB" ]]; then
    echo "ERROR: 10X database not created!"
    exit 1
fi

# Split 10X into heavy and light chains
printf "\n%2d: %-40s %s\n" $((++STEP)) "ParseDb (10X heavy)" "$(date +'%H:%M %D')"
ParseDb.py select \
    -d "$TENX_DB" \
    -f locus \
    -u "IGH" \
    --logic all \
    --regex \
    --outname "${CONTIG_BASE}_heavy" \
    --outdir "$PARSED_DIR"

printf "\n%2d: %-40s %s\n" $((++STEP)) "ParseDb (10X light)" "$(date +'%H:%M %D')"
ParseDb.py select \
    -d "$TENX_DB" \
    -f locus \
    -u "IG[LK]" \
    --logic all \
    --regex \
    --outname "${CONTIG_BASE}_light" \
    --outdir "$PARSED_DIR"

echo "✓ 10X IgBLAST annotation complete"
echo "  Heavy: $PARSED_DIR/${CONTIG_BASE}_heavy_parse-select.tsv"
echo "  Light: $PARSED_DIR/${CONTIG_BASE}_light_parse-select.tsv"

# Check if all required files were created
if [[ ! -f "$IGBLAST_DIR/sorted_cultured_db-pass.tsv" ]] || \
   [[ ! -f "$IGBLAST_DIR/filtered_contig_Susan_Moir_db-pass.tsv" ]]; then
    echo "ERROR: IgBLAST annotation failed!"
    exit 1
fi

# End
echo ""
echo "=========================================="
echo "   IgBLAST Annotation Complete at: $(date '+%F %T')"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Sorted/Cultured: $SORTED_DB"
echo "  10X (all):       $TENX_DB"
echo "  10X (heavy):     $PARSED_DIR/${CONTIG_BASE}_heavy_parse-select.tsv"
echo "  10X (light):     $PARSED_DIR/${CONTIG_BASE}_light_parse-select.tsv"
echo ""
cd ../