#!/usr/bin/env bash
#SBATCH --job-name=merge_datasets
#SBATCH --mem=32g
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=hoehnlab
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

set -euo pipefail

export PATH=$PATH:/dartfs-hpc/rc/home/c/f0070d5/.local/bin/
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate r_phylo

# Define directories
PROJECT_DIR="/dartfs/rc/lab/H/HoehnK/Sherry/hiv_analysis"
DATA_DIR="$PROJECT_DIR/data"
PROCESSED_DIR="$DATA_DIR/processed"
RESULTS_DIR="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/log"
SCRIPTS_DIR="$PROJECT_DIR/scripts"
GERMLINE_DB="$HOME/share/germlines/imgt/human/vdj"

# Create directories
mkdir -p "$RESULTS_DIR" "$LOG_DIR"

# Setup logging
PIPELINE_LOG="$LOG_DIR/02_merge_datasets.log"
exec > >(tee "$PIPELINE_LOG") 2>&1

cd "$PROJECT_DIR"

echo "=========================================="
echo "   Datasets Merging and Processing"
echo "=========================================="
echo "Started at: $(date)"
echo ""

STEP=0

# Check if all required files from IgBLAST were created
if [[ ! -f "$PROCESSED_DIR/01_igblast_output/sorted_cultured_db-pass.tsv" ]] || \
   [[ ! -f "$PROCESSED_DIR/02_parsed_sequences/filtered_contig_Susan_Moir_heavy_parse-select.tsv" ]]; then
    echo "ERROR: IgBLAST annotation failed!"
    exit 1
fi

# Merge datasets
printf "\n%2d: %-40s %s\n" $((++STEP)) "Merging BCR datasets" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/01_data_preprocessing/merge_bcr_datasets.R" \
  --outdir "$PROCESSED_DIR/03_combined_datasets" \
  --use_igblast TRUE

# Check if merge was successful
COMBINED_DATA="$PROCESSED_DIR/03_combined_datasets/all_combined_data.csv"
COMBINED_HEAVY_CHAIN_DATA="$PROCESSED_DIR/03_combined_datasets/all_combined_heavy_chain_data.tsv"
if [[ ! -f "$COMBINED_DATA" ]]; then
    echo "ERROR: Dataset merging failed!"
    exit 1
fi

echo ""
echo "=========================================="
echo "   Datasets Merging and Processing Complete at: $(date '+%F %T')"
echo "=========================================="
echo ""