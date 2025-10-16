#!/usr/bin/env bash
#SBATCH --job-name=bcr_full_pipeline
#SBATCH --mem=32g
#SBATCH --time=24:00:00
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
PIPELINE_LOG="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

cd "$PROJECT_DIR"

echo "=== BCR Analysis Full Pipeline at $(date) ==="

STEP=0

# Step 1: Merge datasets
printf "\n%2d: %-40s %s\n" $((++STEP)) "Merging BCR datasets" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/01_data_preprocessing/merge_bcr_datasets.R" \
  --outdir "$PROCESSED_DIR/03_combined_datasets"

# Check if merge was successful
COMBINED_DATA="$PROCESSED_DIR/03_combined_datasets/all_combined_data.csv"
COMBINED_HEAVY_CHAIN_DATA="$PROCESSED_DIR/03_combined_datasets/all_combined_heavy_chain_data.tsv"
if [[ ! -f "$COMBINED_DATA" ]]; then
    echo "ERROR: Dataset merging failed!"
    exit 1
fi

# Step 2: Determine threshold and prepare Change-O input
printf "\n%2d: %-40s %s\n" $((++STEP)) "Computing clustering threshold" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/02_clonal_analysis/determine_clustering_threshold.R" \
    --heavy_data "$COMBINED_HEAVY_CHAIN_DATA" \
    --outdir "$RESULTS_DIR/clustering_threshold" \
    --nproc 8

# Step 3: Run Change-O DefineClones if threshold was found
THRESHOLD_FILE=$RESULTS_DIR/clustering_threshold/threshold.csv

if [[ -f "$THRESHOLD_FILE" ]]; then
    THRESHOLD=$(awk -F',' 'NR==2 {print $2}' "$THRESHOLD_FILE")
    echo "Using Change-O threshold: $THRESHOLD"
    
    printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "DefineClones"

    CHANGEO_DIR="$RESULTS_DIR/clustering_changeo"
    mkdir -p "$CHANGEO_DIR"

    DefineClones.py \
    -d "$COMBINED_HEAVY_CHAIN_DATA" \
    --mode gene \
    --act set \
    --model ham \
    --norm len \
    --dist "$THRESHOLD" \
    --nproc 4 \
    --outdir "$CHANGEO_DIR" \
    --log "$LOG_DIR/define_clones.log"
    
    # Create germlines
    printf "\n%2d: %-40s %s\n" $((++STEP)) "Creating germlines" "$(date +'%H:%M %D')"
    CreateGermlines.py -d "$CHANGEO_DIR/all_combined_heavy_chain_data_clone-pass.tsv" \
    -r "$GERMLINE_DB" \
    --cloned \
    --outdir "$CHANGEO_DIR"
else
    echo "Warning: No Change-O threshold found, skipping DefineClones"
fi

# Step 4: Run Scoper clonal analysis
SCOPER_DIR="$RESULTS_DIR/clustering_scoper"

printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Scoper clonal analysis"
Rscript scripts/02_clonal_analysis/run_scoper_clonal_clustering.R \
    --all_data "$COMBINED_DATA" \
    --threshold_file "$THRESHOLD_FILE" \
    --outdir "$SCOPER_DIR" \
    --nproc 16 \
    --germline_dir "$GERMLINE_DB"

echo "=== Pipeline Part I Complete at $(date) ==="
echo "Results saved to:"
echo "  - $CHANGEO_DIR"
echo "  - $SCOPER_DIR"
echo "Proceeding to downstream analysis..."