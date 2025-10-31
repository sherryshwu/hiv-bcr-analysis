#!/usr/bin/env bash
#SBATCH --job-name=build_lineage_trees
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
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
RESULTS_DIR="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/log"
SCRIPTS_DIR="$PROJECT_DIR/scripts"

# Setup logging
PIPELINE_LOG="$LOG_DIR/04_build_trees.log"
exec > >(tee "$PIPELINE_LOG") 2>&1

cd "$PROJECT_DIR"

echo "=== BCR Downstream Analysis Pipeline at $(date) ==="

STEP=0

# Use existing Scoper results
SCOPER_DIR="$RESULTS_DIR/clustering_scoper"
ANALYSIS_DIR="$RESULTS_DIR/clone_analysis"
TREE_DIR="$RESULTS_DIR/phylogenetic_trees"
mkdir -p "$ANALYSIS_DIR" "$TREE_DIR"

# Step 5: Pick target clone IDs for tree building
printf "\n%2d: %-40s %s\n" $((++STEP)) "Analyzing clone assignments" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/02_clonal_analysis/compare_clone_splitting_analysis.R" \
  --indir "$SCOPER_DIR" \
  --outdir "$ANALYSIS_DIR"

# Expected outputs 
CSV_NO_SPLIT="$ANALYSIS_DIR/hl_for_kept_clones_no_split_light.csv"
CSV_SPLIT="$ANALYSIS_DIR/hl_for_kept_clones_split_light.csv"

# Verify required files exist
echo "Checking for required input files..."
if [[ ! -f "$CSV_NO_SPLIT" ]]; then
  echo "ERROR: Missing $CSV_NO_SPLIT"; exit 1
fi
if [[ ! -f "$CSV_SPLIT" ]]; then
  echo "WARNING: Missing $CSV_SPLIT (will only build from no_split)"; 
fi

# Step 6a: Build phylogenetic trees (heavy only)
printf "\n%2d: %-40s %s\n" $((++STEP)) "Building phylogenetic trees (heavy)" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/build_trees.R" \
  --input "$CSV_NO_SPLIT" \
  --outdir "$TREE_DIR" \
  --min_size 3 \
  --partition heavy_only \
  --nproc 8

# Step 6b: Build phylogenetic trees (light only)
printf "\n%2d: %-40s %s\n" $((++STEP)) "Building phylogenetic trees (light)" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/build_trees.R" \
  --input "$CSV_SPLIT" \
  --outdir "$TREE_DIR" \
  --min_size 3 \
  --partition light_only \
  --nproc 8

# Step 6c: Build trees (heavy + light)
if [[ -f "$CSV_SPLIT" ]]; then
  printf "\n%2d: %-40s %s\n" $((++STEP)) "Building trees (heavy+light)" "$(date +'%H:%M %D')"
  Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/build_trees.R" \
    --input "$CSV_SPLIT" \
    --outdir "$TREE_DIR" \
    --min_size 3 \
    --partition both \
    --nproc 8
else
  echo "WARNING: Skipping heavy+light trees (files not available)"
fi

echo "=== Downstream Analysis Complete at $(date) ==="
echo "Results saved to:"
echo "  - Clone analysis: $ANALYSIS_DIR"
echo "  - Trees: $TREE_DIR"