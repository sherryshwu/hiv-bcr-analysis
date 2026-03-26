#!/usr/bin/env bash
#SBATCH --job-name=build_lineage_trees
#SBATCH --mem=16g
#SBATCH --time=72:00:00
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
Rscript "$SCRIPTS_DIR/02_clonal_analysis/extract_targeted_clones.R" \
  --indir "$SCOPER_DIR" \
  --outdir "$ANALYSIS_DIR"

# Expected outputs 
CSV_SPLIT="$ANALYSIS_DIR/targeted_clones_split_light.csv"

# Verify required files exist
echo "Checking for required input files..."
if [[ ! -f "$CSV_SPLIT" ]]; then
  echo "ERROR: Missing $CSV_SPLIT"; exit 1
fi

# Step 6a: Build phylogenetic trees (heavy only)
printf "\n%2d: %-40s %s\n" $((++STEP)) "Building phylogenetic trees (heavy)" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/build_trees.R" \
  --input "$CSV_SPLIT" \
  --outdir "$TREE_DIR" \
  --min_size 3 \
  --partition heavy_only \
  --nproc 8

# Step 6b: Build trees (heavy + light)
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

# Step 6c: Build phylogenetic trees (light only)
printf "\n%2d: %-40s %s\n" $((++STEP)) "Building phylogenetic trees (light)" "$(date +'%H:%M %D')"
Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/build_trees.R" \
  --input "$CSV_SPLIT" \
  --outdir "$TREE_DIR" \
  --min_size 3 \
  --partition light_only \
  --nproc 8

# Step 7: Run the combined correlation analysis
printf "\nRunning Final Correlation Analysis %s\n" "$(date +'%H:%M')"
Rscript "$SCRIPTS_DIR/03_phylogenetic_analysis/run_correlation_test.R"

echo "=== Downstream Analysis Complete at $(date) ==="
echo "Results saved to:"
echo "  - Clone analysis: $ANALYSIS_DIR"
echo "  - Trees: $TREE_DIR"