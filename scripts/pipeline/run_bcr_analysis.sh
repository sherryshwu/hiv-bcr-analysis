#!/usr/bin/env bash
#SBATCH --job-name=bcr_full_pipeline
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --account=hoehnlab
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

set -euo pipefail

export PATH=$PATH:/dartfs-hpc/rc/home/c/f0070d5/.local/bin/
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate r_phylo

# Define directories
dir="/dartfs/rc/lab/H/HoehnK/Sherry/hiv_analysis"
outdir="$dir/results"
logdir="$dir/log"
procdir="$dir/data/processed"
pipeline_log="$logdir/full_pipeline_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$outdir" "$logdir"
exec > >(tee -a "$pipeline_log") 2>&1

cd "$dir"

echo "=== BCR Analysis Full Pipeline at $(date) ==="

STEP=0

# Step 1: Merge datasets
printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Merging BCR datasets"
Rscript scripts/01_data_preprocessing/merge_bcr_datasets.R --outdir $procdir

# Check if merge was successful
if [[ ! -f "$procdir/all_combined_bcr_data.csv" ]]; then
    echo "ERROR: Dataset merging failed!"
    exit 1
fi

# Step 2: Determine threshold and prepare Change-O input
printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Threshold determination"
Rscript scripts/02_clonal_analysis/determine_clustering_threshold.R \
    -i "$procdir/all_combined_heavy_chain_data.tsv" \
    --outdir "$outdir/clustering_threshold" \
    --nproc 8

# Step 3: Run Change-O DefineClones if threshold was found
if [[ -f "$outdir/clustering_threshold/threshold.csv" ]]; then
    THRESHOLD=$(awk -F',' 'NR==2 {print $2}' "$outdir/clustering_threshold/threshold.csv")
    echo "Using Change-O threshold: $THRESHOLD"
    
    printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "DefineClones"
    DefineClones.py \
    -d "$procdir/all_combined_heavy_chain_data.tsv" \
    --mode gene \
    --act set \
    --model ham \
    --norm len \
    --dist "$THRESHOLD" \
    --nproc 4 \
    --outdir "$outdir/clustering_changeo" \
    --log "$logdir/define_clones.log"
    
    mkdir -p "$outdir/clustering_changeo"

    printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "CreateGermlines"
    CreateGermlines.py -d "$outdir/clustering_changeo/all_combined_heavy_chain_data_clone-pass.tsv" \
    -r ~/share/germlines/imgt/human/vdj/ \
    --cloned \
    --outdir "$outdir/clustering_changeo"
else
    echo "Warning: No Change-O threshold found, skipping DefineClones"
fi

# Step 4: Run Scoper clonal analysis
printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 30 "Scoper clonal analysis"
Rscript scripts/02_clonal_analysis/run_scoper_clonal_clustering.R \
    -a "$procdir/all_combined_bcr_data.csv" \
    -t "$outdir/clustering_threshold/threshold.csv" \
    --outdir "$outdir/clustering_scoper" \
    --nproc 8

echo "=== Full Pipeline Complete at $(date) ==="
echo "Results in:"
echo "- $outdir/clustering_scoper/ (Threshold determination and Scoper results)"
echo "- $outdir/clustering_changeo/ (Change-O DefineClones results)"