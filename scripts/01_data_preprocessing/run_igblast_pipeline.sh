#!/usr/bin/env bash
#SBATCH --job-name=run_igblast_pipeline
#SBATCH --mem=8g
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --account=hoehnlab
#SBATCH --time=4:00:00
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority

export PATH=$PATH:/dartfs-hpc/rc/home/c/f0070d5/.local/bin/
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
module load python/anaconda3
conda activate r_phylo

# Define run parameters and input files
dir="/dartfs/rc/lab/H/HoehnK/Sherry/hiv_analysis"
indir="$dir/data/raw"
procdir="$dir/data/processed"
logdir="$dir/log"
pipeline_log="$dir/log/preprocessing.log"

# Define naming variables
contig_base="filtered_contig_Susan_Moir"
fasta_file="filtered_contig_Susan_Moir.fasta"
csv_file="filtered_contig_annotations_Susan_Moir.csv"

# Make output directory and empty log files
mkdir -p "$dir" "$procdir" "$logdir" \
         "$procdir/01_igblast_output" \
         "$procdir/02_parsed_sequences"
exec > >(tee -a "$pipeline_log") 2>&1

cd $dir

# Start
echo "OUTPUT DIRECTORY: ${procdir}"
echo -e "START"
STEP=0

# Process filtered_contig Susan Moir.fasta
echo "=== Processing ${fasta_file} ==="

# Assign V,D,J gene annotations
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssignGenes igblast"
AssignGenes.py igblast \
   -s "$indir/10x_data/${fasta_file}" \
   -b ~/share/igblast/ \
   --outdir $procdir/01_igblast_output --nproc 4 \
   --organism human --loci ig --format blast \
   --exec ~/programs/ncbi-igblast-1.22.0/bin/igblastn

# Create AIRR BCR-seq alignment database from IgBLAST output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb igblast"
MakeDb.py igblast \
   -i "$procdir/01_igblast_output/${contig_base}_igblast.fmt7" \
   -s "$indir/10x_data/${fasta_file}" \
   -r ~/share/germlines/imgt/human/vdj/ \
   --outdir $procdir/01_igblast_output \
   --10x "$indir/10x_data/${csv_file}" \
   --extended

# Split into separate light and heavy chain files
printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select heavy ${contig_base}"
   ParseDb.py select -d $procdir/01_igblast_output/${contig_base}_igblast_db-pass.tsv \
   -f locus -u "IGH" \
   --logic all --regex \
   --outname "${contig_base}_heavy" \
   --outdir $procdir/02_parsed_sequences

   printf " %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseDb select light ${contig_base}"
   ParseDb.py select -d $procdir/01_igblast_output/${contig_base}_igblast_db-pass.tsv \
   -f locus -u "IG[LK]" \
   --logic all --regex \
   --outname "${contig_base}_light" \
   --outdir $procdir/02_parsed_sequences

# End
echo "=== DONE $(date '+%F %T') ==="
cd ../
