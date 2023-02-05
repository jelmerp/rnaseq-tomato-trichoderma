#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=STAR_align
#SBATCH --output=slurm-STAR-align-%j.out

# Setup
source ~/.bashrc
source activate star-env

# Bash strict mode
set -e -u -o pipefail

# Command-line args
sample="$1"
fastq_dir="$2"
index_dir="$3"
bam_dir="$4"

# Report
echo "## Starting script star_align.sh"
date
echo "## Sample ID: $sample"
echo "## FASTQ (input dir): $fastq_dir"
echo "## Genome index dir: $index_dir"
echo "## BAM (output) dir: $bam_dir"
echo -e "-------------------\n"

# Process args
shopt -s globstar # Recursive globbing
R1=$(ls "$fastq_dir"/**/*"$sample"*R1*fastq.gz)
R2=$(ls "$fastq_dir"/**/*"$sample"*R2*fastq.gz)

# Report
echo "## R1 input file: $R1"
echo "## R2 input file: $R2"
echo

# Create output dir
mkdir -p "$bam_dir"

# Run STAR
echo "Running STAR...."

STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
   --genomeDir "$index_dir" \
   --readFilesIn "$R1" "$R2" \
   --readFilesCommand zcat \
   --outFileNamePrefix "$bam_dir/$sample"_ \
   --outSAMtype BAM SortedByCoordinate \
   --alignIntronMin 5 --alignIntronMax 350000

# Move logfiles
mkdir -p "$bam_dir"/star_logs  
mv "$bam_dir"/"$sample"*out "$bam_dir"/"$sample"*tab "$bam_dir"/star_logs/

# Report:
echo -e "\n\n----------------"

echo "Printing scontrol job stats:"
scontrol show job "$SLURM_JOB_ID"

echo "Done with script star_align.sh"
date
