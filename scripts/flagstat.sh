#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=flagstat
#SBATCH --output=slurm-bwa-flagstat-%j.out

## Setup
source ~/.bashrc
conda activate samtools-env

## Bash strict mode
set -euo pipefail

## Command-line args
bam="$1"
dir_flagstat="$2"

## Other variables
ncores="$SLURM_CPUS_PER_TASK"

# Report
echo "## Starting script flagstat.sh"
date
echo "## Input BAM file: $bam"
echo "## Output flagstat dir: $dir_flagstat"
echo -e "-------------------\n"

# Process args
sample=$(basename "$bam" .bam)
outfile="$dir_flagstat"/"$sample"_flagstat.txt

## Report
echo "## Sample ID: $sample"
echo "## Output file: $outfile"
echo

## Create output dir
mkdir -p "$dir_flagstat"

## Map with BWA
echo "Running samtools flagstats...."

samtools flagstat -@ "$ncores" "$bam" > "$outfile"

## Report
echo -e "\n\n----------------"

echo "Output file:"
ls -lh "$outfile"

echo "Done with script flagstat.sh"
date
