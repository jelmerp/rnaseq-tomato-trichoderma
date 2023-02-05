#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=bwa_align
#SBATCH --output=slurm-bwa-align-%j.out

## Setup
source ~/.bashrc
conda activate bwa-env
conda activate --stack samtools-env

## Bash strict mode
set -euo pipefail

## Command-line args
fq_R1="$1"
bam="$2"
ref_fa="$3"

## Other variables
ncores="$SLURM_CPUS_PER_TASK"

# Report
echo "## Starting script bwa_align.sh"
date
echo "## R1 input FASTQ file: $fq_R1"
echo "## Reference genome: $ref_fa"
echo "## BAM output file: $bam"
echo -e "-------------------\n"

## Process args
fq_R2=${fq_R1/R1/R2}

mkdir -p "$(dirname "$bam")"

## Report
echo "## R2 input FASTQ file: $fq_R2"
echo

## Create ref genome index if necessary
if [ ! -f "$ref_fa".amb ]; then
    echo "Making BWA genome index..."
    bwa index "$ref_fa"
fi

## Map with BWA
echo "Mapping with BWA...."

bwa mem \
    -t "$ncores" \
    "$ref_fa" \
    "$fq_R1" "$fq_R2" |
    samtools view -@ "$ncores" -bS > "$bam"

## Report
echo -e "\n\n----------------"

echo "Printing scontrol job stats:"
scontrol show job "$SLURM_JOB_ID"

echo "Output BAM file:"
ls -lh "$bam"

echo "Done with script bwa_align.sh"
date
