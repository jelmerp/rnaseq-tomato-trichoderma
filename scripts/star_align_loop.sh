#!/bin/bash

set -e -u -o pipefail

echo "Starting script star_align_loop.sh"
date 
echo -e "-------------------\n"

# Define, get and process variables:
sample_list="$1"
fastq_dir="$2"
index_dir="$3"
bam_dir="$4"

mkdir -p "$bam_dir"

# Report:
echo "FASTQ dir: $fastq_dir"
echo "Sample list: $sample_list"
echo "Genome index dir: $index_dir"
echo "Bam (output) directory: $bam_dir"
echo

# Loop:
samples=($(cat "$sample_list"))

for sample in "${samples[@]}"; do

  echo "Submitting star_align.sh for $sample ..."
  slurm_log=slurm-align-"$sample"-%j.out 
  sbatch -o "$slurm_log" scripts/align/star_align.sh "$sample" "$fastq_dir" "$index_dir" "$bam_dir"

done

# Report:
echo -e "\n\n----------------"
echo "Done with script star_index.sh"
date
