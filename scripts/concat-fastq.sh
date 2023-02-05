#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --out=slurm-concat-fastq-%j.out

set -e -u -o pipefail

# This script will concatenate FASTQ files from the same sample and with the same read direction
# but from different lanes (L001 and L002).

# Save command-line arguments in named variables:
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments: input_dir and output_dir" && exit 1
input_dir="$1"
output_dir="$2"

# Globbing pattern for dir names to search:
DIR_PATTERN="C_BP_*"

# Create output dir, if necessary:
mkdir -p "$output_dir"

# Report before starting the pogram:
echo "## Starting FASTQ concatenation script..."
date
echo "## Input dir:: $input_dir"
echo "## Output dir: $output_dir"
echo -e "---------\n\n"

for sample_dir in "$input_dir"/$DIR_PATTERN; do
    echo "## Directory: $sample_dir"
    
    # Each dir should contain FASTQ files for a single sample.
    # First, extract the sample ID (e.g. "S9") from the 1st FASTQ file name.
    # Then, zero-pad single-digit numbers (e.g. "S9" -> "S09") to ensure proper ordering.
    sample_id=$(basename "$(find "$sample_dir" -name "*fastq.gz" | head -n 1)" | sed -E 's/.*_(S[0-9]+).*/\1/' | sed -E 's/S([0-9])$/S0\1/')
    
    # Count input files - exit if not 2 files per read direction:
    n_R1=$(find "$sample_dir" -name "*_R1_001.fastq.gz" | wc -l)
    n_R2=$(find "$sample_dir" -name "*_R2_001.fastq.gz" | wc -l)
    [[ "$n_R1" != 2 ]] && echo "## ERROR: Did not find two R1 input files" && exit 1
    [[ "$n_R2" != 2 ]] && echo "## ERROR: Did not find two R2 input files" && exit 1

    # List input files:
    echo "## R1 input files:"
    find "$sample_dir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh
    echo "## R2 input files:"
    find "$sample_dir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh

    # Concatenate FASTQ files:
    find "$sample_dir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$output_dir"/"$sample_id"_R1.fastq.gz
    find "$sample_dir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$output_dir"/"$sample_id"_R2.fastq.gz

    # List output files:
    echo "## Output files:"
    du -h "$output_dir"/"$sample_id"_R1.fastq.gz
    du -h "$output_dir"/"$sample_id"_R2.fastq.gz

    echo -e "---------\n"
done

# Make concatenated files read-only to protect against overwriting etc:
echo "## Making concatenated files read-only..."
chmod a-w "$output_dir"/*fastq.gz

# Report:
echo -e "\n## Done with FASTQ concatenation script."
date

# To run:
# input_dir=data/210430_Pearlly_GSL-PY-2114-transfer
# output_dir=data/fastq_concat
# sbatch scripts/concat-fastq.sh "$input_dir" "$output_dir"