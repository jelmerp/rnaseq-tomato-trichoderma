#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=STAR_index
#SBATCH --output=slurm-STAR-index-%j.out

# Setup
source ~/.bashrc
source activate star-env

set -e -u -o pipefail

echo "Starting script star_index.sh"
date 
echo -e "-------------------\n"

# Define, get and process variables:
ref_fa="$1"
index_dir="$2"
index_size="${3-13}"

mkdir -p "$index_dir"

## Report:
echo "Genome FASTA file: $ref_fa"
echo "Genome index dir (output): $index_dir"
echo "genomeSAindexNbases: $index_size"
echo

# Run STAR:
echo "Running STAR...."

# Run STAR:
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
    --runMode genomeGenerate \
    --genomeDir "$index_dir" \
    --genomeFastaFiles "$ref_fa" \
    --genomeSAindexNbases "$index_size" \
    --sjdbOverhang ReadLength-1

# Report:
echo -e "\n\n----------------"
echo "Done with script star_index.sh"
date
