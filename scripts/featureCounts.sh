#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=4
#SBATCH --out=slurm-featureCounts-%j.out

# Software
module load python/3.6-conda5.2
source ~/.bashrc
source activate subread-env

# Bash strict mode
set -e -u -o pipefail

# Command-line args
bam_dir=$1
output_dir=$2
gff_file=$3
t_option=${4-gene}           # "gene" for tomato
g_option=${5-Name}           # "Name" for tomato

# Process args
output_file=$output_dir/counts.txt

mkdir -p "$output_dir"

# Report
echo "## Starting script featureCounts.sh"
date 
echo "## BAM dir: $bam_dir"
echo "## Output file: $output_file"
echo "## GFF file: $gff_file"
echo "## -t option (feature type): $t_option"
echo "## -g option (aggregration ID): $g_option"
echo
echo "## Number of BAM files: $(find "$bam_dir"/*bam | wc -l)"
echo -e "-------------------\n"

# Run featureCounts
featureCounts -s 2 -p -B -C \
    -T "$SLURM_CPUS_ON_NODE" \
    -t "$t_option" -g "$g_option" \
    -a "$gff_file" \
    -o "$output_file" \
    "$bam_dir"/*bam

# Report
echo -e "\n\n----------------"

echo "## Printing scontrol job stats:"
scontrol show job "$SLURM_JOB_ID"

echo -e "\n## Done with script featureCounts.sh"
date