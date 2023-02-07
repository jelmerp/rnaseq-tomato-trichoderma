#!/bin/bash

# ==============================================================================
#                           FILES AND SETTINGS
# ==============================================================================
# Dirs and files
ref_fa=data/ref/tomato/Heinz1706_4.00/S_lycopersicum_chromosomes.4.00.fa
ref_gff=data/ref/tomato/annot/ITAG4.1/ITAG4.1_gene_models.gff
index_dir=data/ref/tomato/Heinz1706_4.00/STAR_index

dir_fq_raw=data/210430_Pearlly_GSL-PY-2114-transfer
dir_fq=data/fastq_concat
dir_fastqc=results/fastqc
dir_multiqc=results/multiqc
dir_bam=results/bam/tomato
dir_count=results/count/tomato
count_table="$dir_count"/counts.txt

# Featurecounts settings
t_option=gene
g_option=Name


# ==============================================================================
#                               PREP SEQUENCES
# ==============================================================================
# Run FastQC
for fq in "$dir_fq_raw"/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o "$dir_fastqc"
done

# Run MultiQC
sbatch scripts/qc/multiqc.sh -i "$dir_fastqc" -o "$dir_multiqc"

# Concatenate FASTQ files
sbatch scripts/concat-fastq.sh -i "$dir_fq_raw" -o "$dir_fq"


# ==============================================================================
#                               MAP AND COUNT
# ==============================================================================
# Index tomato genome
sbatch mcic-scripts/rnaseq/star_index.sh -i "$ref_fa" -o "$index_dir"

# Map to tomato genome
for R1 in "$dir_fq"/*_R1*fastq.gz; do
    scripts/rnaseq/star_align.sh -i "$R1" -r "$index_dir" -o "$dir_bam"
done

# featureCounts
sbatch mcic-scripts/rnaseq/featurecounts.sh \
    -i "$dir_bam" -a "$ref_gff" \
    -o "$count_table" \
    -t "$t_option" -g "$g_option"

# MultiQC
sbatch mcic-scripts/qc/multiqc.sh -i "$dir_count" -o "$dir_multiqc"/featurecounts
