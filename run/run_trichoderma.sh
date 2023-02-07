# MAP TO TRIASP rRNA -----------------------------------------------------------
# following https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/wednesday/contamination.html
ref_fa=metadata/ref/triasp/rRNA.fasta
dir_fastq=data/fastq_concat
dir_bam=results/bam/triasp_rrna
dir_flagstat=results/flagstat

for sample in S{01..40}; do
    echo $sample
    fq_R1=$dir_fastq/"$sample"_R1.fastq.gz
    bam="$dir_bam"/"$sample".bam
    sbatch scripts/align/bwa_align.sh "$fq_R1" "$bam" "$ref_fa"
done

for sample in S{01..40}; do
    echo $sample
    bam="$dir_bam"/"$sample".bam
    sbatch scripts/QC/flagstat.sh "$bam" "$dir_flagstat"
done

grep "mapped (" results/flagstat/*_flagstat.txt |
    sed -E 's@results/flagstat/(S[0-9][0-9])_flagstat.txt:([0-9]+) +.*@\1\t\2@' >results/flagstat/triasp_rrna_counts.txt

# MAP TO TRIASP (ASPERELLOIDES - THIS IS THE FOCAL SPECIES) --------------------
# Variables
ref_fa=metadata/ref/triasp/Triasp1_AssemblyScaffolds.fasta
ref_gff=metadata/ref/triasp/Triasp1_GeneCatalog_genes_20150522.gff
index_dir=metadata/ref/triasp/STAR_index
dir_bam=results/bam/triasp
dir_count=results/count/triasp

# Index triasp (asperelloides) genome
index_size=11 # Index size: log2(35390000)/2 - 1) = 11.5 (Triasp genome size is 35.39 Mb)
sbatch scripts/align/star_index.sh "$ref_fa" "$index_dir" "$index_size"

# Map to triasp genome
scripts/align/star_align_loop.sh "$sample_list" "$dir_fq" "$index_dir" "$dir_bam"

# featureCounts for triasp
t_option=CDS
g_option=proteinId
sbatch scripts/count/featureCounts.sh "$dir_bam" "$dir_count" "$ref_gff" "$t_option" "$g_option"
## Shorter file/sample names in summary table:
sed -E 's_[^\t]+/(S[0-9][0-9])[^\t]+_\1_g' "$dir_count"/counts.txt.summary >"$dir_count"/counts.txt.summary_ed

# MAP TO TRIAS (ASPERELLUM) ----------------------------------------------------
ref_fa=metadata/ref/trias/GCF_003025105.1_Trias_v._1.0_genomic.fna
ref_gff=metadata/ref/trias/GCF_003025105.1_Trias_v._1.0_genomic.gff
index_dir=metadata/ref/trias/STAR_index
dir_bam=results/bam/trias
dir_count=results/count/trias

# Index trias (asperellum) genome
index_size=11
sbatch scripts/align/star_index.sh "$ref_fa" "$index_dir" "$index_size"

# Map to trias genome
scripts/align/star_align_loop.sh "$sample_list" "$dir_fq" "$index_dir" "$dir_bam"

# featureCounts for trias genome
t_option=gene
g_option=Name
sbatch scripts/count/featureCounts.sh "$dir_bam" "$dir_count" "$ref_gff" "$t_option" "$g_option"
## Shorter file/sample names in summary table:
sed -E 's_[^\t]+/(S[0-9][0-9])[^\t]+_\1_g' "$dir_count"/counts.txt.summary >"$dir_count"/counts.txt.summary_ed
