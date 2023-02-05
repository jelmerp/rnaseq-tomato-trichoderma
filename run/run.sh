
## Dirs and files
dir_fastq_raw=data/210430_Pearlly_GSL-PY-2114-transfer
dir_fastq_concat=data/fastq_concat
dir_fastqc=results/fastqc
dir_multiqc=results/multiqc
dir_bam=results/bam/tomato
dir_count=results/count/tomato

sample_list=metadata/samples.txt

ref_fa=data/ref/tomato/Heinz1706_4.00/S_lycopersicum_chromosomes.4.00.fa
ref_gff=data/ref/tomato/annot/ITAG4.1/ITAG4.1_gene_models.gff
index_dir=data/ref/tomato/Heinz1706_4.00/STAR_index

# PREP -------------------------------------------------------------------------
## Run FastQC
scripts/qc/fastqc.sh "$dir_fastq_raw" "$dir_fastqc"

## Run MultiQC
sbatch scripts/qc/multiqc.sh "$dir_fastq_raw" "$dir_multiqc"

## Concatenate FASTQ files
sbatch scripts/concat-fastq.sh "$dir_fastq_raw" "$dir_fastq_concat"

## Create list of samples
find data/fastq_concat -name "*fastq.gz" | sed 's/_R[12].fastq.gz//' | uniq > "$sample_list"


# MAP TO TOMATO ----------------------------------------------------------------
## Index tomato genome
sbatch scripts/align/star_index.sh "$ref_fa" "$index_dir"

## Map to tomato genome
scripts/align/star_align_loop.sh "$sample_list" "$dir_fastq_concat" "$index_dir" "$dir_bam"

## featureCounts
t_option=gene
g_option=Name
sbatch scripts/count/featureCounts.sh "$dir_bam" "$dir_count" "$ref_gff" "$t_option" "$g_option"
sbatch mcic-scripts/qc/multiqc.sh -i "$dir_count" -o "$dir_multiqc"/featurecounts


# KEGG MAP FOR TOMATO ----------------------------------------------------------
## Submit 7.500-record FASTAs to https://www.kegg.jp/blastkoala/ - subset FASTA:

for section in $(seq 0 "$last_section"); do
    start=$((1 + $section * 7500))
    end=$((7500 + $section * 7500))
    echo "Section $section - extraction from $start to $end"
    grep "^>" "$fa_in" | sed -n "${start},${end}p" | sed 's/>//' > "$outdir"/ITAG4.1_proteinIDs_pt"$section".txt
    seqtk subseq "$fa_in" "$outdir"/ITAG4.1_proteinIDs_pt"$section".txt > "$outdir"/ITAG4.1_proteins_pt"$section".fa
done

fa_in=metadata/ref/tomato/annot/ITAG4.1/ITAG4.1_proteins.fasta
outdir=results/kegg/fa_subset && mkdir -p "$outdir"
fa_lin="$outdir"/ITAG4.1_proteins_linearized.fa

n_total=$(grep -c "^>" "$fa_in")
last_section=$(python -c "import math; print(math.ceil($n_total / 7500) - 1)")

awk 'NR==1 { print $0; next } /^>/ { printf("\n%s\n",$0); next } { printf("%s",$0) } END {printf("\n") }' $fa_in > "$fa_lin"

for section in $(seq 0 "$last_section"); do
    outfile="$outdir"/ITAG4.1_proteins_pt"$section".fa
    start=$(( 1 + ($section * 7500 * 2) ))
    end=$(( 15000 + ($section * 7500 * 2 ) ))
    sed -n "${start},${end}p" "$fa_lin" > "$outfile"
    nrecs=$(grep -c "^>" "$outfile")
    echo "Section $section - extraction from $start to $end - Nr records: $nrecs"
done


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
    sed -E 's@results/flagstat/(S[0-9][0-9])_flagstat.txt:([0-9]+) +.*@\1\t\2@' > results/flagstat/triasp_rrna_counts.txt


# MAP TO TRIASP (ASPERELLOIDES - THIS IS THE FOCAL SPECIES) --------------------
## Variables
ref_fa=metadata/ref/triasp/Triasp1_AssemblyScaffolds.fasta
ref_gff=metadata/ref/triasp/Triasp1_GeneCatalog_genes_20150522.gff
index_dir=metadata/ref/triasp/STAR_index
dir_bam=results/bam/triasp
dir_count=results/count/triasp

## Index triasp (asperelloides) genome
index_size=11    # Index size: log2(35390000)/2 - 1) = 11.5 (Triasp genome size is 35.39 Mb)
sbatch scripts/align/star_index.sh "$ref_fa" "$index_dir" "$index_size"

## Map to triasp genome
scripts/align/star_align_loop.sh "$sample_list" "$dir_fastq_concat" "$index_dir" "$dir_bam"

## featureCounts for triasp
t_option=CDS
g_option=proteinId
sbatch scripts/count/featureCounts.sh "$dir_bam" "$dir_count" "$ref_gff" "$t_option" "$g_option"
### Shorter file/sample names in summary table:
sed -E 's_[^\t]+/(S[0-9][0-9])[^\t]+_\1_g' "$dir_count"/counts.txt.summary > "$dir_count"/counts.txt.summary_ed


# MAP TO TRIAS (ASPERELLUM) ----------------------------------------------------
ref_fa=metadata/ref/trias/GCF_003025105.1_Trias_v._1.0_genomic.fna
ref_gff=metadata/ref/trias/GCF_003025105.1_Trias_v._1.0_genomic.gff
index_dir=metadata/ref/trias/STAR_index
dir_bam=results/bam/trias
dir_count=results/count/trias

## Index trias (asperellum) genome
index_size=11
sbatch scripts/align/star_index.sh "$ref_fa" "$index_dir" "$index_size"

## Map to trias genome
scripts/align/star_align_loop.sh "$sample_list" "$dir_fastq_concat" "$index_dir" "$dir_bam"

## featureCounts for trias genome
t_option=gene
g_option=Name
sbatch scripts/count/featureCounts.sh "$dir_bam" "$dir_count" "$ref_gff" "$t_option" "$g_option"
### Shorter file/sample names in summary table:
sed -E 's_[^\t]+/(S[0-9][0-9])[^\t]+_\1_g' "$dir_count"/counts.txt.summary > "$dir_count"/counts.txt.summary_ed


# 2022-12-03 -- Upload to SRA --------------------------------------------------

## Move all FASTQ files to one dir
mkdir -p data/fastq_all
find data/210430_Pearlly_GSL-PY-2114-transfer -name "*fastq.gz" -exec mv -v {} data/fastq_all/ \;

## Create a dir with all files that need to be uploaded
mkdir -p data/fastq_ncbi
cp data/fastq_all/C_BP_[0-9]_* data/fastq_ncbi
cp data/fastq_all/C_BP_1[0-9]_* data/fastq_ncbi
cp data/fastq_all/C_BP_20_* data/fastq_ncbi

ftp -i ftp-private.ncbi.nlm.nih.gov # name: subftp / pw: 1CVwWYg7
cd uploads/rawal.27_osu.edu_2nFEStRf
mkdir jelmer_upload
cd jelmer_upload
mput *
