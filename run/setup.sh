# Check size:
du -sh /fs/scratch/PAS0472/data-transfer/soledad_b/210430_Pearlly_GSL-PY-2114-transfer # 138 GB

# Check number of files using md5-sum doc:
wc -l /fs/scratch/PAS0472/data-transfer/soledad_b/210430_Pearlly_GSL-PY-2114-transfer/md5.txt # 160

# Download *Trichoderma asperelloides* genome files from https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Triasp1
# Genome:  Trichoderma asperelloides TR356 v1.0: Triasp1_AssemblyScaffolds.fasta.gz ("unmasked" file)
# Annotation:  Trichoderma asperelloides TR356 v1.0: Triasp1_GeneCatalog_genes_20150522.gff.gz ("Filtered models" => "Genes")

# Download *Trichoderma asperellum* from https://www.ncbi.nlm.nih.gov/genome/18221?genome_assembly_id=370101
# (Also available at https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Trias1)
cd data/ref/trias
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/025/105/GCF_003025105.1_Trias_v._1.0/GCF_003025105.1_Trias_v._1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/025/105/GCF_003025105.1_Trias_v._1.0/GCF_003025105.1_Trias_v._1.0_genomic.fna.gz
gunzip GCF_*

# ==============================================================================
#                            UPLOAD TO SRA
# ==============================================================================
# Move all FASTQ files to one dir
mkdir -p data/fastq_all
find data/210430_Pearlly_GSL-PY-2114-transfer -name "*fastq.gz" -exec mv -v {} data/fastq_all/ \;

# Create a dir with all files that need to be uploaded
mkdir -p data/fastq_ncbi
cp data/fastq_all/C_BP_[0-9]_* data/fastq_ncbi
cp data/fastq_all/C_BP_1[0-9]_* data/fastq_ncbi
cp data/fastq_all/C_BP_20_* data/fastq_ncbi

ftp -i ftp-private.ncbi.nlm.nih.gov # name: subftp / pw: ...
cd uploads/rawal.27_osu.edu_2nFEStRf
mkdir jelmer_upload
cd jelmer_upload
mput *
