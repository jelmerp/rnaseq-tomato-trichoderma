# SET-UP -----------------------------------------------------------------------

# Install/load packages
if(! 'pacman' %in% installed.packages()) install.packages('pacman')
packages <- c('tidyverse', 'ape')
pacman::p_load(char = packages)

# Source script with functions
source("scripts/GO/GO_fun.R")

# Input files
GFF_file <- "data/ref/annot/ITAG4.1/ITAG4.1_gene_models.gff"
GO_file_Fattel <- "data/ref/annot/ITAG4.1/1.1_GOMAP-output.gaf" # https://www.biorxiv.org/content/10.1101/2021.04.25.441366v2

# Output files
outdir <- "results/GO/map"
if(!dir.exists(outdir)) dir.create(outdir)

GO_map_Fattel_file <- file.path(outdir, "GO-map_Fattel.txt")
gene_length_file <- file.path(outdir, "ITAG4.1_gene-lengths.txt")


# CREATE GO MAP ----------------------------------------------------------------

## Prep GO map -  using GO file Fattel et al.
GO_map <- read_tsv(GO_file_Fattel, skip = 1) %>%
  select(gene_id = db_object_id, GO_term = term_accession) #%>%
  filter(!is.na(GO_term))

## Check number of genes with GO annotation
cat("Number of genes in the GO map:\n")
nrow(GO_map)
#> 367,576
cat("Number of genes in the GO map:\n")
length(unique(GO_map$gene_id))
#> 33,493

write_tsv(GO_map, GO_map_Fattel_file)


# CREATE GENE LENGTH DATAFRAME -------------------------------------------------

## Prep gene lengths
gene_length <- read.gff(GFF_file) %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub(".*Name=(Solyc\\w+\\.\\d).*", "\\1", attributes),
         length = end - start + 1) %>%
  select(gene_id, length) %>%
  arrange(gene_id)

cat("GO map entries matching gene length dataframe (GFF):")
table(GO_map$gene_id %in% gene_length$gene_id)
#> FALSE   TRUE
#> 108     367468
# => All but 108 entries in the GO df are in the GFF

cat("Genes from gene length dataframe (GFF) matching GO map entries:")
table(gene_length$gene_id %in% GO_map$gene_id)
#> FALSE   TRUE
#> 1205    33483
# => All but 1205 genes in the GFF have a GO term

write_tsv(gene_length, gene_length_file)


# OLD --------------------------------------------------------------------------

## Prep GO map - Solnet ITAG4.0
# GO_file_solnet <- "data/ref/annot/ITAG4.1/ITAG4.0_goterms.txt"

# GO_solnet <- read.table(solnet_GO_file,
#                         fill = TRUE,
#                         col.names = c("gene_id", "GO_term"))
# nmax <- max(str_count(GO_solnet$GO_term, ",")) + 1
#
# GO_map <- GO_solnet %>%
#   separate(GO_term, into = paste("term", 1:nmax), sep = ",") %>%
#   pivot_longer(-gene_id, names_to = "term_count", values_to = "GO_term") %>%
#   mutate(GO_term = sub("^$", NA, GO_term)) %>%
#   drop_na() %>%
#   select(gene_id, GO_term) %>%
#   as.data.frame()
#
# GO_map <- GO_map %>%
#   mutate(gene_id = sub("\\.\\d$", "", gene_id)) %>%
#   as.data.frame()
