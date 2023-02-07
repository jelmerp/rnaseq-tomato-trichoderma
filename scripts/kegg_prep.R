# SET-UP -----------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "KEGGREST",
              "clusterProfiler")
pacman::p_load(char = packages)

# Load scripts with functions
source("scripts/GO/GO_fun.R")
source("scripts/kegg/kegg_fun.R")

# Define input files
DE_dir <- "results/DE"
KO_file <- "results/kegg/annot/ITAG4.1_KOlist_matthew.txt"

# Define output files
outdir <- "results/kegg/annot"
outdir_terms <- file.path(outdir, "indiv_terms")
kegg_map_file <- file.path(outdir, "kegg_map.txt")

dir.create(outdir_terms, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# MAKE THE KEGG MAP ------------------------------------------------------------
# Read input files
K_term_df_all <- read.table(KO_file,
                 sep = "\t", header = FALSE, fill = TRUE,
                 col.names = c("gene", "K_term")) %>%
  select(K_term, gene)

# Check and report stats
ngenes_noterm <- sum(grepl("^$", K_term_df_all$K_term))
ngenes_term <- sum(!grepl("^$", K_term_df_all$K_term))
cat("Number of genes _with_ a KEGG term:", ngenes_term, "\n") # 11413
cat("Number of genes _without_ a term:", ngenes_noterm, "\n") # 23275

# Filter
K_term_df <- K_term_df_all %>%
  filter(K_term != "") %>%                   # Remove genes with no associated terms
  mutate(gene = sub("\\.\\d$", "", gene))    # Proteins have extra suffix that genes don"t have

# Get KEGG pathways associated with the different K-terms
K_terms <- sort(unique(K_term_df$K_term))
pathway_df <- do.call(rbind, lapply(K_terms, get_pathway, outdir_terms))

# Create the KEGG map
kegg_map <- merge(K_term_df, pathway_df, by = "K_term") %>%
  select(pathway_id, pathway_description, gene, K_term) %>%
  arrange(pathway_id) # This will take an hour or so

# Write files
write_tsv(kegg_map, kegg_map_file)
