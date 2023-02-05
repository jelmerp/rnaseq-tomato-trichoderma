# SET-UP -----------------------------------------------------------------------

## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "goseq")           # GO statistics
pacman::p_load(char = packages)

## Source script with functions
source("scripts/GO/GO_fun.R")

## Input files
DE_dir <- "results/DE"
gene_length_file <- "results/GO/map/ITAG4.1_gene-lengths.txt"
GO_map_file <- "results/GO/map/GO-map_Fattel.txt"

## Read input files
GO_map <- read.table(GO_map_file, header = TRUE) # (Tibble will cause problems for GOseq)
gene_length <- read_tsv(gene_length_file)

## Output dir
outdir <- "results/GO/pairwise_comps"
if(!dir.exists(outdir)) dir.create(outdir)


# RUN GENE ONTOLOGY ANALYSIS ---------------------------------------------------

## Define comparisons
comps <- list(
  c("leaf_drought_T33", "leaf_irrigated_T33"),
  c("leaf_drought_T0", "leaf_irrigated_T0"),
  c("leaf_irrigated_T0", "leaf_irrigated_T33"),
  c("leaf_drought_T0", "leaf_drought_T33"),
  c("leaf_drought", "leaf_irrigated")
)

## Run GO enrichment analysis for all comparisons
for(comp in comps) {
  #comp <- c("leaf_drought", "leaf_irrigated")

  cat("\n\n--------------------\nRegular DE - All:\n")
  DE_res_both <- get_DE_res(comp, DE_dir, LFC = FALSE)
  GO_res <- run_GO(comp, DE_res_both, GO_map, gene_length,
                   outdir, file_suffix = "_DE_both")

  cat("\n\n--------------------\nRegular DE - overexpressed in level 1:\n")
  DE_res_up <- DE_res_both %>% filter(log2FoldChange > 0)
  GO_res <- run_GO(comp, DE_res_up, GO_map, gene_length,
                   outdir, file_suffix = "_DE_up")

  cat("\n\n--------------------\nRegular DE - underexpressed in level 1:\n")
  DE_res_down <- DE_res_both %>% filter(log2FoldChange < 0)
  GO_res <- run_GO(comp, DE_res_down, GO_map, gene_length,
                   outdir, file_suffix = "_DE_down")
}

## For DE results with LFC-based thresholds:
# for(comp in comps) {
#   ## LFC threshold - all
#   cat("\n\n--------------------\nLFC-based DE - all:\n")
#   DE_res_both <- get_DE_res(comp, DE_dir, LFC = TRUE)
#   GO_res <- run_GO(comp, DE_res_both, GO_map, gene_length,
#                    outdir, file_suffix = "_LFC_both")
#
#   ## LFC threshold - overexpressed in level 1
#   cat("\n\n--------------------\nLFC-based DE - overexpressed in level 1:\n")
#   DE_res_up <- DE_res_both %>% filter(log2FoldChange > 0)
#   GO_res <- run_GO(comp, DE_res_up, GO_map, gene_length,
#                    outdir, file_suffix = "_LFC_up")
#
#   ## LFC threshold - underexpressed in level 1
#   cat("\n\n--------------------\nLFC-based DE - underexpressed in level 1:\n")
#   DE_res_down <- DE_res_both %>% filter(log2FoldChange < 0)
#   GO_res <- run_GO(comp, DE_res_down, GO_map, gene_length,
#                    outdir, file_suffix = "_LFC_down")
# }
