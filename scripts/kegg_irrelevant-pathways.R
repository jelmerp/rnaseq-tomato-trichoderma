# SET-UP -----------------------------------------------------------------------

## Install/load packages
if(! 'pacman' %in% installed.packages()) install.packages('pacman')
packages <- c('tidyverse',       # Misc. data manipulation and plotting
              'here',            # Managing file paths
              "clusterProfiler")
pacman::p_load(char = packages)

## Source script with functions
source("scripts/GO/GO_fun.R")
source("scripts/kegg/kegg_fun.R")

## Input files
DE_dir <- "results/DE"
kegg_annot_dir <- "results/kegg/annot/"
kegg_map_file <- file.path(kegg_annot_dir, "kegg_map.txt")

## Read input files
kegg_map_raw <- read_tsv(kegg_map_file)
length(unique(kegg_map_raw$gene)) # Nr of genes with a KEGG pathway: 7388

## Check human disease pathways
kegg_map_disease <- kegg_map_raw %>%
  filter(grepl("^ko05|^ko07", pathway_id)) # Include only human diseases

unique(kegg_map_disease$pathway_description)

## See which pathways are associated with the same K-term as "Nicotine addiction"
K_terms_nicotine <- kegg_map_disease %>%
  filter(pathway_description == "Nicotine addiction") %>%
  pull(K_term) %>%
  unique()
kegg_map_raw %>%
  filter(K_term %in% K_terms_nicotine) %>%
  pull(pathway_description) %>%
  unique()

## See which pathways are associated with the same K-term as "Legionellosis"
K_terms_legionellosis <- kegg_map_disease %>%
  filter(pathway_description == "Legionellosis") %>%
  pull(K_term) %>%
  unique()
kegg_map_raw %>%
  filter(K_term %in% K_terms_legionellosis) %>%
  pull(pathway_description) %>%
  unique()

## Process input files
kegg_map <- kegg_map_raw %>%
  distinct() %>%                           # Remove duplicate rows
  select(pathway_id, gene) %>%
  filter(!grepl("^ko05", pathway_id)) %>%  # Exclude human diseases
  filter(!grepl("^ko07", pathway_id))      # Exclude human diseases

pathway_descriptions <- kegg_map_raw %>%
  select(pathway_id, pathway_description) %>%
  filter(!grepl("^ko05", pathway_id)) %>%  # Exclude human diseases
  filter(!grepl("^ko07", pathway_id)) %>%       # Exclude human diseases
  distinct()
