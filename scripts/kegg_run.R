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

## Output files
outdir <- "results/kegg"
if(!dir.exists(outdir)) dir.create(outdir)
kegg_results_file <- file.path(outdir, "kegg_results.txt")


# RUN KEGG ENRICHMENT ANALYSIS -------------------------------------------------

## Define pairwise contrasts
comps <- list(
  c("leaf_drought_T33", "leaf_irrigated_T33"),
  c("leaf_drought_T0", "leaf_irrigated_T0"),
  c("leaf_irrigated_T0", "leaf_irrigated_T33"),
  c("leaf_drought_T0", "leaf_drought_T33"),
  c("leaf_drought", "leaf_irrigated")
)
directions <- c("up", "down")
test_types <- c("DE")
### Create a dataframe with all combinations to test:
test_combs <- expand.grid(comp = comps, DE_direction = directions,
                          test_type = test_types,
                          stringsAsFactors = FALSE)

## Run KEGG analysis
kegg_list <- mapply(run_kegg,
                    test_combs$comp, test_combs$DE_direction, test_combs$test_type,
                    MoreArgs = list(DE_dir = DE_dir, kegg_map = kegg_map),
                    SIMPLIFY = FALSE)
kegg_df <- do.call(rbind, kegg_list) %>%
  merge(., pathway_descriptions, by.x = "ID", by.y = "pathway_id") %>%
  select(-Description)

## Check results
kegg_df %>% count(level1, level2, DE_direction)

## Write results to file
write_tsv(kegg_df, kegg_results_file)
