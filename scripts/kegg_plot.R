# SETUP ------------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse")       # Misc. data manipulation and plotting
pacman::p_load(char = packages)

# Source script with functions
source("mcic-scripts/rnaseq/rfuns/enrich_funs.R")

# Input files
kegg_res_file <- "results/kegg/kegg_results.txt"

# Output dir
outdir <- "results/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
figfile_up <- file.path(outdir, "SuppFig3_KEGG_up.png")
figfile_down <- file.path(outdir, "SuppFig4_KEGG_down.png")


# LOAD AND CHECK DATA ----------------------------------------------------------
# Read kegg results
kegg_res <- read_tsv(kegg_res_file, show_col_types = FALSE) |>
  dplyr::rename(category = ID,
                padj = p.adjust,
                description = pathway_description) |>
  mutate(
    numDEInCat = as.integer(sub("/\\d+", "", GeneRatio)),
    sig = ifelse(padj < 0.05, TRUE, FALSE),
    contrast = paste0(level1, "_vs_", level2),
    timepoint = str_extract(level1, "T0|T33"),
    irrigation = str_extract(level1, "drought|irrigated"),
    plant_part = str_extract(level1, "leaf|root")
  )

# Remove some irrelevant pathways
kegg_res <- kegg_res |>
  filter(!grepl("animal|quorum|human|anemia",
               description, ignore.case = TRUE))

# Check number of pathways per contrast
kegg_res %>% count(contrast)

# Settings
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")


# drought > irrigated ----------------------------------------------------------
p_up <- enrich_plot(kegg_res,
                    contrasts = contrast_sel,
                    DE_directions = c("up"),
                    n_tres = 2,
                    padj_tres = 0.01,
                    xlabs = c("T0", "T33"),
                    plot_ontologies = FALSE)
p_up

ggsave(figfile_up, p_up, bg = "white", width = 7, height = 5, dpi = "retina")


# drought < irrigated ----------------------------------------------------------
p_down <- enrich_plot(kegg_res,
                      contrasts = contrast_sel,
                      DE_directions = c("down"),
                      n_tres = 2,
                      padj_tres = 0.01,
                      xlabs = c("T0", "T33"),
                      plot_ontologies = FALSE)
p_down

ggsave(figfile_down, p_down, bg = "white", width = 7, height = 6, dpi = "retina")


# TO vs T33 --------------------------------------------------------------------
# No significant pathways
