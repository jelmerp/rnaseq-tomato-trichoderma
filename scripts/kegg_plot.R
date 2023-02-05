# SETUP ------------------------------------------------------------------------
# Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here")            # Managing file paths
pacman::p_load(char = packages)

# Source script with functions
source("scripts/kegg/kegg_fun.R")
source("scripts/GO/GO_fun.R")

# Input files
kegg_res_file <- "results/kegg/kegg_results.txt"

# Output dir
outdir <- "results/kegg/fig/"
if(!dir.exists(outdir)) dir.create(outdir)


# LOAD AND CHECK DATA ----------------------------------------------------------

## Read kegg results
kegg_res <- read_tsv(kegg_res_file, show_col_types = FALSE) %>%
  mutate(
    contrast = paste0(level1, "_vs_", level2),
    timepoint = str_extract(level1, "T0|T33"),
    irrigation = str_extract(level1, "drought|irrigated"),
    plant_part = str_extract(level1, "leaf|root")
  )

kegg_res %>% count(contrast)

# CREATE PLOTS ----------------------------------------------------------------

# contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
#                   "leaf_drought_T33_vs_leaf_irrigated_T33",
#                   "leaf_drought_vs_leaf_irrigated")
# xlabs <- c("2-factor", "T0", "T33")
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")
xlabs <- c("T0", "NT33")

## Create plot: drought > irrigated at either timepoint
title <- "water deficit > irrigated"

kegg_res_focal <- kegg_res %>%
  filter(contrast %in% contrast_sel,
         DE_direction == "up",             # Overexpressed in drought
         p.adjust < 0.01)

p <- kegg_plot(kegg_res_focal, "timepoint",
          title = title, xlabs = xlabs)
fig_file <- here(outdir, paste0(title, ".png"))
ggsave(fig_file, p, bg = "white", width = 7, height = 6)


## Create plot: drought < irrigated at either timepoint
title <- "water deficit < irrigated"

kegg_res_focal <- kegg_res %>%
  filter(contrast %in% contrast_sel,
         DE_direction == "down",           # Underexpressed in drought
         p.adjust < 0.01)

p <- kegg_plot(kegg_res_focal, "timepoint",
               title = title, xlabs = xlabs)
fig_file <- here(outdir, paste0(title, ".png"))
ggsave(fig_file, p, bg = "white", width = 7, height = 6)


## No significant pathways for T0 v T33

