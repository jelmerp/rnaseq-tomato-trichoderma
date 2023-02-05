# SETUP ------------------------------------------------------------------------
# Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",
              "ggforce")            # Managing file paths
pacman::p_load(char = packages)

# Source script with functions
source("scripts/GO_fun.R")

# Input files
indir <- "results/GO/ranjana/pairwise_comps" # USING RANJANA'S FILES
#indir <- "results/GO/pairwise_comps"         # USING MY FILES
GO_res_files <- list.files(indir, full.names = TRUE, pattern = "_GO.txt")
names(GO_res_files) <- sub("_GO.txt", "", basename(GO_res_files))

# Output dir
outdir <- "results/GO/fig/"
if(!dir.exists(outdir)) dir.create(outdir)
fig3 <- here(outdir, "Fig3.png")
fig5 <- here(outdir, "Fig5.png")
fig6 <- here(outdir, "Fig6.png")


# LOAD AND CHECK DATA ----------------------------------------------------------
## Read all GO results into one dataframe
GO_res <- map_dfr(GO_res_files, read_tsv,
                  col_types = "nnnnccccc",
                  .id = "contrast_full") %>%
  rename(padj = p_adj_over) %>%
  mutate(
    contrast = sub("(.*)_\\w+_\\w+$", "\\1", contrast_full),
    DE_test = sub(".*_(\\w+)_(\\w+)$", "\\1", contrast_full),
    direction = sub(".*_(\\w+)_(\\w+)$", "\\2", contrast_full),
    timepoint = sub(".*_(T0|T33).*", "\\1", level1),
    irrigation = sub(".*(drought|irrigated).*", "\\1", contrast),
    plant_part = sub(".*(leaf|root).*", "\\1", contrast),
    sig = ifelse(padj < 0.05, TRUE, FALSE)
    ) %>%
  filter(!is.na(term),
         ## I am focusing only on the regular DE test results here, not the LFC-based ones
         DE_test == "DE") %>%
  select(-p_adj_under, -level1, -level2)

## Check how many significant GO categories we have
GO_res %>% count(contrast, direction)


# FIG 3 - T0 v T33 -------------------------------------------------------------
f3 <- prep_goplot(GO_res,
                  contrasts = c("leaf_drought_T0_vs_leaf_drought_T33"),
                  directions = c("up", "down"),
                  pivot_by = "direction") %>%
  GO_plot(x_var = "direction", y_var = "cat_term")

ggsave(fig3, f3, bg = "white", width = 8, height = 6)


# FIG 5 - DROUGHT > IRRIGATED --------------------------------------------------
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")
f5 <- prep_goplot(GO_res,
                  contrast = contrast_sel,
                  directions = "up") %>%
  GO_plot(x_var = "contrast", y_var = "cat_term", xlabs = c("T0", "T33"))

ggsave(fig5, f5, bg = "white", width = 8, height = 6)

# FIG 6 - DROUGHT < IRRIGATED --------------------------------------------------
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")
f6 <- prep_goplot(GO_res,
                  contrast = contrast_sel,
                  directions = "down") %>%
  GO_plot(x_var = "contrast", y_var = "cat_term",
          xlabs = c("T0", "T33"), textsize = 8)

ggsave(fig6, f6, bg = "white", width = 8*1.2, height = 9*1.2)
