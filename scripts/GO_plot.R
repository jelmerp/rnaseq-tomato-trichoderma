# SETUP ------------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "ggforce")         # Managing file paths
pacman::p_load(char = packages)

# Source script with functions
source("mcic-scripts/rnaseq/rfuns/enrich_funs.R")

# Define input files
indir <- "results/GO/ranjana/pairwise_comps"
GO_res_files <- list.files(indir, full.names = TRUE, pattern = "_GO.txt")
names(GO_res_files) <- sub("_GO.txt", "", basename(GO_res_files))

# Define output files
outdir <- "results/figures/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
figfile_wd <- file.path(outdir, "Fig3_GO_WD.png")
figfile_wd2 <- file.path(outdir, "Fig3_alt_GO_WD_min2inCat.png")
figfile_up <- file.path(outdir, "SuppFig1_GO_up.png")
figfile_down <- file.path(outdir, "SuppFig2_GO_down.png")


# LOAD AND CHECK DATA ----------------------------------------------------------
# Read all GO results into one dataframe
GO_res <- map_dfr(GO_res_files, read_tsv,
                  col_types = "nnnnccccc",
                  .id = "contrast_full") |>
  rename(padj = p_adj_over,
         description = term) |>
  mutate(
    contrast = sub("(.*)_\\w+_\\w+$", "\\1", contrast_full),
    DE_test = sub(".*_(\\w+)_(\\w+)$", "\\1", contrast_full),
    DE_direction = sub(".*_(\\w+)_(\\w+)$", "\\2", contrast_full),
    timepoint = sub(".*_(T0|T33).*", "\\1", level1),
    irrigation = sub(".*(drought|irrigated).*", "\\1", contrast),
    plant_part = sub(".*(leaf|root).*", "\\1", contrast),
    sig = ifelse(padj < 0.05, TRUE, FALSE)
    ) |>
  filter(!is.na(description),
         # Focusing only on the regular DE test results here, not the LFC-based ones
         DE_test == "DE") |>
  select(-p_adj_under, -level1, -level2)

# Check how many significant GO categories we have
GO_res |> count(contrast, DE_direction)


# T0 vs T33 at WD --------------------------------------------------------------
p_wd <- enrich_plot(GO_res,
                    contrasts = c("leaf_drought_T0_vs_leaf_drought_T33"),
                    DE_directions = c("up", "down"),
                    x_var = "DE_direction",
                    xlabs = c("up in\nNT33", "down in\nNT33"))
p_wd
ggsave(figfile_wd, p_wd, bg = "white", width = 8, height = 6, dpi = "retina")

p_wd2 <- enrich_plot(GO_res,
                    contrasts = c("leaf_drought_T0_vs_leaf_drought_T33"),
                    DE_directions = c("up", "down"),
                    x_var = "DE_direction",
                    xlabs = c("up in\nNT33", "down in\nNT33"),
                    n_tres = 2,
                    xlabsize = 10)
p_wd2
ggsave(figfile_wd2, p_wd2, bg = "white", width = 7, height = 3.5, dpi = "retina")


# DROUGHT > IRRIGATED --------------------------------------------------
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")

p_up <- enrich_plot(GO_res,
                    contrasts = contrast_sel,
                    DE_directions = c("up"),
                    n_tres = 2, padj_tres = 0.01,
                    xlabs = c("T0", "T33"))
p_up

ggsave(figfile_up, p_up, bg = "white", width = 8, height = 6)


# DROUGHT < IRRIGATED --------------------------------------------------
contrast_sel <- c("leaf_drought_T0_vs_leaf_irrigated_T0",
                  "leaf_drought_T33_vs_leaf_irrigated_T33")

p_down <- enrich_plot(GO_res,
                      contrasts = contrast_sel,
                      DE_directions = c("down"),
                      n_tres = 2, padj_tres = 0.01,
                      xlabs = c("T0", "T33"))
p_down

ggsave(figfile_down, p_down, bg = "white", width = 8, height = 8)
