# SET-UP -----------------------------------------------------------------------

## Install/load packages
if(! 'pacman' %in% installed.packages()) install.packages('pacman')
packages <- c('tidyverse',       # Misc. data manipulation and plotting
              'here',            # Managing file paths
              'VennDiagram',     # Venn Diagrams
              "RColorBrewer")    # Plot colors
pacman::p_load(char = packages)

## Source script with functions
source("scripts/GO/GO_fun.R")

## Output dir
outdir <- "results/DE/overlaps"
if(!dir.exists(outdir)) dir.create(outdir)

## Input files
DE_dir <- "results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE,
                       pattern = ".*_T\\d+_DEgenes.txt")
names(DE_files) <- sub("_DEgenes.txt", "", basename(DE_files))

## Read DE results
GO_res <- map_dfr(DE_files, read_tsv,
                  col_types = "cnnnnnncc",
                  .id = "contrast_full")

## Get DE genes for the different comparisons
D0_I0 <- GO_res %>%
  filter(contrast_full == "leaf_drought_T0_vs_leaf_irrigated_T0") %>%
  pull(gene_id)
D33_I33 <- GO_res %>%
  filter(contrast_full == "leaf_drought_T33_vs_leaf_irrigated_T33") %>%
  pull(gene_id)
D0_D33 <- GO_res %>%
  filter(contrast_full == "leaf_drought_T0_vs_leaf_drought_T33") %>%
  pull(gene_id)
# I0_I33 => No DE genes

## Intersect DE genes: get genes that are in multiple comparisons
DI <- generics::intersect(D0_I0, D33_I33)
length(DI)   # 15,255

all <- generics::intersect(D0_D33, DI)
length(all)  # 16

## Make Venn Diagram:
my_cols <- brewer.pal(3, "Pastel2")

venn_plot <- venn.diagram(
  x = list(D0_I0, D33_I33, D0_D33),
  category.names = c('WD vs IR at T0', 'WD vs IR at T33', 'T0 vs T33 at WD'),
  filename = here(outdir, "venn1.tiff"),
  output = TRUE,
  fill = my_cols,
  cat.cex = 1,
  cat.fontface = "bold",
  cat.pos = c(-35, 35, 0),
  cat.dist = c(0.05, 0.05, -0.02)
)
grid::grid.newpage()
grid::grid.draw(venn_plot)

venn.diagram(
  x = list(D0_I0, D33_I33),
  category.names = c('drought0_irr0', 'drought33_irr33'),
  filename = here(outdir, "venn2.tiff"),
  output = TRUE,
  fill = my_cols[1:2],
  cat.cex = 0.8,
  cat.fontface = "bold"
)
