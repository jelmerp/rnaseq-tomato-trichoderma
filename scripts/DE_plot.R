# SETUP ------------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "patchwork",
              "glue",
              "VennDiagram")
pacman::p_load(char = packages)

# Load function from separate script
source("mcic-scripts/rnaseq/rfuns/DE_funs.R")

# Input files
dds_file <- "results/DE/dds_leaf.rds"
gff_file <- "data/ref/tomato/annot/ITAG4.1/ITAG4.1_gene_models.gff"

# Output files
outdir <- here("results/figures/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

figfile_MA <- here(outdir, "Fig2_MA_noshrink.png")
figfile_MA_shrink <- here(outdir, "Fig2_alt_MA_shrink.png")
figfile_volc <- here(outdir, "Fig2_alt_volcano.png")
figfile_volc_freeY <- here(outdir, "Fig2_alt_volcano_freeY.png")
figfile_heatmap <- here(outdir, "SuppFigX_heatmap.png")
figfile_venn <- file.path(outdir, "Fig4_venn.tiff")

# Set ggplot theme
theme_set(theme_bw(base_size = 15))
contrast_cols <- RColorBrewer::brewer.pal(3, "Dark2")

# Load the DESeq object
dds <- readRDS(dds_file)

# Get annotation
annot <- ape::read.gff(gff_file) |>
  filter(type == "mRNA",
         grepl("Note=", attributes)) |>
  transmute(gene = sub(".*Name=(Solyc\\w+\\.\\d).*", "\\1", attributes),
            description = sub(".*Note=([^;]+).*", "\\1", attributes)) |>
  mutate(description = str_trunc(description, width = 50)) |>
  arrange(gene)


# MA-PLOTS (NON-SHRUNKEN LFC) --------------------------------------------------
# Line comparison
res_WD <- extract_DE(comp = c("drought_T0", "drought_T33"), fac = "group",
                    dds = dds, p_tres = 0.1, lfc_tres = 0)

pa <- pMA(res_WD, x_min = 0.1, rm_padj_na = FALSE,
          ptsize_sig = 1.5, alpha_sig = 1, ptcol_sig = contrast_cols[1]) +
  labs(title = "T0 vs T33 for WD") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

# Drought at T0
res_t0 <- extract_DE(comp = c("drought_T0", "irrigated_T0"), fac = "group",
                     dds = dds, p_tres = 0.1, lfc_tres = 0)

pb <- pMA(res_t0, x_min = 0.1, rm_padj_na = FALSE,
          ptcol_sig = contrast_cols[2]) +
  labs(title = "WD vs Irr. at T0") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Drought at T33
res_t33 <- extract_DE(comp = c("drought_T33", "irrigated_T33"), fac = "group",
                     dds = dds, p_tres = 0.1, lfc_tres = 0)

pc <- pMA(res_t33, x_min = 0.1, rm_padj_na = FALSE,
          ptcol_sig = contrast_cols[3]) +
  labs(title = "WD vs Irr. at T33") +
  theme(axis.title.y = element_blank())

# Combined plot
p <- pa / pb / pc
ggsave(figfile_MA, p,
       width = 8, height = 8, dpi = "retina")


# MA-PLOTS (SHRUNKEN LFC) ------------------------------------------------------
# Line comparison
lfc_WD <- shrink_lfc(dds,
                     fac = "group",
                     comp = c("drought_T33", "drought_T0"),
                     p_tres = 0.1, lfc_tres = 0)
pa <- pMA(lfc_WD,
          x_min = 0.5,
          alpha_nonsig = 1, alpha_sig = 1,
          ptsize_nonsig = 2, ptsize_sig = 2)

# Drought at T0
lfc_t0 <- shrink_lfc(dds,
                     fac = "group",
                     comp = c("irrigated_T0", "drought_T0"),
                     p_tres = 0.1, lfc_tres = 0)
pb <- pMA(lfc_t0, x_min = 0.5)

# Drought at T33
lfc_t33 <- shrink_lfc(dds,
                      fac = "group",
                      comp = c("irrigated_T33", "drought_T33"),
                      p_tres = 0.1, lfc_tres = 0)
pc <- pMA(lfc_t33, x_min = 0.5)

# Combined plot
p <- pa / pb / pc
ggsave(figfile_MA_shrink, p, width = 7, height = 7, dpi = "retina")


# VOLCANO PLOTS ----------------------------------------------------------------
res <- bind_rows(res_t0, res_t33, res_WD) |>
  mutate(contrast = case_when(
    contrast == "drought_T0_irrigated_T0" ~ "WD vs Irr. at T0",
    contrast == "drought_T33_irrigated_T33" ~ "WD vs Irr. at T33",
    contrast == "drought_T0_drought_T33" ~ "T0 vs T33 for WD"
  ))

p_volc1 <- pvolc(res, contrasts = "all", facet_scales = "fixed", sig_only = FALSE)
ggsave(figfile_volc, p_volc1, width = 8, height = 7, dpi = "retina")

p_volc2 <- pvolc(res, contrasts = "all", facet_scales = "free_y", sig_only = FALSE)
ggsave(figfile_volc_freeY, p_volc2, width = 8, height = 7, dpi = "retina")


# BOXPLOTS ---------------------------------------------------------------------
genes <- res_WD |> filter(padj < 0.1) |> pull(gene)
counts_norm <- norm_counts(dds) |>
  left_join(annot, by = "gene") |>
  mutate(Irrigation = ifelse(Irrigation == "drought", "WD", "Irr"))

# Plot 1 gene
#pbox(gene = genes[2], count_df = counts_norm, annot = annot,
#     x_by = "Irrigation", col_by = "Treatment")

# Plot all WD-significant genes
walk(genes, pbox,
     count_df = counts_norm, annot = annot,
     x_by = "Irrigation", col_by = "Treatment",
     save_plot = TRUE)


# HEATMAP ----------------------------------------------------------------------
genes <- res_WD |> filter(padj < 0.1) |> pull(gene)
count_mat <- norm_counts(dds, return_matrix = TRUE)
meta <- as.data.frame(colData(dds)) |>
  mutate(Irrigation = ifelse(Irrigation == "drought", "WD", "Irr"))

p <- pheat(genes, count_mat, meta, groups = c("Irrigation", "Treatment"))

ggsave(figfile_heatmap, p, width = 8, height = 7, dpi = "retina")


# VENN DIAGRAM -----------------------------------------------------------------
genes_t0 <- res_t0 |> filter(isDE == TRUE) |> pull(gene)
genes_t33 <- res_t33 |> filter(isDE == TRUE) |> pull(gene)

venn.diagram(
  x = list(genes_t0, genes_t33),
  category.names = c("T0", "T33"),
  filename = figfile_venn,
  output = TRUE,
  fill = contrast_cols[2:3],
  cat.cex = 2.25,
  cat.dist = 0.05,
  cat.pos = c(-40, 40),
  cex = 1.75,
  cat.fontface = "bold",
  ext.text = TRUE,
  ext.percent = c(1, 1, 0.1),
  ext.line.lwd = 2,
  ext.dist = c(-0.1, -0.1),
  ext.length = c(0.75, 0.75),
  ext.pos = c(-30, 30),
  margin = 0.05
)
