# GENERAL SET-UP ---------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "PCAtools",        # Running a PCA and plotting the results
              "patchwork")       # Combining plots
pacman::p_load(char = packages)

# Load functions from separate script
source("scripts/PCA_fun.R")

# Set ggplot theme
theme_set(theme_bw(base_size = 14))

# Define input files
count_table_file <- file.path("results/count/counts.txt")
metadata_file <- file.path("metadata/metadata.txt")

# Define output files
outdir <- file.path("results/figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figfile_all <- file.path(outdir, "Fig1_PCA.png")
figfile_by_irr <- file.path(outdir, "SuppFig5_PCA.png")

# Constants
SAMPLE_SELECTION <- c("S01", "S02", "S03", "S05", "S06", "S07", "S08", "S09",
                      "S10", "S11", "S12", "S14", "S15", "S16", "S17", "S18",
                      "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26",
                      "S27", "S28", "S29", "S30", "S31", "S32", "S33", "S34",
                      "S35", "S36", "S37", "S38", "S39", "S40")


# PREPARE DATA -----------------------------------------------------------------
# Load gene count data
raw_counts <- read.table(count_table_file, sep = "\t", header = TRUE, skip = 1)
my_regex <- ".+tomato.(.+)_Aligned.+"
colnames(raw_counts) <- sub(my_regex, "\\1", colnames(raw_counts))
counts <- raw_counts[, 7:ncol(raw_counts)] # Exclude metadata columns
rownames(counts) <- raw_counts$Geneid
counts <- counts[, SAMPLE_SELECTION]

# Load metadata
metadata <- read.table(metadata_file, sep = "\t", header = TRUE)
metadata <- metadata |> filter(SampleID %in% SAMPLE_SELECTION)

all_names_match <- identical(metadata$SampleID, colnames(counts))
if (all_names_match == FALSE) stop("Sample ID in metadata and count matrix do not match!")

# Create the DESeq object
dds_raw <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ 1)

# Rename irrigation levels
dds_raw$Irrigation <- ifelse(dds_raw$Irrigation == "drought",
                             "WD stress", "Irrigated")
dds_raw$Irrigation <- factor(dds_raw$Irrigation,
                             levels = c("Irrigated", "WD stress"))

# Create subsets
dds <- dds_raw[, dds_raw$Part == "Leaf"]
dds_drought <- dds[, dds$Irrigation == "WD stress"]
dds_irrigated <- dds[, dds$Irrigation == "Irrigated"]


# PCA WITH ALL LEAF SAMPLES ---------------------------------------------------
res_all <- pca_run(dds)
p <- pca_plot(res_all, color_var = "Treatment", shape_var = "Irrigation")

ggsave(figfile_all, p, width = 8, height = 7)


# SEPARATE PCAS FOR WD AND IRR -------------------------------------------------
res_irr <- pca_run(dds_irrigated)
res_wd <- pca_run(dds_drought)

p_irr <- pca_plot(res_irr, color_var = "Treatment")
p_wd <- pca_plot(res_wd, color_var = "Treatment")

p <- p_irr + p_wd +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20),
        legend.position = "top")

ggsave(figfile_by_irr, p, width = 10, height = 6.5)
