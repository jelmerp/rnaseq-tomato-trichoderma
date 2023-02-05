# GENERAL SET-UP ---------------------------------------------------------------

## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "PCAtools")        # Running a PCA and plotting the results
pacman::p_load(char = packages)

## Load functions from separate script
source("scripts/PCA_fun.R")

## Set ggplot theme
theme_set(theme_bw())

## Define input files
count_table_file <- here("results/count/counts.txt")
metadata_file <- here("metadata/metadata.txt")

## Define output files
outdir <- here("results/PCA/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# PREPARE DATA -----------------------------------------------------------------

## Constants
SAMPLE_SELECTION <- c("S01", "S02", "S03", "S05", "S06", "S07", "S08", "S09",
                      "S10", "S11", "S12", "S14", "S15", "S16", "S17", "S18",
                      "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26",
                      "S27", "S28", "S29", "S30", "S31", "S32", "S33", "S34",
                      "S35", "S36", "S37", "S38", "S39", "S40")

## Load gene count data
raw_counts <- read.table(count_table_file, sep = "\t", header = TRUE, skip = 1)
my_regex <- ".+tomato.(.+)_Aligned.+"
colnames(raw_counts) <- sub(my_regex, "\\1", colnames(raw_counts))
counts <- raw_counts[, 7:ncol(raw_counts)] # Exclude metadata columns
rownames(counts) <- raw_counts$Geneid

## Load metadata
metadata <- read.table(metadata_file, sep = "\t", header = TRUE)
metadata <- filter(metadata, SampleID %in% SAMPLE_SELECTION)
counts <- counts[, SAMPLE_SELECTION]

matching_names <- identical(metadata$SampleID, colnames(counts))
if(matching_names == FALSE) stop("Sample ID in metadata and count matrix do not match!")

## Create the DESeq object
dds_raw <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ 1)


## Create subsets
dds_drought <- dds_raw[, dds_raw$Irrigation == "drought"]
dds_irrigated <- dds_raw[, dds_raw$Irrigation == "irrigated"]

dds_leaf <- dds_raw[, dds_raw$Part == "Leaf"]
dds_root <- dds_raw[, dds_raw$Part == "Root"]

dds_drought_leaf <- dds_drought[, dds_drought$Part == "Leaf"]
dds_drought_root <- dds_drought[, dds_drought$Part == "Root"]

dds_irrigated_leaf <- dds_irrigated[, dds_irrigated$Part == "Leaf"]
dds_irrigated_root <- dds_irrigated[, dds_irrigated$Part == "Root"]


# RUN AND PLOT THE PCA ---------------------------------------------------------

## Select dataset
dds <- dds_leaf

## Run the PCA
vsd <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
pca_res <- pca(vsd, metadata = colData(dds), removeVar = 0.1)
percent_var <- round(pca_res$variance, 2)

## Create a regular plot with PC1 and PC2
p <- biplot(pca_res,
            lab = NULL,
            pointSize = 4,
            colby = "Treatment",
            shape = "Irrigation") +
  theme_bw()
p

## Create a regular plot with PC1 and PC3 (could do any combination of PCs)
p <- biplot(pca_res,
            x = "PC1",
            y = "PC3",
            lab = NULL,
            pointSize = 4,
            colby = "Treatment",
            shape = "Irrigation") +
  theme_bw()
p

## Create a biplot
p <- biplot(pca_res,
            showLoadings = TRUE,
            ntopLoadings = 5,    # Number of genes to plot (multiplied by 2)
            colby = "Treatment",
            shape = "Irrigation") +
  xlab(paste0("PC1 (", percent_var[1], "% of variation)")) +
  ylab(paste0("PC2 (", percent_var[2], "% of variation)")) +
  theme_bw()
p

## Create a screeplot
p <- screeplot(pca_res) +
  ggtitle(NULL)
p

# Save a plot
pca_plotfile <- paste0(outdir, "testplot.png")
ggsave(pca_plotfile, p, width = 8, height = 7)


# ALTERNATIVELY, USE CUSTOM FUNCTIONS FROM PCA_fun.R ---------------------------

# Defaults for the run_pca() function:
# X-axis shows PC1 (`x = "PC1"`) => e.g. use `x = "PC3"` to show PC3
# Y-axis shows PC2 (`y = "PC2"`)
# Point color varies by Treatment (`color_var = "Treatment"`) => to omit varying colors, use `color_var = NULL`
# Point shape varies by Irrigation (`shape_var = "Irrigation"`) => to omit varying shapes, use `color_var = NULL`
# No plot title (`title = NULL`) => e.g. use `title = "My plot title"` to change
# No sample names added to points (`add_names = FALSE`) => Use `add_names = TRUE` to add sample names

# Use all defaults:
pca_res <- pca_run(dds_leaf)
p <- pca_plot(pca_res)
p <- pca_plot(pca_res, x = "PC3", y = "PC4")  # Plot PC4 and PC4
p <- pca_plot(pca_res, add_names = TRUE)
p <- pca_biplot(pca_res, n_genes = 10)
p <- pca_biplot(pca_res, n_genes = 20)

# Plot irrigated leaf only, PC3 and PC4, with no shape variation in points (since it"s leaf only)
pca_res <- pca_run(dds_irrigated_leaf)
p <- pca_plot(pca_res, shape_var = NULL, add_names = TRUE)
p <- pca_biplot(pca_res, shape_var = NULL)
