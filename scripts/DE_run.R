# SETUP ------------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "apeglm")
pacman::p_load(char = packages)

# Load the data
count_table_file <- here("results/count/counts.txt")
metadata_file <- here("metadata/metadata.txt")

# Create output dirs
outdir <- here("results/DE/")
plotdir <- here("results/DE/fig/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)


# FUNCTIONS --------------------------------------------------------------------
sig_contrast <- function(my_contrast, dds, outdir, level_prefix) {

  contrast_full <- paste0(level_prefix, "_", my_contrast)

  ## Write all results:
  res <- results(dds,
                 contrast = c("group", my_contrast)) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    mutate(level1 = contrast_full[1],
           level2 = contrast_full[2])

  contrast_pasted <- paste0(contrast_full, collapse = "_vs_")
  outfile <- file.path(outdir, paste0(contrast_pasted, '_all-genes.txt'))
  write_tsv(res, outfile)

  ## Write significant results:
  res_sig <- res %>% dplyr::filter(padj < 0.1)

  cat(my_contrast[1], "versus", my_contrast[2], ":", nrow(res_sig), "significant\n")

  if(nrow(res_sig > 0)) {
    contrast_pasted <- paste0(contrast_full, collapse = "_vs_")
    outfile <- file.path(outdir, paste0(contrast_pasted, '_DEgenes.txt'))
    write_tsv(res_sig, outfile)
  }

  return(res_sig)
}


# PREPARE DATA -----------------------------------------------------------------
# Load input data
raw_counts <- read.table(count_table_file,
                         sep = "\t", header = TRUE, skip = 1)
my_regex <- ".+tomato.(.+)_Aligned.+"
colnames(raw_counts) <- sub(my_regex, "\\1", colnames(raw_counts))
counts <- raw_counts[, 7:ncol(raw_counts)]
rownames(counts) <- raw_counts$Geneid

# Load the metadata
metadata <- read.table(file = metadata_file, sep = "\t", header = TRUE)

matching_names <- identical(metadata$SampleID, colnames(counts))
if (matching_names == FALSE) stop("Sample ID in metadata and count matrix do not match!")

# Create the DESeq object
dds_raw <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ 1)

# Remove S04 which has very few reads
dds_raw <- dds_raw[, dds_raw$SampleID != "S04"]
dds_raw <- dds_raw[, dds_raw$SampleID != "S13"]

# Subset to leaf
dds_leaf <- dds_raw[, dds_raw$Part == "Leaf"]
dds <- dds_leaf
level_prefix <- "leaf"

# Set the analysis design:
dds$group <- factor(paste(dds$Irrigation , dds$Treatment, sep = "_"))
dds$group <- relevel(dds$group, ref = "drought_T0")
design(dds) <- formula(~ group)
dds <- DESeq(dds)

# Save RDS
dds_leaf_file <- "results/DE/dds_leaf.rds"
saveRDS(dds, dds_leaf_file)


# RUN DE FOR ALL COMPARISONS ---------------------------------------------------
comps <- combn(levels(dds@colData$group), 2, simplify = FALSE)

# Default comparison with adjusted p-value <0.1
sig_all_contrasts <- map_dfr(.x = comps, .f = sig_contrast,
                             dds, outdir, level_prefix)

# LFC-based threshold after LFC shrinkage
sig_all_contrasts <- map_dfr(.x = comps, .f = sig_contrast,
                             dds, shrunken_lfc = TRUE, lfc_threshold = 5,
                             outdir, level_prefix)
