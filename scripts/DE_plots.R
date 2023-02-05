# SETUP ------------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "patchwork")
pacman::p_load(char = packages)

# Load function from separate script
source("scripts/DE_fun.R")

# Set ggplot theme
theme_set(theme_bw(base_size = 15))

# Load the data
dds_file <- "results/DE/old/dds_leaf.rds"

# Create output dirs
outdir <- here("results/figures_2023/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load the data
dds <- readRDS(dds_file)
resultsNames(dds)

# Run another analysis so we can access the T33 comparison for LFC shrinkage
dds2 <- dds
dds2$group <- relevel(dds2$group, ref = "drought_T33")
dds2 <- DESeq(dds2)


# EXPLORE OTHER RESULTS --------------------------------------------------------
# No significant genes for T0 vs 33 Irr
results(dds, contrast = c("group", "irrigated_T0", "irrigated_T33"), tidy = TRUE) |>
  filter(padj < 0.1)

# Check interaction - No significant results
ddsi <- dds
ddsi$Irrigation <- factor(ddsi$Irrigation, levels = c("irrigated", "drought"))
ddsi$Treatment <- factor(ddsi$Treatment, levels = c("T0", "T33"))
design(ddsi) <- formula(~ Irrigation + Treatment + Irrigation:Treatment)
ddsi <- DESeq(ddsi)
resultsNames(ddsi)
results(ddsi, name = "Irrigationdrought.TreatmentT33", tidy = TRUE) |>
  filter(padj < 0.1)


# MA-PLOTS (NON-SHRUNKEN LFC) --------------------------------------------------
# Line comparison
res_WD <- results(dds, contrast = c("group", c("drought_T0", "drought_T33")))
pa <- plot_MA(res_WD, x_min = 0.1,
              ptsize_sig = 2, alpha_sig = 1) +
  labs(title = "T0 vs T33 for WD") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Drought at T0
res_t0 <- results(dds, contrast = c("group", c("drought_T0", "irrigated_T0")))
pb <- plot_MA(res_t0, x_min = 0.1) +
  labs(title = "WD vs Irr. at T0") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Drought at T33
res_t33 <- results(dds, contrast = c("group", c("drought_T33", "irrigated_T33")))
pc <- plot_MA(res_t33, x_min = 0.1) +
  labs(title = "WD vs Irr. at T33")

# Combined plot
p <- pa / pb / pc
MA_file <- here(outdir, "MA_noshrink.png")
ggsave(MA_file, p,
       width = 7, height = 7, dpi = "retina")


# MA-PLOTS (SHRUNKEN LFC) ------------------------------------------------------
# Line comparison
lfc_WD <- shrink_lfc(dds, contrast = c("drought_T33", "drought_T0"))
pa <- plot_MA(lfc_WD, x_min = 0.5,
              alpha_nonsig = 1, alpha_sig = 1, ptsize_nonsig = 2, ptsize_sig = 2)

# Drought at T0
lfc_t0 <- shrink_lfc(dds, contrast = c("irrigated_T0", "drought_T0"))
pb <- plot_MA(lfc_t0, x_min = 0.5)

# Drought at T33
lfc_t33 <- shrink_lfc(dds2, contrast = c("irrigated_T33", "drought_T33"))
pc <- plot_MA(lfc_t33, x_min = 0.5)

# Combined plot
p <- pa / pb / pc
MA_shrink_file <- here(outdir, "MA_shrink.png")
ggsave(MA_shrink_file, p,
       width = 7, height = 7, dpi = "retina")


# VOLCANO PLOTS ----------------------------------------------------------------




# PLOT INDIV GENES -------------------------------------------------------------

