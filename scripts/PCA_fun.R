## Function to run the PCA using the PCAtools `pca()` function
pca_run <- function(dds,
                    remove_prop = 0.1) {   # Remove this proportion of variables (genes) with the lowest variance
  vsd <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  pca_res <- pca(vsd, metadata = colData(dds), removeVar = remove_prop)
  return(pca_res)
}

## Function to create a regular PCA plot
pca_plot <- function(
  pca_res,                      # PCA results after running `PCA_run()`
  color_var = "Treatment",      # Vary point color by this variable from the metadata
  shape_var = "Irrigation",     # Vary point shape by this variable from the metadata
  x = "PC1",                    # Principal component to plot on the x-axis
  y = "PC2",                    # Principal component to plot on the y-axis
  add_names = TRUE,             # Add sample names to points TRUE/FALSE
  title = NULL)                 # Add a title as a string; "NULL" is no title
{

  percent_var <- round(pca_res$variance, 2)

  if (add_names == TRUE) names <- pca_res$yvars
  if (add_names == FALSE) names <- 0

  x_nr <- as.integer(sub("PC", "", x))
  y_nr <- as.integer(sub("PC", "", y))

  p <- biplot(pca_res,               # `biplot()` function from PCAtools
              x = x,
              y = y,
              selectLab = names,
              colby = color_var,
              shape = shape_var) +
    xlab(paste0(x, " (", percent_var[x_nr], "% of variance)")) +
    ylab(paste0(y, " (", percent_var[y_nr], "% of variance)")) +
    ggtitle(title) +
    theme_bw()

  print(p)
  return(p)
}

## Function to create a PCA biplot
pca_biplot <- function(
  pca_res,                      # PCA results after running `PCA_run()`
  color_var = "Treatment",      # Vary point color by this variable from the metadata
  shape_var = "Irrigation",     # Vary point shape by this variable from the metadata
  x = "PC1",                    # Principal component to plot on the x-axis
  y = "PC2",                    # Principal component to plot on the y-axis
  n_genes = 5,                  # Number of genes to plot loadings for
  add_names = TRUE,             # Add sample names to points TRUE/FALSE
  title = NULL)                 # Add a title as a string; "NULL" is no title
{
  percent_var <- round(pca_res$variance, 2)

  if (add_names == TRUE) name_idx <- 1:nrow(pca_res$metadata)
  if (add_names == FALSE) name_idx <- 0

  x_nr <- as.integer(sub("PC", "", x))
  y_nr <- as.integer(sub("PC", "", y))

  p <- biplot(pca_res,         # `biplot()` function from PCAtools
              x = x,
              y = y,
              selectLab = name_idx,
              showLoadings = TRUE,
              ntopLoadings = n_genes / 2,
              colby = color_var,
              shape = shape_var) +
    xlab(paste0(x, " (", percent_var[x_nr], "% of variance)")) +
    ylab(paste0(y, " (", percent_var[y_nr], "% of variance)")) +
    ggtitle(title) +
    theme_bw()

  print(p)
  return(p)
}

## Function corresponding to the old way of making the PCA
## (The DEseq function used `prcomp()` under the hood)
pca_prcomp <- function(
  dds,
  color_var = 'Treatment',
  shape_var = 'Irrigation',
  x_PC = 1,
  y_PC = 2,
  ntop = 500,
  add_names = FALSE,
  title = NULL)
{

  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))

  percent_var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)

  d <- as.data.frame(cbind(colData(vsd), pca$x))

  p <- ggplot(d,
              aes_string(x = paste0("PC", x_PC),
                         y = paste0("PC", y_PC),
                         color = color_var,
                         shape = shape_var)) +
    geom_point(size = 6) +
    xlab(paste0("PC", x_PC, " (", percent_var[x_PC], "% of variance)")) +
    ylab(paste0("PC", y_PC, " (", percent_var[y_PC], "% of variance)")) +
    ggtitle(title)

  if(add_names == TRUE) {
    p <- p + geom_text_repel(aes(label = SampleID),
                             max.overlaps = 20,
                             point.padding = 3)
  }

  print(p)
  return(p)
}
