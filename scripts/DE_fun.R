plotgene <- function(geneID, dds, my_contrast) {

  d <- plotCounts(dds,
                  gene = geneID,
                  intgroup = "group",
                  returnData = TRUE)

  d$focal <- d$group %in% my_contrast

  p <- ggplot(d,
              aes(x = group, y = count, color = focal)) +
    scale_color_manual(values = c("black", "red")) +
    geom_point(position = position_jitter(width = 0.1, height = 0)) +
    labs(title = geneID) +
    theme_bw() +
    theme(legend.position = "none")

  print(p)
}


shrink_lfc <- function(dds,
                       contrast,
                       for_MA = TRUE,
                       lfc_thres = 0,
                       p_thres = 0.1,
                       s_thres = 0.005) {

  coef <- paste0("group_", paste0(contrast, collapse = "_vs_"))

  res <- lfcShrink(dds,
                   coef = coef,
                   type = "apeglm",
                   lfcThreshold = lfc_thres)

  if (for_MA == FALSE) {
    res <- res |>
      as.data.frame() |>
      rownames_to_column("gene_id") |>
      mutate(level1 = contrast[1],
            level2 = my_contrast[2])
  }

  # Report
  if (lfc_thres == 0) {
    n_sig <- sum(res$padj < p_thres, na.rm = TRUE)
    cat("Number of significant genes (p-value):", n_sig, "\n")
  } else {
    n_sig <- sum(res$svalue < s_thres, na.rm = TRUE)
    cat("Number of significant genes (s-value):", n_sig, "\n")
  }

  return(res)
}

plot_MA <- function(deseq_results,
                    p_thres = 0.1,
                    ptsize_nonsig = 1,
                    ptsize_sig = 1,
                    alpha_nonsig = 0.3,
                    alpha_sig = 0.3,
                    x_min = NA,
                    interactive = FALSE) {

  d <- plotMA(deseq_results,
              alpha = p_thres,
              returnData = TRUE) |>
    rownames_to_column("gene")

  p <- ggplot() +
    geom_point(data = filter(d, isDE == FALSE),
               aes(x = mean, y = lfc),
               shape = 21, fill = "grey50", color = "grey50",
               size = ptsize_nonsig, alpha = alpha_nonsig) +
    geom_point(data = filter(d, isDE == TRUE),
               aes(x = mean, y = lfc),
               shape = 21, fill = "blue", color = "blue",
               size = ptsize_sig, alpha = alpha_sig) +
    scale_x_log10(labels = scales::comma,
                  breaks = c(1, 10, 100, 1000, 10000, 100000)) +
    coord_cartesian(xlim = c(x_min, NA)) +
    geom_hline(yintercept = 0) +
    guides(color = FALSE) +
    labs(x = "Mean of normalized counts",
         y = "Log2-fold Change") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))

  if (interactive == TRUE) p <- ggplotly(p, tooltip = "text")

  return(p)
}

plot_heatmap <- function(geneIDs, dds) {

  ntd <- assay(normTransform(dds))

  ntd_sel <- ntd[match(geneIDs, rownames(ntd)), ]
  df_meta <- as.data.frame(colData(dds)[, c("Irrigation", "Treatment")])

  pheatmap(ntd_sel,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           annotation_col = df_meta)
}

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
