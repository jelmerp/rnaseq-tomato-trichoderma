## Read and process the DE results for a given comparison/contrast
get_DE_res <- function(comp, DE_dir, LFC = FALSE) {

  comp <- paste0(comp, collapse = "_vs_")

  if (LFC == FALSE) {

    DE_file <- file.path(DE_dir, paste0(comp, "_all-genes.txt"))
    DE_res <- read_tsv(DE_file, col_types = cols()) %>%
      filter(!is.na(padj)) %>%
      arrange(gene_id)

  } else {

    DE_file <- file.path(DE_dir, paste0(comp, "_LFC5_all-genes.txt"))
    DE_res <- read_tsv(DE_file, col_types = cols()) %>%
      filter(!is.na(svalue)) %>%
      arrange(gene_id)
  }

  return(DE_res)
}

## Create a "DE vector" for the goseq function: named vector of 0s and 1s indicating significance
get_DE_vec <- function(DE_res, LFC = FALSE,
                       p_threshold = 0.1, s_threshold = 0.005) {

  comp <- paste0(comp, collapse = "_vs_")

  if (LFC == FALSE) {

    DE_vec <- ifelse(DE_res$padj < p_threshold, 1, 0)
    names(DE_vec) <- DE_res$gene_id

  } else {

    DE_vec <- ifelse(DE_res$svalue < s_threshold, 1, 0)
    names(DE_vec) <- DE_res$gene_id
  }

  cat("Number of significant/total genes:", sum(DE_vec), "/", length(DE_vec), "\n")

  return(DE_vec)
}

## Get the top n (max_genes) DE genes in each GO category as a string:
get_GO_genes <- function(GO_cat, DE_res, GO_map, max_genes = 50) {
  # GO_cat <- "GO:0005886"

  ## Get all DE genes ordered by significance:
  DE_genes <- DE_res %>%
    arrange(padj) %>%
    filter(padj < 0.1) %>%
    pull(gene_id)

  ## Get DE genes in focal GO category:
  DE_in_cat <- GO_map %>%
    filter(GO_term == GO_cat,
           gene_id %in% DE_genes) %>%
    pull(gene_id)

  ## Get top n ("max_genes") most significantly DE genes:
  n_genes_take <- ifelse(length(DE_in_cat) > max_genes, max_genes, length(DE_in_cat))
  DE_in_cat_sel <- DE_in_cat[order(match(DE_in_cat, DE_genes))][1:n_genes_take]

  ## Convert to a single string:
  DE_in_cat_string <- paste0(DE_in_cat_sel, collapse = " / ")

  return(DE_in_cat_string)
}

## Function to run GO analysis
run_GO <- function(comp_name, DE_res, GO_map, gene_length,
                   outdir, file_suffix = NULL) {

  DE_vec <- get_DE_vec(DE_res, LFC = FALSE)

  if (sum(DE_vec > 0)) {
    ## Remove rows from gene length df not in the DE_vec
    gene_length <- gene_length %>%
      filter(gene_id %in% names(DE_vec))

    ## Remove elements from DE_vec not among the gene lengths
    DE_vec <- DE_vec[names(DE_vec) %in% gene_length$gene_id]

    ## Probability weighting function based on gene lengths
    pwf <- nullp(
      DEgenes = DE_vec,
      bias.data = gene_length$length,
      plot.fit = FALSE
    )

    ## Run GO test
    GO_res <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")

    ## Process GO results
    GO_res <- GO_res %>%
      filter(numDEInCat > 0) %>%    # P-adjustment only for genes that were actually tested
      mutate(
        p_adj_over = p.adjust(over_represented_pvalue, method = "BH"),
        p_adj_under = p.adjust(under_represented_pvalue, method = "BH"),
        level1 = comp_name[1],
        level2 = comp_name[2]
      ) %>%
      select(p_adj_over, p_adj_under, numDEInCat, numInCat,
             category, ontology, term, level1, level2) %>%
      filter(p_adj_over < 0.05) %>%      # Only keep significant categories
      as_tibble()

    ## Get DE genes in each GO category
    GO_res$genes_DE_in_cat <- sapply(GO_res$category,
                                     get_GO_genes,
                                     DE_res, GO_map)

    ## Report
    cat("\n--------------\nRan comparison for:", comp_name, "\n")
    cat("Number of significant DE genes:", sum(DE_vec), "\n")
    cat("Number of significant GO categories:", nrow(GO_res), "\n")

    ## Write results to file
    comp_name <- paste0(comp_name, collapse = "_vs_")
    outfile <- file.path(outdir, paste0(comp_name, file_suffix, "_GO.txt"))
    write_tsv(GO_res, outfile)

    return(GO_res)
  } else {
    cat("No significant GO categories\n")
  }
}
