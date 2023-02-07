## Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
get_pathway <- function(K_term, outdir) {
  # K_term <- "K13449"

  cat("K_term:", K_term, " ")

  tryCatch(
    {
      kegg_info <- keggGet(K_term)
      pathway_df <- data.frame(pathway_description = kegg_info[[1]]$PATHWAY) %>%
        rownames_to_column("pathway_id") %>%
        mutate("K_term" = K_term)

      cat(" Nr of pathways:", nrow(pathway_df), "\n")

      if(nrow(pathway_df) > 0) {
        pathway_df_file <- file.path(outdir, paste0(K_term, ".txt"))
        write_tsv(pathway_df, pathway_df_file)
        return(pathway_df)
      } else {
        return(NULL)
      }
    },
    error = function(cond) {
      message("keggGet failure")
      message(cond)
      return(NULL)
    }
  )
}

## Function to run a KEGG enrichment analysis for a certain DE output file
run_kegg <- function(comp, DE_direction, test_type = "DE",
                     DE_dir, kegg_map) {

  if (test_type == "DE") LFC_arg <- FALSE else LFC_arg <- TRUE
  DE_res <- get_DE_res(comp, DE_dir, LFC = LFC_arg)

  if (DE_direction == "up") DE_res <- DE_res %>% filter(log2FoldChange > 0)
  if (DE_direction == "down") DE_res <- DE_res %>% filter(log2FoldChange < 0)

  if (test_type == "DE") DE_res <- DE_res %>% filter(padj < 0.1)
  if (test_type == "LFC") DE_res <- DE_res %>% filter(svalue < 0.005)

  DE_genes <- DE_res$gene_id

  cat("---------------------\n")
  cat(comp[1], "vs", comp[2], ":", DE_direction, "in", comp[1],
      "- test type:", test_type, "\n")
  cat("Number of DE genes:", length(DE_genes), "\n")

  if (length(DE_genes) == 0) return(NULL)

  kegg_res <- enricher(DE_genes, TERM2GENE = kegg_map)
  kegg_res <- as.data.frame(kegg_res) %>%
    mutate(level1 = comp[1],
           level2 = comp[2],
           DE_direction = DE_direction,
           test_type = test_type) %>%
    relocate(geneID, .after = test_type)

  cat("Number of enriched pathways:", nrow(kegg_res), "\n")
  return(kegg_res)
}
