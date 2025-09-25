#' @title enrichGroup
#' @description Enrichment for genes in specific pattern
#' @details This function generates an enrichResult instance
#' @param gene.group a data frame indicate genes in each pattern
#' @param species select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param ... pass to the function enrichGO
#' @return A list of enrichResult instance
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export enrichGroup
#' @import org.Hs.eg.db org.Mm.eg.db
#' @importFrom clusterProfiler enrichGO
enrichGroup <- function(gene.group, species="mouse", ont="BP", ...) {
  OrgDb <- ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db")
  enrich_result <- list()
  for (pattern in unique(gene.group$pattern)){
    genes <- rownames(gene.group)[which(gene.group$pattern==pattern)]
    enrich_result[[pattern]] <- clusterProfiler::enrichGO(gene=genes, OrgDb=OrgDb, keyType="SYMBOL", ont=ont, ...)
  }
  return(enrich_result)
}

#' @title enrichPattern
#' @description Flexible GO enrichment for one, multiple, or all temporal patterns
#' @details This function performs enrichment analysis for user-specified patterns.
#'          If no pattern is provided, all patterns in gene.group will be analyzed
#' @param gene.group A data frame indicating genes and their assigned temporal patterns
#' @param patterns Character vector of patterns to analyze (e.g., c("I","D")),
#'        or NULL to analyze all patterns
#' @param species "mouse" or "human"
#' @param ont One of "BP", "MF", "CC", or "ALL"
#' @param universe Background set of genes (default: NULL → auto-detected)
#' @param ... Additional arguments passed to enrichment function
#' @return A named list of enrichResult objects, one per pattern
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export enrichPattern
#' @import org.Hs.eg.db org.Mm.eg.db
#' @importFrom clusterProfiler enrichGO
enrichPattern <- function(gene.group, patterns = NULL, species = "mouse", ont = "BP", universe = NULL, ...) {
  OrgDb <- ifelse(species == "mouse", "org.Mm.eg.db", "org.Hs.eg.db")

  # If no pattern is specified, use all
  if (is.null(patterns)) {
    patterns <- unique(gene.group$pattern)
  }

  enrich_result_list <- list()
  for (pattern in patterns) {
    genes <- rownames(gene.group)[gene.group$pattern == pattern]
    enrich_result <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = OrgDb,
      keyType = "SYMBOL",
      ont = ont,
      universe = universe,
      ...
    )
    # Compute EnrichRatio
    if (!is.null(enrich_result@result) && nrow(enrich_result@result) > 0) {
      enrich_result@result$EnrichRatio <- with(enrich_result@result, {
        (as.numeric(sub("/.*", "", GeneRatio)) * as.numeric(sub(".*/", "", BgRatio))) /
          (as.numeric(sub(".*/", "", GeneRatio)) * as.numeric(sub("/.*", "", BgRatio)))
      })
    }
    enrich_result_list[[pattern]] <- enrich_result
  }
  return(enrich_result_list)
}


#' @title compareEnrichBin
#' @description  Bin-based enrichment for ordered genes in a specific pattern
#' @details This function partitions genes within a temporal expression pattern into bins
#'          along pseudotime and performs GO enrichment on each bin. The output contains
#'          enrichment results for all bins
#' @param gene.group A data frame indicating genes in each pattern
#' @param pattern The expression pattern for enrichment analysis
#' @param bin.width The width of each bin
#' @param stride Stride of each step
#' @param species Select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param universe Pass to the universe paramenter of enrichGO
#' @param ... Pass to the function enrichGO
#' @return A compareClusterResult object summarizing enrichment results across bins
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export compareEnrichBin
#' @import org.Hs.eg.db org.Mm.eg.db dplyr
#' @importFrom clusterProfiler compareCluster
compareEnrichBin <- function(gene.group, pattern, bin.width=0.2, stride=0.1, species="human", ont="BP", universe=FALSE, ...){
  OrgDb <- ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db")
  genes <- rownames(gene.group)[which(gene.group$pattern==pattern)]
  bin.width <- ifelse(bin.width<1, as.integer(bin.width*length(genes)), bin.width)
  stride <- ifelse(stride<1, as.integer(stride*length(genes)), stride)
  gene_list <- list()
  gene_pos <- c()
  pos_start <- 1
  pos_end <- 1
  if (length(universe) == 1) {
    universe <- genes
  }
  while (pos_end < length(genes)) {
    pos_end <- pos_start + bin.width - 1
    if(pos_end > length(genes)) {
      if (pos_end - length(genes) < 0.5*stride){
        ### if the last bin is long enough, then isolate it out
        pos_end <- length(genes)
        pos_start <- length(genes) - bin.width + 1
      }
      else {
        ### if the last bin is not long enough, merge
        pos_end <- length(genes)
        pos_start <- pos_start - stride
      }
    }
    gene_list[[as.character(pos_start)]] <- genes[pos_start:pos_end]
    gene_pos[[as.character(pos_start)]] <- paste0(pos_start, "-", pos_end)
    pos_start <- pos_start + stride
  }
  names(gene_list) <- unlist(gene_pos)
  genes_bin_enrich <- clusterProfiler::compareCluster(gene_list, fun = "enrichGO", OrgDb = OrgDb, universe=universe, keyType = "SYMBOL", ont="BP", pvalueCutoff = 1, qvalueCutoff = 1, ...)
  genes_bin_enrich@compareClusterResult[, "Cluster"] <- factor(genes_bin_enrich@compareClusterResult[, "Cluster"], levels=unlist(gene_pos))
  ## Compute Enrichratio
  genes_bin_enrich@compareClusterResult$EnrichRatio <- with(genes_bin_enrich@compareClusterResult, {
    (as.numeric(sub("/.*", "", GeneRatio)) * as.numeric(sub(".*/", "", BgRatio))) /
      (as.numeric(sub(".*/", "", GeneRatio)) * as.numeric(sub("/.*", "", BgRatio)))
  })
  return(genes_bin_enrich)
}

