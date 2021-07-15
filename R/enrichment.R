#' @title enrichGroup
#' @description Enrichment for genes in specific pattern
#' @details This function generates an enrichResult instance
#' @param gene.group a data frame indicate genes in each pattern
#' @param species select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param ... pass to the function enrichGO
#' @return A list of enrichResult instance
#' @author Zhicheng Ji, Changxin Wan
#' @export enrichGroup
#' @import org.Hs.eg.db org.Mm.eg.db
#' @importFrom clusterProfiler enrichGO
enrichGroup <- function(gene.group, species="mouse", ont="BP", ...) {
  OrgDb <- ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db")
  enrich_result <- list()
  for (pattern in unique(gene.group$pattern)){
    genes <- rownames(gene.group)[which(gene.group$pattern==pattern)]
    enrich_result[[pattern]] <- enrichGO(gene=genes, OrgDb=OrgDb, keyType="SYMBOL", ont=ont, ...)
  }
  return(enrich_result)
}

#' @title enrichPattern
#' @description Enrichment for genes in specific pattern
#' @details This function generates an enrichResult instance
#' @param gene.group a data frame indicate genes in each pattern
#' @param pattern the pattern for enrichment analysis
#' @param species select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param ... pass to the function enrichGO
#' @return An enrichResult instance
#' @author Zhicheng Ji, Changxin Wan
#' @export enrichPattern
#' @import org.Hs.eg.db org.Mm.eg.db
#' @importFrom clusterProfiler enrichGO

enrichPattern <- function(gene.group, pattern, species="mouse", ont="BP", ...) {
  OrgDb <- ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db")
  genes <- rownames(gene.group)[which(gene.group$pattern==pattern)]
  enrich_result <- enrichGO(gene=genes, OrgDb=OrgDb, keyType="SYMBOL", ont=ont, ...)
  return(enrich_result)
}


#' @title compareEnrichBin
#' @description Enrichment for ordered genes in specific pattern
#' @details This function generates a list with genes and enrichResult instance
#' @param gene.group a data frame indicate genes in each pattern
#' @param pattern the pattern for enrichment analysis
#' @param bin.width the width of each bin
#' @param stride stride of each step
#' @param species select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param universe pass to the universe paramenter of enrichGO
#' @param ... pass to the function enrichGO
#' @return compareClusterResult instance
#' @author Zhicheng Ji, Changxin Wan
#' @export compareEnrichBin
#' @import org.Hs.eg.db org.Mm.eg.db dplyr clusterProfiler
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
  genes_bin_enrich <- compareCluster(gene_list, fun = "enrichGO", OrgDb = OrgDb, universe=universe, keyType = "SYMBOL", ont="BP", pvalueCutoff = 1, qvalueCutoff = 1, ...)
  genes_bin_enrich@compareClusterResult[, "Cluster"] <- factor(genes_bin_enrich@compareClusterResult[, "Cluster"], levels=unlist(gene_pos))
  return(genes_bin_enrich)
}



