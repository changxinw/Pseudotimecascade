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

#' @title enrichBin
#' @description Enrichment for ordered genes in specific pattern
#' @details This function generates a list with genes and enrichResult instance
#' @param gene.group a data frame indicate genes in each pattern
#' @param pattern the pattern for enrichment analysis
#' @param bin.width the width of each bin
#' @param stride stride of each step
#' @param species select from human or mouse
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
#' @param ... pass to the function enrichGO
#' @return A list with genes and enrichResult instance
#' @author Zhicheng Ji, Changxin Wan
#' @export enrichBin
#' @import org.Hs.eg.db org.Mm.eg.db
#' @importFrom clusterProfiler enrichGO
enrichBin <- function(gene.group, pattern, bin.width=0.2, stride=0.1, species="mouse", ont="BP", ...){
  OrgDb <- ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db")
  genes <- rownames(gene.group)[which(gene.group$pattern==pattern)]
  res_bin <- list()
  bin.width <- ifelse(bin.width<1, as.integer(bin.width*length(genes)), bin.width)
  stride <- ifelse(stride<1, as.integer(stride*length(genes)), stride)
  res_enrich <- list()
  pos_start <- 1
  while (pos_start < length(genes)) {
    pos_end <- pos_start + bin.width - 1
    if(pos_end > length(genes)) {
      pos_end <- length(genes)
      pos_start <- length(genes) - bin.width + 1
    }
    res_enrich[[paste0(pos_start, "-", pos_end)]] <- enrichGO(gene=genes[pos_start:pos_end], OrgDb=OrgDb, keyType="SYMBOL", ont=ont, universe=genes, ...)
    pos_start <- pos_start + stride
  }
  return(res_enrich)
}
