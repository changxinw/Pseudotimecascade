#' @title HeatmapSTIP
#' @description  Generate heatmap for STIP result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan
#' @export
#' @import ComplexHeatmap circlize dplyr grDevices

HeatmapSTIP <- function(x, gl, annotation, ...){
  paletteLength <- 1000
  myColor <- colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = colorRamp2(myBreaks, myColor)

  gl <- intersect(rownames(x), gl)
  gene_labels = rowAnnotation(Genes = anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot = rowAnnotation(Genes = annotation[rownames(x)], col=list(Genes = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")))

  p = Heatmap(x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F, right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot, ...)
  return(p)
}

