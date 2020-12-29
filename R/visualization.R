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


#' @title ScatterPlotSTIP
#' @description  Generate scatter plot for gene in fitted matrix
#' @details Input fitted expression, output a scatter plot
#' @param data A fitted gene expression matrix
#' @param gene Plotted gene in the matrix
#' @return A ggplot object
#' @author Zhicheng Ji, Changxin Wan
#' @export ScatterPlotSTIP
#' @import ggplot2
ScatterPlotSTIP <- function(data, gene) {
  plot_df <- data.frame(t(data[gene, ]), row.names = colnames(data))
  plot_df$cell <- 1:nrow(plot_df)
  colnames(plot_df) <- c("gene", "cell")
  zp <- plot_df$cell[which(sapply(1:(nrow(plot_df)-1), function(x) plot_df[x, "gene"]*plot_df[x+1, "gene"]<=0))]

  p <- ggplot(plot_df, aes_string(x="cell", y="gene")) +
    geom_point(size=0.5) +
    labs(x="Cells", y=gene) +
    geom_vline(xintercept = zp, color="red", linetype="dashed", size=0.5) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", size=0.5) +
    theme(legend.position="none",
          legend.title = element_text(size = 7),
          legend.text = element_text(size=5),
          legend.key.width=unit(1, "lines"),
          plot.title = element_text(face="bold", hjust=0.5),
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(face="bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(face="bold"),
          axis.text.y = element_text(face="bold"),
          panel.background = element_blank(),
          panel.grid = element_blank())
  return(p)
}

