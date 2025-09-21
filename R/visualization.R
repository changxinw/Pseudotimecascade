#' @title HeatmapSTIP
#' @description  Generate heatmap for STIP result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export HeatmapSTIP
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @import dplyr
HeatmapSTIP <- function(x, gl, annotation, ...){
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- grDevices::colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = circlize::colorRamp2(myBreaks, myColor)
  gl <- intersect(rownames(x), gl)
  gene_labels = ComplexHeatmap::rowAnnotation(Genes = ComplexHeatmap::anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot = ComplexHeatmap::rowAnnotation(Genes = annotation[rownames(x)], col=list(Genes = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")))
  p = ComplexHeatmap::Heatmap(x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F, right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot, ...)
  return(p)
}


#' @title MSHeatmapSTIP
#' @description  Generate heatmap for multi-sample STIP result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param interval A list contains zero points for each gene in each sample
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export MSHeatmapSTIP
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark restore_matrix
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom grid grid.points grid.segments gpar unit
#' @importFrom stats sd qnorm na.omit
#' @import dplyr
MSHeatmapSTIP <- function(x, gl, annotation, interval, ...){
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- grDevices::colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = circlize::colorRamp2(myBreaks, myColor)
  gl <- intersect(rownames(x), gl)
  gene_labels <- ComplexHeatmap::rowAnnotation(Genes = ComplexHeatmap::anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot <- ComplexHeatmap::rowAnnotation(Genes = annotation[rownames(x)], col=list(Genes = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")))
  p <- ComplexHeatmap::Heatmap(x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = FALSE, show_row_names=FALSE, show_column_names = FALSE,
                               right_annotation = gene_labels, left_annotation = gene_annot,
                               layer_fun = function(j, i, x, y, w, h, fill){
                                 ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
                                 pindex <- c()
                                 sindex_s <- c()
                                 sindex_e <- c()
                                 for (rpn in seq(1, length(interval))){
                                   tmp_interval <- interval[[rpn]][names(annotation), ]
                                   for (m in 1:nrow(tmp_interval)){
                                     tmp_index <- stats::na.omit(unlist(tmp_interval[m, ], use.names=FALSE))
                                     if (length(tmp_index) != 0){
                                       pindex <- c(pindex, ind_mat[m, tmp_index])
                                       start <- max(1, round(mean(tmp_index) + (stats::qnorm(0.025)*stats::sd(tmp_index)/sqrt(length(tmp_index)))))
                                       sindex_s <- c(sindex_s, ind_mat[m, start])
                                       end <- round(mean(tmp_index) + (stats::qnorm(0.975)*stats::sd(tmp_index)/sqrt(length(tmp_index))))
                                       sindex_e <- c(sindex_e, ind_mat[m, end])
                                     }
                                   }
                                 }
                                 grid::grid.points(x[pindex], y[pindex], pch = 16, gp = grid::gpar(col = "black", alpha=0.5), size = grid::unit(2, "mm"))
                                 grid::grid.segments(x[sindex_s], y[sindex_s], x[sindex_e], y[sindex_e], gp = grid::gpar(col = "black", lwd = 0.5, alpha=0.5))}, ...)
  return(p)
}

#' #' @title ScatterPlotSTIP
#' #' @description  Generate scatter plot for gene in fitted matrix
#' #' @details Input fitted expression, output a scatter plot
#' #' @param data A fitted gene expression matrix
#' #' @param gene Plotted gene in the matrix
#' #' @param count set true if data is count matrix
#' #' @return A ggplot object
#' #' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' #' @export ScatterPlotSTIP
#' #' @import ggplot2
#' ScatterPlotSTIP <- function(data, gene, count=FALSE) {
#'   plot_df <- data.frame(t(data[gene, ]), row.names = colnames(data))
#'   plot_df$cell <- 1:nrow(plot_df)
#'   colnames(plot_df) <- c("gene", "cell")
#'   zp <- plot_df$cell[which(sapply(1:(nrow(plot_df)-1), function(x) plot_df[x, "gene"]*plot_df[x+1, "gene"]<=0))]
#'
#'   p <- ggplot(plot_df, aes_string(x="cell", y="gene")) +
#'     geom_point(size=0.5) +
#'     labs(x="Cells", y=gene) +
#'     theme(legend.position="none",
#'           legend.title = element_text(size = 7),
#'           legend.text = element_text(size=5),
#'           legend.key.width=unit(1, "lines"),
#'           plot.title = element_text(face="bold", hjust=0.5),
#'           plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"),
#'           axis.line = element_line(colour = "black"),
#'           axis.text.x = element_text(face="bold"),
#'           axis.title.x=element_blank(),
#'           axis.title.y=element_text(face="bold"),
#'           axis.text.y = element_text(face="bold"),
#'           panel.background = element_blank(),
#'           panel.grid = element_blank())
#'   if (!count){
#'     p <- p + geom_vline(xintercept = zp, color="red", linetype="dashed", size=0.5) +
#'              geom_hline(yintercept = 0, color="red", linetype="dashed", size=0.5)
#'   }
#'   return(p)
#' }
#'
#'
