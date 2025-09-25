#' @title plotEnrichGroup
#' @description Visualize GO enrichment results for a specific expression pattern
#' @details This function takes an enrichment result from \code{enrichPattern()}
#'          or group-enrichment list, selects user-defined GO terms (or top-ranked ones),
#'          and generates a bubble plot summarizing the enrichment.
#'
#' @param enrich_obj An \code{enrichResult} object (from \code{enrichPattern()})
#'                   or one element of a group enrichment list
#' @param terms Character vector of GO IDs to display. If NULL, top \code{n} terms are shown.
#' @param n Number of top terms to show if \code{terms = NULL} (default: 10)
#' @return A ggplot2 object of the bubble plot
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export plotEnrichGroup
#' @import ggplot2
#' @importFrom dplyr arrange
plotEnrichGroup <- function(enrich_obj, terms = NULL, n = 10) {
  df <- enrich_obj@result

  # Select the terms or take the first n
  if (!is.null(terms)) {
    df <- df[df$ID %in% terms, ]
  } else {
    df <- df %>% dplyr::arrange(qvalue) %>% head(n)
  }

  # Convert GeneRatio to a numerical value
  df$GeneRatio_num <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
  df <- df[order(df$qvalue), ]

  # Plot
  p <- ggplot(df, aes(GeneRatio_num, reorder(Description, GeneRatio_num))) +
    geom_point(aes(size = Count, color = qvalue)) +
    scale_color_gradient(low = "red", high = "blue", name = "q-value") +
    scale_size(range = c(2, 8), name = "Count") +
    labs(x = "Gene Ratio", y = NULL) +
    theme_bw()

  return(p)
}

#' @title plotEnrichBin
#' @description Visualize bin-based GO enrichment results
#' @details Given the output of \code{compareEnrichBin()}, this function generates
#'          a bubble plot of enriched GO terms across pseudotime bins. The top
#'          terms within each bin are selected based on q-value cutoff and ranking.
#'
#' @param bin_enrich A \code{compareClusterResult} object, typically from \code{compareEnrichBin()}.
#' @param n Number of top GO terms to select per bin (default = 5).
#' @param qval_cutoff q-value cutoff for filtering enriched terms (default = 0.05).
#' @param font.size Font size for the plot theme (default = 12).
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @return A ggplot2 object representing the bubble plot
#' @export plotEnrichBin
#' @import ggplot2 dplyr DOSE
plotEnrichBin <- function(bin_enrich, n = 5, qval_cutoff = 0.05, font.size = 12) {

  # select top enriched terms per bin
  tmp_enrich <- bin_enrich@compareClusterResult %>%
    dplyr::group_by(Cluster) %>%
    dplyr::filter(qvalue <= qval_cutoff) %>%
    dplyr::slice_min(order_by = qvalue, n = n) %>%
    dplyr::ungroup()

  # update results for plotting
  bin_enrich@compareClusterResult <- bin_enrich@compareClusterResult[
    bin_enrich@compareClusterResult$ID %in% unique(tmp_enrich$ID), ]
  bin_enrich@compareClusterResult$Description <- factor(
    bin_enrich@compareClusterResult$Description,
    levels = rev(unique(tmp_enrich$Description))
  )
  bin_enrich@compareClusterResult[bin_enrich@compareClusterResult$qvalue > 2*qval_cutoff, "qvalue"] <- 2*qval_cutoff

  # bubble plot
  p <- ggplot(bin_enrich@compareClusterResult, aes(x = Cluster, y = Description, size = Count)) +
    geom_point(aes(color = qvalue)) +
    scale_color_gradient2(low = "red", high = "blue", midpoint = qval_cutoff, name = "q-value") +
    DOSE::theme_dose(font.size = font.size) +
    labs(x = "Gene ranking", y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  return(p)
}


#' @title PseudotimeHeatmap
#' @description  Generate heatmap for Pseudotimecascade result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export PseudotimeHeatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @import dplyr
PseudotimeHeatmap <- function(x, gl, annotation, ...){
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- grDevices::colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = circlize::colorRamp2(myBreaks, myColor)
  gl <- intersect(rownames(x), gl)
  gene_labels = ComplexHeatmap::rowAnnotation(Pattern = ComplexHeatmap::anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot = ComplexHeatmap::rowAnnotation(pattern = annotation[rownames(x)], col=list(pattern = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")))
  p = ComplexHeatmap::Heatmap(x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F, right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot, ...)
  return(p)
}


#' @title PseudotimeHeatmapMS
#' @description  Generate heatmap for multi-sample Pseudotimecascade result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param interval A list contains zero points for each gene in each sample
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export PseudotimeHeatmapMS
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark restore_matrix
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom grid grid.points grid.segments gpar unit
#' @importFrom stats sd qnorm na.omit
#' @import dplyr
PseudotimeHeatmapMS <- function(x, gl, annotation, interval, ...){
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- grDevices::colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = circlize::colorRamp2(myBreaks, myColor)
  gl <- intersect(rownames(x), gl)
  gene_labels <- ComplexHeatmap::rowAnnotation(Pattern = ComplexHeatmap::anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot <- ComplexHeatmap::rowAnnotation(pattern = annotation[rownames(x)], col=list(pattern = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")))
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
                                 grid::grid.points(x[pindex], y[pindex], pch = 16, gp = grid::gpar(col = "black", alpha=0.2), size = grid::unit(0.5, "mm"))
                                 grid::grid.segments(x[sindex_s], y[sindex_s], x[sindex_e], y[sindex_e], gp = grid::gpar(col = "black", lwd = 0.5, alpha=0.5))}, ...)
  return(p)
}

#' #' @title ScatterPlotPseudotime
#' #' @description  Generate scatter plot for gene in fitted matrix
#' #' @details Input fitted expression, output a scatter plot
#' #' @param data A fitted gene expression matrix
#' #' @param gene Plotted gene in the matrix
#' #' @param count set true if data is count matrix
#' #' @return A ggplot object
#' #' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' #' @export ScatterPlotPseudotime
#' #' @import ggplot2
#' ScatterPlotPseudotime <- function(data, gene, count=FALSE) {
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
