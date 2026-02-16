#' @title plotEnrichGroup
#' @description Visualize GO enrichment results for a specific expression pattern
#' @details This function takes an enrichment result from \code{enrichPattern()}
#'          or group-enrichment list, selects user-defined GO terms (or top-ranked ones),
#'          and generates a bubble plot summarizing the enrichment. Gene ratio is shown on the x-axis, enriched terms on the y-axis,
#'          point size represents gene count, and color represents adjusted q-value.
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
  p <- ggplot(df, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
    geom_point(aes(color = qvalue, size = Count)) +
    scale_color_gradient2(
      low = "red",
      high = "blue",
      mid = "purple",
      midpoint = 0.015,
      name = "qvalue"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "Count",
      guide = guide_legend(
        override.aes = list(color = "black", shape = 16),
        order = 1
      )
    ) +
    labs(x = "GeneRatio", y = "") +
    theme_bw() +
    theme(
      text = element_text(family = "Helvetica"),
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    ) +
    guides(
      color = guide_colorbar(order = 2)
    )
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
  tmp_enrich <- data.frame(
    bin_enrich@compareClusterResult %>%
      dplyr::group_by(Cluster) %>%
      dplyr::filter(qvalue <= qval_cutoff) %>%
      dplyr::arrange(qvalue) %>%
      dplyr::top_n(n = n, wt = -qvalue)
  )
  tmp_enrich <- tmp_enrich[order(tmp_enrich$Cluster, tmp_enrich$qvalue), ]

  # update results for plotting
  bin_enrich@compareClusterResult <- bin_enrich@compareClusterResult[
    bin_enrich@compareClusterResult$ID %in% unique(tmp_enrich$ID), ]
  bin_enrich@compareClusterResult$Description <- factor(
    bin_enrich@compareClusterResult$Description,
    levels = rev(unique(tmp_enrich$Description))
  )
  bin_enrich@compareClusterResult[bin_enrich@compareClusterResult$qvalue > 2*qval_cutoff, "qvalue"] <- 2*qval_cutoff

  # bubble plot
  p <- ggplot(
    bin_enrich@compareClusterResult,
    aes_string(x = "Cluster", y = "Description", size = "Count")
  ) +
    geom_point(aes_string(color = "qvalue")) +
    scale_color_gradient2(
      low = "red",
      high = "blue",
      midpoint = qval_cutoff,
      breaks = c(0, qval_cutoff),
      labels = c("0", as.character(qval_cutoff)),
      limits = c(0, max(bin_enrich@compareClusterResult$qvalue, na.rm = TRUE))
    ) +
    guides(
      size  = guide_legend(order = 1),
      color = guide_colorbar(order = 2, frame.colour = NA)
    ) +
    DOSE::theme_dose(font.size = font.size) +
    labs(y = "", x = "Gene ranking") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  return(p)
}


#' @title PseudotimeHeatmap
#' @description  Generate heatmap for Pseudotimecascade result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap.
#'          Optionally, switch points can be overlaid on the heatmap as dots when \code{show_sp = TRUE}.
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param show_sp Logical; whether to overlay switch points on the heatmap (default: FALSE).
#' @param switch_point Named character vector of switch points for each gene (e.g., \code{"12,35"}).
#'                    Names should match \code{rownames(x)}; values are comma-separated column indices.
#' @param ... parameters passed to Heatmap
#' @return A ComplexHeatmap object
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export PseudotimeHeatmap
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @import dplyr
PseudotimeHeatmap <- function(x, gl, annotation, show_sp = FALSE, switch_point = NULL, ...){
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- grDevices::colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = circlize::colorRamp2(myBreaks, myColor)
  font_family <- "Helvetica"
  mark_fontsize <- 12
  legend_fontsize <- 9
  legend_title_size <- 10

  gl <- intersect(rownames(x), gl)
  gene_labels = ComplexHeatmap::rowAnnotation(
    Pattern = ComplexHeatmap::anno_mark(
      at = which(rownames(x) %in% gl),
      labels = gl,
      labels_gp = grid::gpar(fontsize = mark_fontsize, fontfamily = font_family),
      link_gp   = grid::gpar(lwd = 0.6),
      padding   = grid::unit(3, "mm"),
      extend    = grid::unit(4, "mm")
    ),
    annotation_name_gp = grid::gpar(fontsize = legend_title_size, fontfamily = font_family)
  )

  gene_annot = ComplexHeatmap::rowAnnotation(
    Pattern = annotation[rownames(x)],
    col = list(Pattern = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")),
    annotation_name_gp = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
    annotation_legend_param = list(
      title_gp  = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
      labels_gp = grid::gpar(fontsize = legend_fontsize,  fontfamily = font_family)
    )
  )

  ## precompute switch point positions
  sp_flag <- NULL
  if (isTRUE(show_sp) && !is.null(switch_point)) {

    ## align by gene names: switch_point must be a named vector
    sp_vec <- switch_point[rownames(x)]
    sp_vec[is.na(sp_vec)] <- ""

    sp_flag <- matrix(FALSE, nrow = nrow(x), ncol = ncol(x))

    for (i in seq_along(sp_vec)) {
      s <- sp_vec[[i]]
      if (s == "") next
      idx <- suppressWarnings(as.integer(strsplit(s, ",", fixed = TRUE)[[1]]))
      idx <- idx[!is.na(idx)]
      idx <- idx[idx >= 1 & idx <= ncol(x)]
      if (length(idx) == 0) next
      sp_flag[i, idx] <- TRUE
    }
  }

  if (isTRUE(show_sp) && !is.null(sp_flag)) {
    p = ComplexHeatmap::Heatmap(
      x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F,
      right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot,
      heatmap_legend_param = list(
        title_gp  = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
        labels_gp = grid::gpar(fontsize = legend_fontsize,  fontfamily = font_family)
      ),
      layer_fun = function(j, i, x, y, w, h, fill) {
        hit <- sp_flag[cbind(i, j)]
        if (any(hit)) {
          grid::grid.points(
            x[hit], y[hit],
            pch = 16,
            size = grid::unit(0.7, "mm"),
            gp = grid::gpar(col = "grey60")
          )
        }
      },
      ...
    )
  } else {
    p = ComplexHeatmap::Heatmap(
      x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F,
      right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot,
      heatmap_legend_param = list(
        title_gp  = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
        labels_gp = grid::gpar(fontsize = legend_fontsize,  fontfamily = font_family)
      ),
      ...
    )
  }

  return(p)
}


#' @title PseudotimeHeatmapMS
#' @description  Generate heatmap for multi-sample Pseudotimecascade result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap.
#'          Sample-level switch points are overlaid as points, and the corresponding confidence intervals
#'          across samples are shown as line segments.
#' @param x A gene expression matrix
#' @param gl Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @param interval A list of switch point index data.frames (one per sample),
#'                or a single data.frame. Each element should have rows as genes and columns as samples,
#'                containing switch point indices (with \code{NA} allowed).
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
  #   interval <- interval[names(annotation), ]
  if (is.data.frame(interval) | is.matrix(interval)) {
    interval <- list(interval)
  }
  x <- x[names(annotation), ]
  paletteLength <- 1000
  myColor <- colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = colorRamp2(myBreaks, myColor)
  font_family <- "Helvetica"
  mark_fontsize <- 12
  legend_fontsize <- 9
  legend_title_size <- 10

  gl <- intersect(rownames(x), gl)
  gene_labels <- rowAnnotation(
    Pattern = anno_mark(
      at = which(rownames(x) %in% gl),
      labels = gl,
      labels_gp = grid::gpar(fontsize = mark_fontsize, fontfamily = font_family),
      link_gp   = grid::gpar(lwd = 0.6),
      padding   = grid::unit(3, "mm"),
      extend    = grid::unit(2, "mm")
    ),
    annotation_name_gp = grid::gpar(fontsize = legend_title_size, fontfamily = font_family)
  )

  gene_annot <- rowAnnotation(
    Pattern = annotation[rownames(x)],
    col = list(Pattern = c(I="#d73027", ID="#f46d43", IDI="#fdae61", IDID="#fee090", D="#4575b4", DI="#74add1", DID="#abd9e9", DIDI="#e0f3f8")),
    annotation_name_gp = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
    annotation_legend_param = list(
      title_gp  = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
      labels_gp = grid::gpar(fontsize = legend_fontsize,  fontfamily = font_family)
    )
  )

  p <- Heatmap(
    x, col=col, name="Expression",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names=FALSE, show_column_names = FALSE,
    right_annotation = gene_labels, left_annotation = gene_annot,
    heatmap_legend_param = list(
      title_gp  = grid::gpar(fontsize = legend_title_size, fontfamily = font_family),
      labels_gp = grid::gpar(fontsize = legend_fontsize,  fontfamily = font_family)
    ),
    layer_fun = function(j, i, x, y, w, h, fill){
      ind_mat <- restore_matrix(j, i, x, y)
      pindex <- c()
      sindex_s <- c()
      sindex_e <- c()

      for (rpn in seq(1, length(interval))){
        tmp_interval <- interval[[rpn]][names(annotation), ]
        for (m in 1:nrow(tmp_interval)){
          tmp_index <- na.omit(unlist(tmp_interval[m, ], use.names=FALSE))
          if (length(tmp_index) != 0){
            pindex <- c(pindex, ind_mat[m, tmp_index])
            start <- round(mean(tmp_index) + (qnorm(0.025)*sd(tmp_index)/sqrt(length(tmp_index))))
            sindex_s <- c(sindex_s, ind_mat[m, start])
            end <- round(mean(tmp_index) + (qnorm(0.975)*sd(tmp_index)/sqrt(length(tmp_index))))
            sindex_e <- c(sindex_e, ind_mat[m, end])
          }
        }
      }
      grid.points(x[pindex], y[pindex], pch = 16, gp = gpar(col = "black", alpha=0.5), size = unit(2, "mm"))
      grid.segments(x[sindex_s], y[sindex_s], x[sindex_e], y[sindex_e], gp = gpar(col = "black", lwd = 0.5, alpha=0.5))
    },
    ...
  )
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
