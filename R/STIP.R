#' @title PreprocessSTIP
#' @description State Transition Inference Prediction Preprocess
#' @details This function generates a table that performs (STIP) State Transition Inference Prediction
#' @param x a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' @return A scaled gene expression matrix
#' @export
#' @import
#' @examples
#' PreprocessSTIP(fit_df)
PreprocessSTIP <- function(x, gl){
  colnames(x) <- 1: ncol(x)
  genes <- row.names(x)
  dn <- dimnames(x)
  x <- t(apply(x, 1, scale))
  dimnames(x) <- dn

  zpdirection <- x[, 1] < x[, ncol(x)]
  ### find the number of switch point
  zp <- apply(x, 1, function(sf) {
    names(which(sapply(1: (length(sf)-1), function(i) sf[i]*sf[i+1] < 0)))
  })
  zpnum <- sapply(zp, length)
  inczp <- names(which(zpdirection[zpnum==1]))
  deczp <- names(which(!zpdirection[zpnum==1]))
  multipoint <- names(zpnum)[zpnum > 1]

  gene_group = do.call("rbind", list(data.frame(a=inczp, b="inczp"), data.frame(a=deczp, b="deczp"), data.frame(a=multipoint, b="multipoint")))
  # gene_annot = gene_annot[, 1]
  gene_annot = gene_group[, 2]
  names(gene_annot) = gene_group[, 1]


  ### get each gene ranked
  geneorder <- NULL
  if (length(deczp) > 0) {
    geneorder <- c(geneorder, names(sort(unlist(zp[deczp]), decreasing = F)))
  }
  if (length(multipoint) > 0) {
    geneorder <- c(geneorder, names(sort(sapply(zp[multipoint], function(i) i[1]))))
  }
  if (length(inczp) > 0) {
    geneorder <- c(geneorder, names(sort(unlist(zp[inczp]))))
  }

  geneorder <- rev(geneorder)
  plotdata <- x[geneorder, ]

  # plotdata <- melt(plotdata)
  # colnames(plotdata) <- c("Gene", "Pseudotime", "Expression")
  # return(gene_annot)

  p = HeatmapSTIP(plotdata, gl, gene_annot)
  return(p)
}


#' @title HeatmapSTIP
#' @description  Generate heatmap for STIP result
#' @details Input a gene expression matrix and annotation matrixes, output a heatmap
#' @param x A gene expression matrix
#' @param genes Marked gene list
#' @param annotation Annotation matrix for genes in expression matrix
#' @return A ComplexHeatmap object
#' @export
#' @import ComplexHeatmap

HeatmapSTIP <- function(x, gl, annotation){
  paletteLength <- 1000
  myColor <- colorRampPalette(c("darkblue", "#6baed6", "#bdd7e7", "white", "#fcae91", "#fb6a4a", "darkred"))(paletteLength)
  myBreaks <- c(seq(min(x), 0, length.out=paletteLength/2), seq(max(x)/paletteLength, max(x), length.out=paletteLength/2))
  col = colorRamp2(myBreaks, myColor)

  gene_labels = rowAnnotation(Genes = anno_mark(at=which(rownames(x) %in% gl), labels=gl))
  gene_annot = rowAnnotation(Genes = annotation[rownames(x)], col=list(Genes = c(inczp="#e41a1c", deczp="#377eb8", multipoint="#4daf4a")))

  p = Heatmap(x, col=col, name="Expression", cluster_rows = FALSE, cluster_columns = F, right_annotation = gene_labels, show_row_names=F, show_column_names = F, left_annotation = gene_annot, use_raster=FALSE)
  return(p)
}
