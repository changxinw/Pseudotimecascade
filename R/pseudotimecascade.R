#' #' @title PreprocessSTIP
#' #' @description State Transition Inference Prediction Preprocess
#' #' @details This function generates a table that performs (STIP) State Transition Inference Prediction
#' #' @param data a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' #' @param gl marked gene list
#' #' @return A Heatmap-class object
#' #' @export PreprocessSTIP
#' PreprocessSTIP <- function(data, gl){
#'   gene_group <- genePattern(data)
#'   plotdata <- data[rownames(gene_group), ]
#'   p <- HeatmapSTIP(plotdata, gl, as.matrix(gene_group)[, "pattern"])
#'   return(p)
#' }

#' @title genePattern
#' @description State transition pattern of each gene
#' @details This function generates the state transition pattern of input gene
#' @param data a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' @return A dataframe of state transtion pattern
#' @export genePattern
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
genePattern <- function(data){
  ### separation based on zero points
  ### find the number of switch point, only specify the name of mono
  zp <- apply(data, 1, function(sf){
    names(which(sapply(1: (length(sf)-1), function(i) sf[i]*sf[i+1] < 0)))
  })

  zp_direction <- sapply(names(zp), function(gene){
    direction <- lapply(zp[[gene]], function(y){
      ### fix the bug when switch point is 1
      ifelse(data[gene, as.numeric(sub("V", "", y))] < data[gene, as.numeric(sub("V", "", y)) + 1], "I", "D")
    })
    return(paste0(unlist(direction), collapse=""))
  })

  gene_group <- data.frame(zp_direction, row.names = rownames(data))
  gene_group$rank_point <- sapply(rownames(gene_group), function(x) as.numeric(sub("V", "", zp[[x]][1])))
  gene_group$switch_point <- sapply(rownames(gene_group), function(x) paste0(sub("V", "", zp[[x]]), collapse=","))
  gene_group$switch_point_number <- sapply(zp, length)
  colnames(gene_group) <- c("pattern", "rank_point", "switch_point", "switch_point_number")
  gene_group <- gene_group[with(gene_group, order(pattern, rank_point)), ]

  return(gene_group)
}
