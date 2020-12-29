#' @title PreprocessSTIP
#' @description State Transition Inference Prediction Preprocess
#' @details This function generates a table that performs (STIP) State Transition Inference Prediction
#' @param x a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' @param gl marked gene list
#' @return A scaled gene expression matrix
#' @export PreprocessSTIP
PreprocessSTIP <- function(x, gl){
  gene_group <- genePattern(x)
  plotdata <- x[rownames(gene_group), ]
  p <- HeatmapSTIP(plotdata, gl, as.matrix(gene_group)[, "pattern"])
  return(p)
}

#' @title genePattern
#' @description State transition pattern of each gene
#' @details This function generates the state transition pattern of input gene
#' @param x a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' @return A dataframe of state transtion pattern
#' @author Zhicheng Ji, Changxin Wan
genePattern <- function(x){
  ### separation based on zero points
  ### find the number of switch point, only specify the name of mono
  zp <- apply(x, 1, function(sf){
    names(which(sapply(1: (length(sf)-1), function(i) sf[i]*sf[i+1] < 0)))
  })

  zp_direction <- sapply(names(zp), function(gene){
    direction <- lapply(zp[[gene]], function(y){
      ifelse(x[gene, as.numeric(sub("V", "", y)) - 1] < 0, "I", "D")
    })
    return(paste0(unlist(direction), collapse=""))
  })

  gene_group <- data.frame(zp_direction, row.names = rownames(x))
  gene_group$zp <- sapply(rownames(gene_group), function(x) as.numeric(sub("V", "", zp[[x]][1])))
  gene_group$zpnum <- sapply(zp, length)
  colnames(gene_group) <- c("pattern", "zp", "zpnum")
  gene_group <- gene_group[with(gene_group, order(pattern, zp)), ]

  return(gene_group)
}
