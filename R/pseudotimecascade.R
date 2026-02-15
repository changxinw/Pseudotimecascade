#' #' @title PreprocessPseudotime
#' #' @description Pseudotimecascade Preprocess
#' #' @details This function generates a table that performs Pseudotimecascade
#' #' @param data a single cell expression matrix or data.frame with rows as genes and columns as cells. Cells should follow the order of pseudo time
#' #' @param gl marked gene list
#' #' @return A Heatmap-class object
#' #' @export PreprocessPseudotime
#' PreprocessPseudotime <- function(data, gl){
#'   gene_group <- genePattern(data)
#'   plotdata <- data[rownames(gene_group), ]
#'   p <- PseudotimeHeatmap(plotdata, gl, as.matrix(gene_group)[, "pattern"])
#'   return(p)
#' }

#' @title genePattern
#' @description State transition pattern of each gene with flexible switch point definitions
#' @details  This function generates the state transition pattern of input gene
#' along pseudotime. Switch points can be defined by:
#' \itemize{
#'   \item \code{"zero_crossing"}: sign change around 0 (default).
#'   \item \code{"threshold_crossing"}: sign change around a user-defined \code{cutoff}.
#'   \item \code{"max_first_derivative"}: position with the largest absolute first derivative.
#' }
#' @param data a single cell expression matrix or data.frame with rows as genes and columns as cells.
#' Cells should follow the order of pseudo time.
#' @param switch_method character specifying switch point definition.
#' Default is \code{"zero_crossing"}.
#' @param cutoff numeric baseline used only when \code{switch_method = "threshold_crossing"}.
#' @return A dataframe of state transition pattern
#' @export genePattern
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
genePattern <- function(data,
                         switch_method = c("zero_crossing",
                                             "threshold_crossing",
                                             "max_first_derivative"),
                           cutoff = 0) {

  switch_method <- match.arg(switch_method)


  ## compute switch points list zp
  get_zp <- function(data) {

    # for zero_crossing
    if (switch_method == "zero_crossing") {
      zp <- apply(data, 1, function(sf){
        names(which(sapply(1:(length(sf)-1), function(i) sf[i]*sf[i+1] < 0)))
      })
      return(zp)
    }

    # for threshold_crossing
    if (switch_method == "threshold_crossing") {
      zp <- apply(data, 1, function(sf){
        sf2 <- sf - cutoff
        names(which(sapply(1:(length(sf2)-1), function(i) sf2[i]*sf2[i+1] < 0)))
      })
      return(zp)
    }


    if (switch_method == "max_first_derivative") {
      # Note: apply will simplify to a matrix/array if lengths differ; here length is always 1.
      zp <- apply(data, 1, function(sf){
        d <- diff(sf)
        if (all(is.na(d))) return(character(0))
        j <- which.max(abs(d))   # segment index
        paste0("V", j)
      })
      return(zp)
    }

    stop("Unsupported switch_method.")
  }

  zp <- get_zp(data)
  zp_direction <- sapply(names(zp), function(gene){
    direction <- lapply(zp[[gene]], function(y){
      ifelse(
        data[gene, as.numeric(sub("V", "", y))] <
          data[gene, as.numeric(sub("V", "", y)) + 1],
        "I", "D"
      )
    })
    return(paste0(unlist(direction), collapse=""))
  })

  gene_group <- data.frame(zp_direction, row.names = rownames(data))
  gene_group$rank_point <- sapply(rownames(gene_group), function(x)
    as.numeric(sub("V", "", zp[[x]][1]))
  )
  gene_group$switch_point <- sapply(rownames(gene_group), function(x)
    paste0(sub("V", "", zp[[x]]), collapse=",")
  )
  gene_group$switch_point_number <- sapply(zp, length)
  colnames(gene_group) <- c("pattern", "rank_point", "switch_point", "switch_point_number")
  gene_group <- gene_group[with(gene_group, order(pattern, rank_point)), ]

  return(gene_group)
}
