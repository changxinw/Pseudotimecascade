#' @title fitData
#' @description Fit single cell gene expression data according to cell pseudotime
#' @details This function generates a fitted expression matrix of single cell RNA-seq
#' @param data a single cell expression matrix with rows as genes and columns as cells.
#' @param expr.cut cutoff of lowerest expression
#' @param expr.cut.rate cutoff of cells lower than lowest expression
#' @param pseudo.time cells ranked according to pseudo time
#' @param p.adjust.method method for multiple hypothesis test
#' @param new_data input matrix for model prediction
#' @param verbose show message of running process
#' @param mc.cores number of cores for parallel computing
#' @return A list contains scaled fitted gene expression matrix
#' @author Zhicheng Ji, Changxin Wan
#' @export fitData
#' @import VGAM parallel
#' @importFrom stats p.adjust

fitData <- function(data, pt=1:ncol(data), expr.cut=0.1, expr.cut.rate=0.05, pseudo.time=colnames(data), p.adjust.method="BH", new_data=data.frame(pt=seq(1, ncol(data))), verbose=TRUE, mc.cores=1){
  data <- data[, pseudo.time]
  data <- data[rowMeans(data>expr.cut)>expr.cut.rate, ,drop = FALSE]
  ### output percentage of process
#   pt <- 1:ncol(data)
  if (mc.cores > 1) {
    model <- mclapply(1:nrow(data), function(igene){
      if (igene %% 100 == 0 && verbose){
        message(paste0(igene, " genes processed!"))
      }
      vgam(data[igene, ] ~ s(pt, df=2), family=tobit(Lower=0.1, type.fitted = "mean.obs"), control=vgam.control(maxit=50))
    }, mc.cores=mc.cores)
  }
  else {
    model <- lapply(1:nrow(data), function(igene){
      if (igene %% 100 == 0 && verbose){
        message(paste0(igene, " genes processed!"))
      }
      vgam(data[igene, ] ~ s(pt, df=2), family=tobit(Lower=0.1, type.fitted = "mean.obs"), control=vgam.control(maxit=50))
    })
  }
  names(model) <- row.names(data)

  ### get the scaled fitted value for each gene
  fit_data <- t(sapply(names(model),function(sg) {
#     scale(fitted(model[[sg]])) #
      scale(predict(model[[sg]], new_data, type="response"))
  }))
  ### get p values for each gene
  if (mc.cores > 1){
    pval_data <- mclapply(names(model), function(sg) {
      lrtest(model[[sg]])@Body[2,5]
    }, mc.cores = mc.cores)
  } else {
    pval_data <- lapply(names(model), function(sg) {
      lrtest(model[[sg]])@Body[2,5]
    })
  }
  ### p value adjustment method
  qval_data <- p.adjust(unlist(pval_data), method = p.adjust.method)
  ### set cutoff for fdr here
  # fit_data_sig <- as.data.frame(fit_data[qval_data<=adj.p.cutoff & unlist(pval_data)<=p.cutoff, ])

  res_list <- list()
  res_list[["model"]] <- model
  res_list[["data"]] <- fit_data
  res_list[["pval"]] <- unlist(pval_data)
  res_list[["qval"]] <- qval_data

  return(res_list)
}
