#' @title fitData
#' @description Fit single cell gene expression data according to cell pseudotime
#' @details This function generates a fitted expression matrix of single cell RNA-seq
#' @param data a single cell expression matrix with rows as genes and columns as cells.
#' @param pseudo.time cells ranked according to pseudo time
#' @param zero.rate cutoff of zero rate among cells to filter gene
#' @param p.adjust.method method for multiple hypothesis test
#' @param p.cutoff cutoff for p value
#' @param adj.p.cutoff cutoff for adjusted p value
#' @param res column number of output fitted matrix
#' @param verbose show message of running process
#' @param mc.cores number of cores for parallel computing
#' @return A scaled gene expression matrix
#' @author Zhicheng Ji, Changxin Wan
#' @export fitData
#' @import VGAM parallel
#' @importFrom stats p.adjust

fitData <- function(data, pseudo.time=colnames(data), zero.rate=0.9, p.adjust.method="BH", p.cutoff=0.05, adj.p.cutoff=0.05, res=ncol(data), verbose=TRUE, mc.cores=1){
  data <- data[, pseudo.time]
  data <- data[rowMeans(data == 0) <= zero.rate, ]
  ### pt should be a parameter to set
  pt <- 1:res
  ### output percentage of process
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
    scale(fitted(model[[sg]]))
  }))
  ### get p values for each gene
  pval_data <- lapply(names(model), function(sg) {
    lrtest(model[[sg]])@Body[2,5]
  })
  ### p value adjustment method
  qval_data <- p.adjust(unlist(pval_data), method = p.adjust.method)
  ### set cutoff for fdr here
  fit_data_sig <- as.data.frame(fit_data[qval_data<=adj.p.cutoff & unlist(pval_data)<=p.cutoff, ])

  return(fit_data_sig)
}

