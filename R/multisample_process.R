#' @title PseudotimeMSProcess
#' @description Perform multi-sample pattern integration for Pseudotimecascade results
#'
#' @details
#' This function integrates single-sample Pseudotimecascade results stored in a list and
#' performs consensus pattern analysis across samples. For each gene, the
#' most frequent expression pattern is identified, and a confidence score
#' is computed as the proportion of samples supporting the dominant pattern.
#' Genes with confidence greater than or equal to \code{conf} and passing
#' the q-value threshold are retained.
#'
#' For retained genes, the function:
#' \itemize{
#'   \item Computes the average fitted expression trajectory across valid samples.
#'   \item Reassigns gene patterns based on the mean trajectory.
#'   \item Collects sample-level switch points and organizes them into interval matrices.
#' }
#'
#' @param integrate_list A named list of single-sample Pseudotimecascade results. Each element
#' must contain \code{gene_group} and \code{fit_data_list}.
#' @param conf Minimum confidence threshold for consensus pattern selection (default = 0.75).
#' @param qval q-value cutoff for filtering significant genes (default = 0.05).
#'
#' @return A list (hereafter referred to as \code{gene_mean_list}) containing:
#' \describe{
#'   \item{mean_expr}{Matrix of averaged fitted expression across samples.}
#'   \item{mean_pattern}{Gene pattern result based on mean trajectory.}
#'   \item{df_pattern}{Consensus pattern table with confidence scores.}
#'   \item{df_qvalue}{Filtered q-value matrix aligned to consensus genes.}
#'   \item{df_switch_point}{List of switch point interval matrices across samples.}
#' }
#'
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export PseudotimeMSProcess

PseudotimeMSProcess <- function(integrate_list, conf=0.75, qval=0.05){
  stip_genes <- unique(unlist(lapply(names(integrate_list), function(bm){rownames(integrate_list[[bm]][["gene_group"]])})))
  df_pattern <- data.frame(row.names=stip_genes)
  df_qvalue <- data.frame(row.names=stip_genes)

  for (pt in names(integrate_list)){
    gene_group <- integrate_list[[pt]][["gene_group"]][, "pattern"]
    names(gene_group) <- rownames(integrate_list[[pt]][["gene_group"]])
    df_pattern[, pt] <- gene_group[stip_genes]
    vec_qval <- integrate_list[[pt]][["fit_data_list"]][["qval"]]
    names(vec_qval) <- rownames(integrate_list[[pt]][["fit_data_list"]][["data"]])
    df_qvalue[, pt] <- vec_qval[stip_genes]
  }
  df_qvalue_mtx <- (df_qvalue <= qval)
  df_qvalue_mtx[!df_qvalue_mtx] <- NA
  df_pattern <- as.data.frame(lapply(colnames(df_pattern), function(x){df_pattern[df_qvalue_mtx[, x], x]}), row.names=rownames(df_pattern), col.names=colnames(df_pattern))
  df_pattern <- df_pattern[rowSums(is.na(df_pattern)) != ncol(df_pattern), ]
  pattern <- unlist(apply(df_pattern, 1, function(x){names(sort(table(x), decreasing = TRUE))[1]}))
  confidence <- unlist(apply(df_pattern, 1, function(x){sort(table(x), decreasing = TRUE)[1]})) / ncol(df_pattern)
  ### select genes with high confidence
  df_pattern$pattern <- pattern
  df_pattern$confidence <- confidence
  df_pattern <- df_pattern[df_pattern$confidence >= conf, ]

  df_pattern_mtx <- apply(df_pattern[, 1:(ncol(df_pattern)-2)], 2, function(x) {x==df_pattern[, "pattern"]})
  df_pattern_mtx[!df_pattern_mtx] <- NA

  ### update df_qvalue matrix to eliminate those genes with an inconsistant pattern
  df_qvalue <- df_qvalue[rownames(df_pattern), ] * df_pattern_mtx

  ### get gene expression from each valid sample
  gene_mean_list <- list()
  genes <- rownames(df_pattern_mtx)
  for (gene in genes){
    samples <- na.omit(colnames(df_pattern_mtx)[df_pattern_mtx[gene, ]])
    gene_mean_list[[gene]] <- rowMeans(data.frame(lapply(samples, function(x) {as.vector(integrate_list[[x]][['fit_data_list']][["data"]][gene, ])})), na.rm=TRUE)
  }
  gene_mean_expr <- do.call("rbind", gene_mean_list)
  gene_mean_pattern <- genePattern(as.data.frame(gene_mean_expr))

  ### get gene switch point list
  ### generate confidence interval dataframe list, q005
  maxsp <- max(gene_mean_pattern[, 4])
  genes <- rownames(df_pattern_mtx)
  interval_list <- list()

  for (i in seq(1, maxsp)){
    rank_point_list <- lapply(names(integrate_list), function(x) {unlist(lapply(strsplit(integrate_list[[x]][["gene_group"]][genes, "switch_point"], ","), function(y){as.numeric(y[i])}))})
    interval_list[[i]] <- as.data.frame(rank_point_list, row.names=genes, col.names=names(integrate_list))
    interval_list[[i]] <- interval_list[[i]] * df_pattern_mtx
  }

  ### output result
  gene_mean_list <- list()
  gene_mean_list[["mean_expr"]] <- gene_mean_expr
  gene_mean_list[["mean_pattern"]] <- gene_mean_pattern
  gene_mean_list[["df_pattern"]] <- df_pattern
  gene_mean_list[["df_qvalue"]] <- df_qvalue
  gene_mean_list[["df_switch_point"]] <- interval_list
  return(gene_mean_list)
}



#' @title SelectTopGenes
#' @description Select top-ranked genes from multi-sample Pseudotimecascade results
#'
#' @details
#' This function ranks genes based on the average significance across samples.
#' Specifically, q-values are transformed using \eqn{-\log_{10}(q)}, and the mean
#' transformed value across samples is used for ranking. The top \code{top}
#' genes are selected and all corresponding components (mean expression,
#' pattern assignment, q-values, and switch points) are subset accordingly.
#'
#' @param gene_mean_list A multi-sample result list returned by
#' \code{PseudotimeMSProcess()}.
#' @param top Number of top-ranked genes to select (default = 1000).
#'
#' @return A list with the same structure as \code{gene_mean_list},
#' but restricted to the selected top genes.
#'
#' @author Zhicheng Ji, Changxin Wan, Beijie Ji
#' @export SelectTopGenes
SelectTopGenes <- function(gene_mean_list, top=1000){
  df_qvalue <- gene_mean_list[["df_qvalue"]][rownames(gene_mean_list[["df_pattern"]]), ]
  df_qvalue[df_qvalue==0] <- 10
  df_qvalue[df_qvalue==10] <- min(df_qvalue, na.rm = TRUE)
  df_qvalue <- -log10(df_qvalue)
  df_qvalue$mean <- rowMeans(df_qvalue, na.rm=TRUE)

  top_genes <- rownames(df_qvalue[order(-df_qvalue$mean), ])[1:top]
  gene_mean_list[["df_qvalue"]][top_genes, ]

  ### output matrix for top genes
  top_mean_list <- list()
  top_mean_list[["mean_expr"]] <- gene_mean_list[["mean_expr"]][top_genes, ]
  top_mean_list[["mean_pattern"]] <- gene_mean_list[["mean_pattern"]][top_genes, ]
  top_mean_list[["mean_pattern"]] <- top_mean_list[["mean_pattern"]][order(top_mean_list[["mean_pattern"]]$pattern, top_mean_list[["mean_pattern"]]$rank_point), ]
  top_mean_list[["df_pattern"]] <- gene_mean_list[["df_pattern"]][top_genes, ]
  top_mean_list[["df_qvalue"]] <- gene_mean_list[["df_qvalue"]][top_genes, ]
  top_mean_list[["df_switch_point"]] <- list()
  for (i in seq(1, length(gene_mean_list[["df_switch_point"]]))){
    top_mean_list[["df_switch_point"]][[i]] <- gene_mean_list[["df_switch_point"]][[i]][top_genes, ]
  }
  return(top_mean_list)
}
