# genePattern

State transition pattern of each gene with flexible switch point
definitions

## Usage

``` r
genePattern(
  data,
  switch_method = c("zero_crossing", "threshold_crossing", "max_first_derivative"),
  cutoff = 0
)
```

## Arguments

- data:

  a single cell expression matrix or data.frame with rows as genes and
  columns as cells. Cells should follow the order of pseudo time.

- switch_method:

  character specifying switch point definition. Default is
  `"zero_crossing"`.

- cutoff:

  numeric baseline used only when
  `switch_method = "threshold_crossing"`.

## Value

A dataframe of state transition pattern

## Details

\#' @title PreprocessPseudotime \#' @description Pseudotimecascade
Preprocess \#' @details This function generates a table that performs
Pseudotimecascade \#' @param data a single cell expression matrix or
data.frame with rows as genes and columns as cells. Cells should follow
the order of pseudo time \#' @param gl marked gene list \#' @return A
Heatmap-class object \#' @export PreprocessPseudotime
PreprocessPseudotime \<- function(data, gl) gene_group \<-
genePattern(data) plotdata \<- data\[rownames(gene_group), \] p \<-
PseudotimeHeatmap(plotdata, gl, as.matrix(gene_group)\[, "pattern"\])
return(p)

This function generates the state transition pattern of input gene along
pseudotime. Switch points can be defined by:

- `"zero_crossing"`: sign change around 0 (default).

- `"threshold_crossing"`: sign change around a user-defined `cutoff`.

- `"max_first_derivative"`: position with the largest absolute first
  derivative.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
