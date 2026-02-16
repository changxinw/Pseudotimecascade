# PseudotimeHeatmap

Generate heatmap for Pseudotimecascade result

## Usage

``` r
PseudotimeHeatmap(x, gl, annotation, show_sp = FALSE, switch_point = NULL, ...)
```

## Arguments

- x:

  A gene expression matrix

- gl:

  Marked gene list

- annotation:

  Annotation matrix for genes in expression matrix

- show_sp:

  Logical; whether to overlay switch points on the heatmap (default:
  FALSE).

- switch_point:

  Named character vector of switch points for each gene (e.g.,
  `"12,35"`). Names should match `rownames(x)`; values are
  comma-separated column indices.

- ...:

  parameters passed to Heatmap

## Value

A ComplexHeatmap object

## Details

Input a gene expression matrix and annotation matrixes, output a
heatmap. Optionally, switch points can be overlaid on the heatmap as
dots when `show_sp = TRUE`.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
