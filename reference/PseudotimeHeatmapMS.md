# PseudotimeHeatmapMS

Generate heatmap for multi-sample Pseudotimecascade result

## Usage

``` r
PseudotimeHeatmapMS(x, gl, annotation, interval, ...)
```

## Arguments

- x:

  A gene expression matrix

- gl:

  Marked gene list

- annotation:

  Annotation matrix for genes in expression matrix

- interval:

  A list of switch point index data.frames (one per sample), or a single
  data.frame. Each element should have rows as genes and columns as
  samples, containing switch point indices (with `NA` allowed).

- ...:

  parameters passed to Heatmap

## Value

A ComplexHeatmap object

## Details

Input a gene expression matrix and annotation matrixes, output a
heatmap. Sample-level switch points are overlaid as points, and the
corresponding confidence intervals across samples are shown as line
segments.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
