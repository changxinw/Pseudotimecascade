# plotEnrichGroup

Visualize GO enrichment results for a specific expression pattern

## Usage

``` r
plotEnrichGroup(enrich_obj, terms = NULL, n = 10)
```

## Arguments

- enrich_obj:

  An `enrichResult` object (from
  [`enrichPattern()`](https://changxinw.github.io/Pseudotimecascade/reference/enrichPattern.md))
  or one element of a group enrichment list

- terms:

  Character vector of GO IDs to display. If NULL, top `n` terms are
  shown.

- n:

  Number of top terms to show if `terms = NULL` (default: 10)

## Value

A ggplot2 object of the bubble plot

## Details

This function takes an enrichment result from
[`enrichPattern()`](https://changxinw.github.io/Pseudotimecascade/reference/enrichPattern.md)
or group-enrichment list, selects user-defined GO terms (or top-ranked
ones), and generates a bubble plot summarizing the enrichment. Gene
ratio is shown on the x-axis, enriched terms on the y-axis, point size
represents gene count, and color represents adjusted q-value.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
