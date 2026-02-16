# plotEnrichBin

Visualize bin-based GO enrichment results

## Usage

``` r
plotEnrichBin(bin_enrich, n = 5, qval_cutoff = 0.05, font.size = 12)
```

## Arguments

- bin_enrich:

  A `compareClusterResult` object, typically from
  [`compareEnrichBin()`](https://changxinw.github.io/Pseudotimecascade/reference/compareEnrichBin.md).

- n:

  Number of top GO terms to select per bin (default = 5).

- qval_cutoff:

  q-value cutoff for filtering enriched terms (default = 0.05).

- font.size:

  Font size for the plot theme (default = 12).

## Value

A ggplot2 object representing the bubble plot

## Details

Given the output of
[`compareEnrichBin()`](https://changxinw.github.io/Pseudotimecascade/reference/compareEnrichBin.md),
this function generates a bubble plot of enriched GO terms across
pseudotime bins. The top terms within each bin are selected based on
q-value cutoff and ranking.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
