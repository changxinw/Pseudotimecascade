# enrichPattern

Flexible GO enrichment for one, multiple, or all temporal patterns

## Usage

``` r
enrichPattern(
  gene.group,
  patterns = NULL,
  species = "mouse",
  ont = "BP",
  universe = NULL,
  ...
)
```

## Arguments

- gene.group:

  A data frame indicating genes and their assigned temporal patterns

- patterns:

  Character vector of patterns to analyze (e.g., c("I","D")), or NULL to
  analyze all patterns

- species:

  "mouse" or "human"

- ont:

  One of "BP", "MF", "CC", or "ALL"

- universe:

  Background set of genes (default: NULL → auto-detected)

- ...:

  Additional arguments passed to enrichment function

## Value

A named list of enrichResult objects, one per pattern

## Details

This function performs enrichment analysis for user-specified patterns.
If no pattern is provided, all patterns in gene.group will be analyzed

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
