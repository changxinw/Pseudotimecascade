# compareEnrichBin

Bin-based enrichment for ordered genes in a specific pattern

## Usage

``` r
compareEnrichBin(
  gene.group,
  pattern,
  bin.width = 0.2,
  stride = 0.1,
  species = "human",
  ont = "BP",
  universe = FALSE,
  ...
)
```

## Arguments

- gene.group:

  A data frame indicating genes in each pattern

- pattern:

  The expression pattern for enrichment analysis

- bin.width:

  The width of each bin

- stride:

  Stride of each step

- species:

  Select from human or mouse

- ont:

  One of "BP", "MF", and "CC" subontologies, or "ALL" for all three

- universe:

  Pass to the universe paramenter of enrichGO

- ...:

  Pass to the function enrichGO

## Value

A compareClusterResult object summarizing enrichment results across bins

## Details

This function partitions genes within a temporal expression pattern into
bins along pseudotime and performs GO enrichment on each bin. The output
contains enrichment results for all bins

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
