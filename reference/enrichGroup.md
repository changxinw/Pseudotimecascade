# enrichGroup

Enrichment for genes in specific pattern

## Usage

``` r
enrichGroup(gene.group, species = "mouse", ont = "BP", ...)
```

## Arguments

- gene.group:

  a data frame indicate genes in each pattern

- species:

  select from human or mouse

- ont:

  One of "BP", "MF", and "CC" subontologies, or "ALL" for all three

- ...:

  pass to the function enrichGO

## Value

A list of enrichResult instance

## Details

This function generates an enrichResult instance

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
