# SelectTopGenes

Select top-ranked genes from multi-sample Pseudotimecascade results

## Usage

``` r
SelectTopGenes(gene_mean_list, top = 1000)
```

## Arguments

- gene_mean_list:

  A multi-sample result list returned by
  [`PseudotimeMSProcess()`](https://changxinw.github.io/Pseudotimecascade/reference/PseudotimeMSProcess.md).

- top:

  Number of top-ranked genes to select (default = 1000).

## Value

A list with the same structure as `gene_mean_list`, but restricted to
the selected top genes.

## Details

This function ranks genes based on the average significance across
samples. Specifically, q-values are transformed using
\\-\log\_{10}(q)\\, and the mean transformed value across samples is
used for ranking. The top `top` genes are selected and all corresponding
components (mean expression, pattern assignment, q-values, and switch
points) are subset accordingly.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
