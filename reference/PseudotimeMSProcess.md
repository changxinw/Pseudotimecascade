# PseudotimeMSProcess

Perform multi-sample pattern integration for Pseudotimecascade results

## Usage

``` r
PseudotimeMSProcess(integrate_list, conf = 0.75, qval = 0.05)
```

## Arguments

- integrate_list:

  A named list of single-sample Pseudotimecascade results. Each element
  must contain `gene_group` and `fit_data_list`.

- conf:

  Minimum confidence threshold for consensus pattern selection (default
  = 0.75).

- qval:

  q-value cutoff for filtering significant genes (default = 0.05).

## Value

A list (hereafter referred to as `gene_mean_list`) containing:

- mean_expr:

  Matrix of averaged fitted expression across samples.

- mean_pattern:

  Gene pattern result based on mean trajectory.

- df_pattern:

  Consensus pattern table with confidence scores.

- df_qvalue:

  Filtered q-value matrix aligned to consensus genes.

- df_switch_point:

  List of switch point interval matrices across samples.

## Details

This function integrates single-sample Pseudotimecascade results stored
in a list and performs consensus pattern analysis across samples. For
each gene, the most frequent expression pattern is identified, and a
confidence score is computed as the proportion of samples supporting the
dominant pattern. Genes with confidence greater than or equal to `conf`
and passing the q-value threshold are retained.

For retained genes, the function:

- Computes the average fitted expression trajectory across valid
  samples.

- Reassigns gene patterns based on the mean trajectory.

- Collects sample-level switch points and organizes them into interval
  matrices.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
