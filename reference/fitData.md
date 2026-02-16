# fitData

Fit single cell gene expression data along pseudotime with flexible
smoothing parameters

## Usage

``` r
fitData(
  data,
  expr.cut = 0.1,
  expr.cut.rate = 0.05,
  pseudo.time = colnames(data),
  p.adjust.method = "BH",
  pt = 1:ncol(data),
  new_data = data.frame(pt = seq(1, ncol(data))),
  verbose = TRUE,
  mc.cores = 1
)
```

## Arguments

- data:

  a single cell expression matrix with rows as genes and columns as
  cells.

- expr.cut:

  cutoff of lowest expression.

- expr.cut.rate:

  cutoff proportion of cells above expression threshold.

- pseudo.time:

  cells ranked according to pseudotime.

- p.adjust.method:

  method for multiple hypothesis test.

- pt:

  numeric pseudotime index used in model fitting.

- new_data:

  input matrix for model prediction.

- verbose:

  show message of running process.

- mc.cores:

  number of cores for parallel computing.

- df:

  degrees of freedom for spline smoothing (default = 2).

- tobit_lower:

  lower bound for tobit model (default = 0.1).

- maxit:

  maximum number of iterations for model fitting (default = 50).

## Value

A list containing fitted models, scaled fitted expression matrix, raw
p-values and adjusted q-values.

## Details

This function fits gene expression trajectories using a VGAM tobit model
with smoothing spline on pseudotime. Model parameters such as degrees of
freedom, tobit lower bound, and maximum iterations can be adjusted.

## Author

Zhicheng Ji, Changxin Wan, Beijie Ji
