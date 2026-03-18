# SparCC p-value matrix (Rcpp implementation)

Fast C++-accelerated SparCC p-values via bootstrap and permutation. Uses
`sparcc_matrix_rcpp` internally for speed. Input: samples x taxa matrix.
Output: taxa x taxa p-value matrix.

## Usage

``` r
sparcc_pvalue_rcpp(
  data,
  iter = 20,
  inner_iter = 10,
  th = 0.1,
  R = 20,
  ncpus = 1
)
```

## Arguments

- data:

  Numeric matrix. Rows = samples, columns = taxa (ASVs/OTUs).

- iter:

  Integer. Number of outer iterations (Dirichlet resampling). Default
  20.

- inner_iter:

  Integer. Max inner iterations for pair exclusion. Default 10.

- th:

  Numeric. Threshold for excluding highly correlated pairs. Default 0.1.

- R:

  Integer. Number of bootstrap/permutation replicates. Default 20.

- ncpus:

  Integer. Number of CPUs for parallel boot. Default 1.

## Value

Numeric matrix of taxa x taxa p-values. Diagonal = 0. NaN indicates
observed correlation outside bootstrap CI.

## Examples

``` r
if (FALSE) { # \dontrun{
p_mat <- sparcc_pvalue_rcpp(asv_mat, R = 20, ncpus = 4)
} # }
```
