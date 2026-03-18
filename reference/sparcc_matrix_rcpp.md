# SparCC correlation matrix (Rcpp implementation)

Fast C++ implementation of SparCC for inferring correlations from
compositional count data. Input: samples x taxa matrix. Output: taxa x
taxa correlation matrix.

## Usage

``` r
sparcc_matrix_rcpp(data, iter = 20, inner_iter = 10, th = 0.1, nthreads = 0L)
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

- nthreads:

  Integer. Number of OpenMP threads for parallel `sparccinner` (0 = use
  default). Only effective when OpenMP is available.

## Value

Numeric matrix of taxa x taxa correlations (median across iterations).

## Examples

``` r
if (FALSE) { # \dontrun{
# data: samples x taxa count matrix
cor_mat <- sparcc_matrix_rcpp(asv_mat, iter = 20, inner_iter = 10, th = 0.1)
# With 4 threads (when OpenMP available)
cor_mat <- sparcc_matrix_rcpp(asv_mat, iter = 20, nthreads = 4)
} # }
```
