# SPIEC-EASI matrix (Rcpp implementation)

Fast C++-accelerated SPIEC-EASI. Uses Rcpp for CLR transformation; Uses
`huge` for estimation and `huge.select` for StARS model selection.
Input: samples x taxa matrix.

## Usage

``` r
spieceasi_matrix_rcpp(
  data,
  method = "glasso",
  output = c("partial_correlation", "adjacency", "covariance", "correlation",
    "precision", "beta", "stability"),
  ...
)
```

## Arguments

- data:

  Numeric matrix. Rows = samples, columns = taxa (ASVs/OTUs).

- method:

  Character. "glasso" or "mb".

- output:

  Character. One of "partial_correlation", "adjacency", "covariance",
  "correlation", "precision", "beta", "stability".

- ...:

  Passed to `spieceasi_fit_rcpp` (e.g. `pulsar.select`, `pulsar.params`
  or `stars.params`, `verbose`).

## Value

Matrix as requested by `output`.

## Examples

``` r
if (FALSE) { # \dontrun{
pcor_mat <- spieceasi_matrix_rcpp(asv_mat, method = "glasso", output = "partial_correlation")
adj_mat  <- spieceasi_matrix_rcpp(asv_mat, method = "glasso", output = "adjacency")
} # }
```
