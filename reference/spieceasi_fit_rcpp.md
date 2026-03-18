# SPIEC-EASI fit (Rcpp-accelerated CLR)

Same as `spieceasi_fit` but uses fast C++ CLR transformation.

## Usage

``` r
spieceasi_fit_rcpp(
  data,
  method = "glasso",
  sel.criterion = "stars",
  verbose = TRUE,
  pulsar.select = TRUE,
  pulsar.params = list(),
  stars.params = list(),
  lambda.log = TRUE,
  ...
)
```
