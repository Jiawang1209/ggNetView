# Build a graph object from multiple numeric matrices

Build a graph object from multiple numeric matrices

## Usage

``` r
build_graph_from_multi_mat(
  mat1,
  mat2,
  ...,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  directed = FALSE,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- mat1:

  Numeric matrix. A numeric matrix with samples in colums and variables
  in rows

- mat2:

  Numeric matrix. A numeric matrix with samples in colums and variables
  in rows

- ...:

  Additional numeric matrices with samples in columns and variables in
  rows.

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- node_annotation:

  Data frame. Optional node annotation table, containing metadata such
  as taxonomy or functional categories.

- directed:

  Logical (default: `FALSE`). Whether edges between nodes are directed.

- top_modules:

  Integer. Number of top-ranked modules to retain for downstream
  visualization or analysis.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

A graph object representing the correlation-based multi numeric matrix.

## Examples

``` r
# \donttest{
set.seed(1)
mat1 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
mat2 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
mat3 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
rownames(mat1) <- paste0("A", seq_len(10))
rownames(mat2) <- paste0("B", seq_len(10))
rownames(mat3) <- paste0("C", seq_len(10))
colnames(mat1) <- colnames(mat2) <- colnames(mat3) <-
  paste0("sample", seq_len(20))
obj <- build_graph_from_multi_mat(
  mat1          = mat1,
  mat2          = mat2,
  mat3,
  module.method = "Fast_greedy"
)
#> The max module in network is 5 we use the 5  modules for next analysis
obj
#> # A tbl_graph: 30 nodes and 300 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 30 × 8 (active)
#>    name  modularity modularity2 modularity3 Modularity Degree Segree Strength
#>    <chr> <fct>      <fct>       <chr>       <fct>       <dbl>  <dbl>    <dbl>
#>  1 A1    1          1           1           1              20     20     8.54
#>  2 A6    1          1           1           1              20     20     8.30
#>  3 A9    1          1           1           1              20     20     7.06
#>  4 B2    1          1           1           1              20     20     6.69
#>  5 B6    1          1           1           1              20     20     8.04
#>  6 B8    1          1           1           1              20     20     8.36
#>  7 C3    1          1           1           1              20     20     5.57
#>  8 C5    1          1           1           1              20     20     5.39
#>  9 C6    1          1           1           1              20     20     5.98
#> 10 C10   1          1           1           1              20     20     7.25
#> # ℹ 20 more rows
#> #
#> # Edge Data: 300 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1     1    18  0.693      -0.693 Negative      
#> 2     1     4  0.313      -0.313 Negative      
#> 3     1    29  0.904      -0.904 Negative      
#> # ℹ 297 more rows
# }
```
