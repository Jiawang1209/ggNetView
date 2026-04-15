# Build a graph object from Two numeric matrix

Build a graph object from Two numeric matrix

## Usage

``` r
build_graph_from_double_mat(
  mat1,
  mat2,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  directed = F,
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

A graph object representing the correlation-based two numeric matrix.

## Examples

``` r
# \donttest{
set.seed(1)
mat1 <- matrix(stats::rnorm(15 * 20), nrow = 15, ncol = 20)
mat2 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
rownames(mat1) <- paste0("A", seq_len(15))
rownames(mat2) <- paste0("B", seq_len(10))
colnames(mat1) <- colnames(mat2) <- paste0("sample", seq_len(20))
obj <- build_graph_from_double_mat(
  mat1          = mat1,
  mat2          = mat2,
  module.method = "Fast_greedy"
)
#> The max module in network is 6 we use the 6  modules for next analysis
obj
#> # A tbl_graph: 25 nodes and 150 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 25 × 8 (active)
#>    name  modularity modularity2 modularity3 Modularity Degree Segree Strength
#>    <chr> <fct>      <fct>       <chr>       <fct>       <dbl>  <dbl>    <dbl>
#>  1 B5    1          1           1           1              15     15     2.66
#>  2 B10   1          1           1           1              15     15     3.14
#>  3 A1    1          1           1           1              10     10     2.38
#>  4 A8    1          1           1           1              10     10     1.85
#>  5 B2    2          2           2           2              15     15     3.14
#>  6 B7    2          2           2           2              15     15     2.01
#>  7 A2    2          2           2           2              10     10     1.96
#>  8 A11   2          2           2           2              10     10     2.05
#>  9 A12   2          2           2           2              10     10     1.61
#> 10 A13   2          2           2           2              10     10     1.30
#> # ℹ 15 more rows
#> #
#> # Edge Data: 150 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1     3    21  0.118       0.118 Positive      
#> 2     3     5  0.101      -0.101 Negative      
#> 3     3    24  0.148       0.148 Positive      
#> # ℹ 147 more rows
# }
```
