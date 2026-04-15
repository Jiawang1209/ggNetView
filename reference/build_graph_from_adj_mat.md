# Build a correlation-based network from a adjacency matrix

Build a correlation-based network from a adjacency matrix

## Usage

``` r
build_graph_from_adj_mat(
  adjacency_matrix,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- adjacency_matrix:

  Numeric matrix. A numeric matrix with adjacency matrix.

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- node_annotation:

  Data frame. Optional node annotation table, containing metadata such
  as taxonomy or functional categories.

- top_modules:

  Integer. Number of top-ranked modules to retain for downstream
  visualization or analysis.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

A graph object representing the correlation-based microbial network.
Node/edge attributes include correlation statistics and (optionally)
module labels.

## Examples

``` r
data(adjacency_matrix_example)
set.seed(1)
idx <- sample(ncol(adjacency_matrix_example), 80)
adj <- adjacency_matrix_example[idx, idx]
obj <- build_graph_from_adj_mat(
  adjacency_matrix = adj,
  module.method    = "Fast_greedy",
  top_modules      = 5
)
#> The max module in network is 3 we use the 3  modules for next analysis
obj
#> # A tbl_graph: 6 nodes and 3 edges
#> #
#> # An unrooted forest with 3 trees
#> #
#> # Node Data: 6 × 7 (active)
#>   name     modularity modularity2 modularity3 Modularity Degree Strength
#>   <chr>    <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl>
#> 1 ASV_1029 1          1           1           1               1    0.808
#> 2 ASV_1655 1          1           1           1               1    0.808
#> 3 ASV_1848 2          2           2           2               1    0.809
#> 4 ASV_1120 2          2           2           2               1    0.809
#> 5 ASV_2230 3          3           3           3               1    0.857
#> 6 ASV_2437 3          3           3           3               1    0.857
#> #
#> # Edge Data: 3 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1     1     2  0.808       0.808 Positive      
#> 2     3     4  0.809       0.809 Positive      
#> 3     5     6  0.857       0.857 Positive      
```
