# Build a graph object from a data frame

Build a graph object from a data frame

## Usage

``` r
build_graph_from_df(
  df,
  node_annotation = NULL,
  directed = F,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- df:

  Data frame. Edge list with columns `from`, `to`, and optionally
  `weight`. If `weight` is absent, an unweighted graph is constructed.

- node_annotation:

  Data Frame The annotation file of nodes in network

- directed:

  Logical (default: `FALSE`). Whether edges between nodes are directed.

- module.method:

  Character Module analysis methods contains "Fast_greedy", "Walktrap",
  "Edge_betweenness", "Spinglass"

- top_modules:

  Integer. Number of top-ranked modules to select.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

An graph object representing the correlation network. Node/edge
attributes include correlation statistics and (optionally) module
labels.

## Examples

``` r
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation,
  directed        = FALSE,
  module.method   = "Fast_greedy",
  top_modules     = 5
)
obj
#> # A tbl_graph: 100 nodes and 50 edges
#> #
#> # An unrooted forest with 50 trees
#> #
#> # Node Data: 100 × 9 (active)
#>    name  group modularity modularity2 modularity3 Modularity Degree Segree
#>    <chr> <chr> <fct>      <fct>       <chr>       <fct>       <dbl>  <dbl>
#>  1 C13   C     1          1           1           1               1      1
#>  2 C28   C     1          1           1           1               1      1
#>  3 C2    C     10         10          10          10              1      1
#>  4 D9    D     10         10          10          10              1      1
#>  5 A3    A     11         11          11          11              1      1
#>  6 D38   D     11         11          11          11              1      1
#>  7 B12   B     12         12          12          12              1      1
#>  8 D19   D     12         12          12          12              1      1
#>  9 A1    A     13         13          13          13              1      1
#> 10 D40   D     13         13          13          13              1      1
#> # ℹ 90 more rows
#> # ℹ 1 more variable: Strength <dbl>
#> #
#> # Edge Data: 50 × 4
#>    from    to weight correlation
#>   <int> <int>  <dbl>       <dbl>
#> 1     9    10   45.2        45.2
#> 2    11   100   50.6        50.6
#> 3     5     6   37.8        37.8
#> # ℹ 47 more rows
```
