# Build a pie graph object from a data frame

Build a pie graph object from a data frame

## Usage

``` r
build_graph_from_pie(df, node_annotation = NULL, directed = F, seed = 1115)
```

## Arguments

- df:

  Data frame. Edge list with columns `from`, `to`, and optionally
  `weight`. If `weight` is absent, an unweighted graph is constructed.

- node_annotation:

  Data Frame The annotation file of nodes in network

- directed:

  Logical (default: `FALSE`). Whether edges between nodes are directed.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

An graph object representing the correlation network. Node/edge
attributes include correlation statistics and (optionally) module
labels.

## Examples

``` r
data(ppi_example)
obj <- build_graph_from_pie(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
obj
#> # A tbl_graph: 100 nodes and 50 edges
#> #
#> # An unrooted forest with 50 trees
#> #
#> # Node Data: 100 × 2 (active)
#>    name  group
#>    <chr> <chr>
#>  1 A1    A    
#>  2 A2    A    
#>  3 A3    A    
#>  4 A4    A    
#>  5 A5    A    
#>  6 A6    A    
#>  7 A7    A    
#>  8 A8    A    
#>  9 A9    A    
#> 10 A10   A    
#> # ℹ 90 more rows
#> #
#> # Edge Data: 50 × 3
#>    from    to weight
#>   <int> <int>  <dbl>
#> 1     1   100   45.2
#> 2     2    99   50.6
#> 3     3    98   37.8
#> # ℹ 47 more rows
```
