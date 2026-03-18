# Build a graph object from a adjacency matrix and module Info

Build a graph object from a adjacency matrix and module Info

## Usage

``` r
build_graph_from_adj_mat_module(
  adjacency_matrix,
  node_annotation = NULL,
  directed = F,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- adjacency_matrix:

  Numeric matrix. A numeric matrix with adjacency matrix.

- node_annotation:

  Data Frame \#' The annotation file of nodes in network contain `node`
  and `Modularity`

- directed:

  Logical (default: `FALSE`). Whether edges between nodes are directed.

- top_modules:

  Integer. Number of top-ranked modules to select.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

An graph object representing the correlation network.

## Examples

``` r
NULL
#> NULL
```
