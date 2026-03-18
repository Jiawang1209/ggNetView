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
NULL
#> NULL
```
