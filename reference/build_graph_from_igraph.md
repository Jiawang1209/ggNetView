# Build a graph object from a graph object

Build a graph object from a graph object

## Usage

``` r
build_graph_from_igraph(
  igraph,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- igraph:

  a igraph object

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- node_annotation:

  Data Frame The annotation file of nodes in network

- top_modules:

  Integer. Number of top-ranked modules to select.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

An graph object

## Examples

``` r
NULL
#> NULL
```
