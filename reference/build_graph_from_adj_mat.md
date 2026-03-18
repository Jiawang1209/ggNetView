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
NULL
#> NULL
```
