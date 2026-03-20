# Build a graph object from a graph object

Build a graph object from a graph object

## Usage

``` r
build_graph_from_igraph(
  igraph,
  use_existing_modules = TRUE,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  module_attr = NULL,
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- igraph:

  a igraph object

- use_existing_modules:

  Logical (default = TRUE). If `TRUE`, reuse module information already
  stored on igraph vertices when available.

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- module_attr:

  Character or NULL (default = NULL). Optional vertex attribute name
  containing precomputed module labels. When `NULL`, the function will
  auto-detect one of `"Modularity"`, `"modularity2"`, `"modularity3"`,
  or `"modularity"`.

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
