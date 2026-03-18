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
NULL
#> NULL
```
