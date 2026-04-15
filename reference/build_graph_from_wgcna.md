# Build a graph object from WGCNA result

Build a graph object from WGCNA result

## Usage

``` r
build_graph_from_wgcna(
  wgcna_tom,
  module = NULL,
  node_annotation = NULL,
  directed = F,
  seed = 1115
)
```

## Arguments

- wgcna_tom:

  WGCNA TOM matrix.

- module:

  Module data frame

- node_annotation:

  Data Frame The annotation file of nodes in network

- directed:

  Logical (default: `FALSE`). Whether edges between nodes are directed.

- seed:

  nteger (default = 1115). Random seed for reproducibility.

## Value

An graph object representing the correlation network. Node/edge
attributes from WGCNA TOM matrix (optionally) module labels.

## Examples

``` r
if (FALSE) { # \dontrun{
# `wgcna_tom` is a long-format TOM edge list (from, to, weight)
# produced by `trans_TOM_in_WGCNA()`, and `module` is the WGCNA
# module-assignment data frame with columns (ID, Module).
obj <- build_graph_from_wgcna(
  wgcna_tom = wgcna_tom,
  module    = module
)
} # }
```
