# Extract adjacency matrix from graph object

Returns the adjacency matrix of the graph with node names as rownames
and colnames. Non-zero entries indicate edges; values correspond to edge
weights when the graph is weighted.

## Usage

``` r
get_graph_adjacency(graph_obj)
```

## Arguments

- graph_obj:

  A graph object from `build_graph_from_mat` or `build_graph_from_df`.

## Value

A numeric matrix with rownames and colnames set to node IDs (`name`).

## Examples

``` r
NULL
#> NULL
```
