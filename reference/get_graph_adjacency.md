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
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
adj <- get_graph_adjacency(obj)
dim(adj)
#> [1] 100 100
adj[1:3, 1:3]
#>     C13 C28 C2
#> C13   0   1  0
#> C28   1   0  0
#> C2    0   0  0
```
