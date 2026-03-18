# Extract node table from graph object

Returns the node (vertex) table as a data frame, including attributes
such as `name`, `Modularity`, `Degree`, `Strength`, etc.

## Usage

``` r
get_graph_nodes(graph_obj)
```

## Arguments

- graph_obj:

  A graph object from `build_graph_from_mat` or `build_graph_from_df`.

## Value

A data frame with one row per node; columns include node attributes.

## Examples

``` r
NULL
#> NULL
```
