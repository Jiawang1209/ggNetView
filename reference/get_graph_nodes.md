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
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
head(get_graph_nodes(obj))
#>   name group modularity modularity2 modularity3 Modularity Degree Segree
#> 1  C13     C          1           1           1          1      1      1
#> 2  C28     C          1           1           1          1      1      1
#> 3   C2     C         10          10          10         10      1      1
#> 4   D9     D         10          10          10         10      1      1
#> 5   A3     A         11          11          11         11      1      1
#> 6  D38     D         11          11          11         11      1      1
#>   Strength
#> 1 26.73367
#> 2 26.73367
#> 3 37.37186
#> 4 37.37186
#> 5 37.80905
#> 6 37.80905
```
