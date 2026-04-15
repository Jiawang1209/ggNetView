# Extract subgrah from graph object

Extract subgrah from graph object

## Usage

``` r
get_subgraph(graph_obj, select_module = NULL)
```

## Arguments

- graph_obj:

  An graph object from build_graph_from_mat or build_graph_from_df. The
  network object to be visualized.

- select_module:

  a character vectors Select the module name in graph object

## Value

list

## Examples

``` r
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
sg <- get_subgraph(obj, select_module = "1")
#>    Module Number
#> 1       1      2
#> 2      10      2
#> 3      11      2
#> 4      12      2
#> 5      13      2
#> 6      14      2
#> 7      15      2
#> 8      16      2
#> 9      17      2
#> 10     18      2
#> 11     19      2
#> 12      2      2
#> 13     20      2
#> 14     21      2
#> 15     22      2
#> 16 Others     70
names(sg)
#> [1] "sub_graph_all"    "stat_module"      "sub_graph_select"
sg$stat_module
#>    Module Number
#> 1       1      2
#> 2      10      2
#> 3      11      2
#> 4      12      2
#> 5      13      2
#> 6      14      2
#> 7      15      2
#> 8      16      2
#> 9      17      2
#> 10     18      2
#> 11     19      2
#> 12      2      2
#> 13     20      2
#> 14     21      2
#> 15     22      2
#> 16 Others     70
```
