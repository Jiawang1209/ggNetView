# Extract node and edge information from graph object

Extract node and edge information from graph object

## Usage

``` r
get_info_from_graph(graph_obj)
```

## Arguments

- graph_obj:

  An graph object from build_graph_from_mat or build_graph_from_df. The
  network object to be visualized.

## Value

a list

## Examples

``` r
data(adjacency_matrix_example)
set.seed(1)
idx <- sample(ncol(adjacency_matrix_example), 80)
obj <- build_graph_from_adj_mat(adjacency_matrix_example[idx, idx])
#> The max module in network is 3 we use the 3  modules for next analysis
info <- get_info_from_graph(obj)
names(info)
#> [1] "node_info" "edge_info"
head(info$node_info)
#> # A tibble: 6 × 4
#>   name     Modularity Degree Strength
#>   <chr>    <ord>       <dbl>    <dbl>
#> 1 ASV_1029 1               1    0.808
#> 2 ASV_1655 1               1    0.808
#> 3 ASV_1848 2               1    0.809
#> 4 ASV_1120 2               1    0.809
#> 5 ASV_2230 3               1    0.857
#> 6 ASV_2437 3               1    0.857
head(info$edge_info)
#> # A tibble: 3 × 5
#>   from     to       weight correlation corr_direction
#>   <chr>    <chr>     <dbl>       <dbl> <chr>         
#> 1 ASV_1029 ASV_1655  0.808       0.808 Positive      
#> 2 ASV_1848 ASV_1120  0.809       0.809 Positive      
#> 3 ASV_2230 ASV_2437  0.857       0.857 Positive      
```
