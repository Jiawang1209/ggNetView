# Reorder the Modularity from graph object

Reorder the Modularity from graph object

## Usage

``` r
order_graph(graph_obj, order)
```

## Arguments

- graph_obj:

  An graph object from build_graph_from_mat or build_graph_from_df. The
  network object to be visualized.

- order:

  a character vectors

## Value

An graph object with reorder

## Examples

``` r
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
current <- levels(get_graph_nodes(obj)$Modularity)
obj2 <- order_graph(obj, order = rev(current))
levels(get_graph_nodes(obj2)$Modularity)
#>  [1] "Others" "22"     "21"     "20"     "2"      "19"     "18"     "17"    
#>  [9] "16"     "15"     "14"     "13"     "12"     "11"     "10"     "1"     
```
