# Trans adjacency matrix to edge data frame

Trans adjacency matrix to edge data frame

## Usage

``` r
trans_adjacency_matrix_to_df(adjacency_matrix)
```

## Arguments

- adjacency_matrix:

  Numeric matrix. A numeric matrix with adjacency matrix.

## Value

a data frame of edge

## Examples

``` r
data(adjacency_matrix_example)
set.seed(1)
idx <- sample(ncol(adjacency_matrix_example), 50)
edge_df <- trans_adjacency_matrix_to_df(adjacency_matrix_example[idx, idx])
head(edge_df)
#> [1] from to  
#> <0 rows> (or 0-length row.names)
```
