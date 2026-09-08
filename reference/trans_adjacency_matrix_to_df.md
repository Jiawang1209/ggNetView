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

A data frame with columns \`from\`, \`to\`, \`weight\` describing the
edges of the adjacency matrix (one row per edge of the undirected
graph).

## Examples

``` r
data(adjacency_matrix_example)
set.seed(1)
idx <- sample(ncol(adjacency_matrix_example), 50)
edge_df <- trans_adjacency_matrix_to_df(adjacency_matrix_example[idx, idx])
head(edge_df)
#> [1] from   to     weight
#> <0 rows> (or 0-length row.names)
```
