# Two-block node table for double-matrix examples

Node-level information used to build cross-block ("double matrix")
networks. Each row corresponds to one node and identifies which of the
two blocks it belongs to.

## Usage

``` r
double_mat_node_df
```

## Format

A \`data.frame\` with 100 rows and 2 columns.

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(double_mat_node_df)
head(double_mat_node_df)
#>    name      type
#> 1 BASV1 Bacterial
#> 2 BASV2 Bacterial
#> 3 BASV3 Bacterial
#> 4 BASV4 Bacterial
#> 5 BASV5 Bacterial
#> 6 BASV6 Bacterial
```
