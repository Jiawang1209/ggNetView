# Two-block node table with pre-computed modularity

Same as \[double_mat_node_df\] with an extra \`Module\` / modularity
column attached so that examples can skip community detection.

## Usage

``` r
double_mat_node_df_with_modularity
```

## Format

A \`data.frame\` with 100 rows and 3 columns.

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(double_mat_node_df_with_modularity)
head(double_mat_node_df_with_modularity)
#>    name      type Modularity
#> 1 BASV1 Bacterial  Bacterial
#> 2 BASV2 Bacterial  Bacterial
#> 3 BASV3 Bacterial  Bacterial
#> 4 BASV4 Bacterial  Bacterial
#> 5 BASV5 Bacterial  Bacterial
#> 6 BASV6 Bacterial  Bacterial
```
