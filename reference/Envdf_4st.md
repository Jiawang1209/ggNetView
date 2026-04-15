# Environmental variables across four sampling stages

Wide-form environmental data covering four stages (\`Env01\`-\`Env04\`),
used by \[gglink_heatmaps()\] examples.

## Usage

``` r
Envdf_4st
```

## Format

A \`data.frame\` with 24 rows (samples) and 56 columns.

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(Envdf_4st)
Envdf_4st[1:3, 1:6]
#>     N_1  P_1   K_1  Ca_1  Mg_1  S_1
#> 18 19.8 42.1 139.9 519.4  90.0 32.3
#> 15 13.4 39.1 167.3 356.7  70.7 35.2
#> 24 20.2 67.7 207.1 973.3 209.1 58.1
```
