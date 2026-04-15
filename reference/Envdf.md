# Environmental variables (single block)

Environmental measurements per sample used together with the species/OTU
tables for correlation and Mantel-style examples.

## Usage

``` r
Envdf
```

## Format

A \`data.frame\` with 24 rows (samples) and 14 columns (environmental
variables).

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(Envdf)
head(Envdf)
#>       N    P     K    Ca    Mg    S    Al   Fe    Mn   Zn  Mo Baresoil Humdepth
#> 18 19.8 42.1 139.9 519.4  90.0 32.3  39.0 40.9  58.1  4.5 0.3     43.9      2.2
#> 15 13.4 39.1 167.3 356.7  70.7 35.2  88.1 39.0  52.4  5.4 0.3     23.6      2.2
#> 24 20.2 67.7 207.1 973.3 209.1 58.1 138.0 35.4  32.1 16.8 0.8     21.2      2.0
#> 27 20.6 60.8 233.7 834.0 127.2 40.7  15.4  4.4 132.0 10.7 0.2     18.7      2.9
#> 23 23.8 54.5 180.6 777.0 125.8 39.5  24.2  3.0  50.1  6.6 0.3     46.0      3.0
#> 19 22.8 40.9 171.4 691.8 151.4 40.8 104.8 17.6  43.6  9.1 0.4     40.5      3.8
#>     pH
#> 18 2.7
#> 15 2.8
#> 24 3.0
#> 27 2.8
#> 23 2.7
#> 19 2.7
```
