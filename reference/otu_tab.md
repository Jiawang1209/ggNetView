# Example OTU abundance table

A microbial Operational Taxonomic Unit (OTU) abundance matrix used by
the package examples. Rows are OTUs and columns are samples.

## Usage

``` r
otu_tab
```

## Format

A \`data.frame\` with 2,859 rows (OTUs) and 18 columns (samples). Cell
values are raw read counts.

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(otu_tab)
otu_tab[1:5, 1:5]
#>        KO1  KO2  KO3  KO4  KO5
#> ASV_1 1113 1968  816 1372 1062
#> ASV_2 1922 1227 2355 2218 2885
#> ASV_3  568  460  899  902 1226
#> ASV_4 1433  400  535  759 1287
#> ASV_6  882  673  819  888 1475
```
