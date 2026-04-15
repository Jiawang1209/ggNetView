# Taxonomy table for \[otu_NatureWater\]

Taxonomic annotation matching \[otu_NatureWater\].

## Usage

``` r
tax_NatureWater
```

## Format

A \`tbl_df\` with 51,998 rows and 18 columns.

## Source

Adapted from the same published study as \[otu_NatureWater\].

## Examples

``` r
data(tax_NatureWater)
head(tax_NatureWater)
#> # A tibble: 6 × 18
#>   ID    Feature_ID     Kindom Phylum Class Order Family Genus Species Confidence
#>   <chr> <chr>          <chr>  <chr>  <chr> <chr> <chr>  <chr> <chr>        <dbl>
#> 1 ASV1  feb89dc9caf19… Bacte… Prote… Gamm… Xant… Rhoda… Dokd… Unknown      0.999
#> 2 ASV2  08bd45c0807b6… Bacte… Prote… Gamm… Burk… Rhodo… Zoog… Unknown      0.999
#> 3 ASV3  2a2af3995b876… Bacte… Bacte… Bact… Chit… Sapro… mida… midas_…      0.970
#> 4 ASV4  f759b2828f1c6… Bacte… Bacte… Bact… Chit… Sapro… OLB8  midas_…      1.000
#> 5 ASV5  adfb67f460d65… Bacte… Prote… Gamm… Burk… Rhodo… Unkn… Unknown      1.000
#> 6 ASV6  e38bafa6321bd… Bacte… Nitro… Nitr… Nitr… Nitro… Nitr… Nitros…      0.974
#> # ℹ 8 more variables: frequency <dbl>, niche_width <dbl>, Taxa <chr>,
#> #   Group <chr>, Partition <chr>, Function <chr>, Process <chr>, MRA <dbl>
```
