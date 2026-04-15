# Rarefied OTU abundance table

A rarefied version of \[otu_tab\] in which every sample has been
subsampled to a common sequencing depth.

## Usage

``` r
otu_rare
```

## Format

A \`data.frame\` with 2,859 rows (OTUs) and 18 columns (samples).

## Source

Derived from \[otu_tab\] by rarefaction.

## Examples

``` r
data(otu_rare)
otu_rare[1:5, 1:5]
#>        KO1  KO2  KO3  KO4  KO5
#> ASV_1  992 1636  604 1084  806
#> ASV_2 1725 1018 1814 1743 2196
#> ASV_3  520  389  687  701  932
#> ASV_4 1280  328  425  580 1004
#> ASV_6  794  557  633  706 1142
```
