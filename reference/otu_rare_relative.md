# Relative-abundance OTU table

A relative-abundance transformation of \[otu_rare\]: each column sums
approximately to 1.

## Usage

``` r
otu_rare_relative
```

## Format

A \`data.frame\` with 2,859 rows and 18 columns.

## Source

Derived from \[otu_rare\].

## Examples

``` r
data(otu_rare_relative)
otu_rare_relative[1:5, 1:5]
#>              KO1        KO2        KO3        KO4        KO5
#> ASV_1 0.03306667 0.05453333 0.02013333 0.03613333 0.02686667
#> ASV_2 0.05750000 0.03393333 0.06046667 0.05810000 0.07320000
#> ASV_3 0.01733333 0.01296667 0.02290000 0.02336667 0.03106667
#> ASV_4 0.04266667 0.01093333 0.01416667 0.01933333 0.03346667
#> ASV_6 0.02646667 0.01856667 0.02110000 0.02353333 0.03806667
```
