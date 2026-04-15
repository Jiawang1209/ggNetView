# Protein-protein interaction example with module assignment

Same shape as \[ppi_example\] but the annotation table additionally
carries a \`Module\` column.

## Usage

``` r
ppi_module
```

## Format

A list with two elements (\`ppi\`, \`annotation\`).

## Source

Subset of the STRING database (https://string-db.org/).

## Examples

``` r
data(ppi_module)
names(ppi_module)
#> [1] "ppi"        "annotation"
```
