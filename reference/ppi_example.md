# Protein-protein interaction example

A small protein-protein interaction (PPI) network used by the
\`build_graph_from_pie\` and STRING-DB examples. The list contains the
interaction edge table and a node annotation table.

## Usage

``` r
ppi_example
```

## Format

A list with two elements:

- \`ppi\`:

  Edge \`data.frame\` (\`from\`, \`to\`, \`weight\`).

- \`annotation\`:

  Node annotation \`data.frame\`.

## Source

Subset of the STRING database (https://string-db.org/).

## Examples

``` r
data(ppi_example)
names(ppi_example)
#> [1] "ppi"        "annotation"
```
