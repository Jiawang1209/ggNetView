# Pre-computed adjacency matrix

A symmetric numeric adjacency matrix derived from \[otu_rare_relative\]
for tests that need a deterministic input.

## Usage

``` r
adjacency_matrix_example
```

## Format

A 2,859 x 2,859 numeric \`matrix\`.

## Source

Computed from \[otu_rare_relative\] using a Pearson correlation
threshold.

## Examples

``` r
data(adjacency_matrix_example)
adjacency_matrix_example[1:5, 1:5]
#>       ASV_1     ASV_2 ASV_3 ASV_4     ASV_6
#> ASV_1     0 0.0000000     0     0 0.0000000
#> ASV_2     0 0.0000000     0     0 0.8947427
#> ASV_3     0 0.0000000     0     0 0.0000000
#> ASV_4     0 0.0000000     0     0 0.0000000
#> ASV_6     0 0.8947427     0     0 0.0000000
```
