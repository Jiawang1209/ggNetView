# Translate TOM matrin to create graph object

Translate TOM matrin to create graph object

## Usage

``` r
trans_TOM_in_WGCNA(TOM, mat, threshold = NULL)
```

## Arguments

- TOM:

  matrix TOM matrix from WGCNA result

- mat:

  matrix matrox to WGCNA analysis

- threshold:

  numeric the threshold in weight

## Value

A Data frame contain from, to and weight

## Examples

``` r
if (FALSE) { # \dontrun{
# `TOM` is a topological overlap matrix from WGCNA and `mat` is the
# expression matrix used to compute it.
edge_df <- trans_TOM_in_WGCNA(TOM = TOM, mat = mat, threshold = 0.1)
head(edge_df)
} # }
```
