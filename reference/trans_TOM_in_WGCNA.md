# Translate TOM matrin to create graph object

Translate TOM matrin to create graph object

## Usage

``` r
trans_TOM_in_WGCNA(TOM, mat, threshold = NULL, top_k = NULL)
```

## Arguments

- TOM:

  matrix TOM matrix from WGCNA result

- mat:

  matrix matrox to WGCNA analysis

- threshold:

  numeric Absolute-weight cutoff. Edges with `abs(weight) <= threshold`
  are dropped. Set to `NULL` (default) to skip.

- top_k:

  integer Per-node top-`k` nearest-neighbour sparsification. For each
  node, only its `k` strongest TOM neighbours are eligible to contribute
  edges. An undirected edge (i, j) survives if j is in i's top-`k`
  **or** i is in j's top-`k` (mutual-kNN *union*, the standard choice
  for WGCNA visualization). Set to `NULL` (default) to skip. Typical
  values for WGCNA plots are 10 - 30.

## Value

A Data frame contain from, to and weight

## Examples

``` r
if (FALSE) { # \dontrun{
# `TOM` is a topological overlap matrix from WGCNA and `mat` is the
# expression matrix used to compute it.

# Simple absolute-weight threshold:
edge_df <- trans_TOM_in_WGCNA(TOM = TOM, mat = mat, threshold = 0.1)

# Per-node top-20 neighbours (recommended for visualization):
edge_df <- trans_TOM_in_WGCNA(TOM = TOM, mat = mat, top_k = 20)

# Combine: top-20 per node, then drop anything below 0.05 globally.
edge_df <- trans_TOM_in_WGCNA(TOM = TOM, mat = mat,
                              threshold = 0.05, top_k = 20)
head(edge_df)
} # }
```
