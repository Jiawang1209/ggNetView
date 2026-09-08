# Mantel test utilities for species-environment distance matrix correlation

These functions compute Mantel statistics between species and
environmental distance matrices. The implementation uses
[`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
and
[`vegan::mantel.partial`](https://vegandevs.github.io/vegan/reference/mantel.html)
directly, with a workflow designed for ggNetView's heatmap-link
visualization.

For each species column and each environmental column, builds a distance
matrix and runs Mantel test. Output format matches
[`psych::corr.test`](https://rdrr.io/pkg/psych/man/corr.test.html) for
drop-in use in `gglink_heatmaps`.

Runs Mantel test between each (spec_block, env_block) pair. Each block
is a subset of columns. Uses full distance matrices per block. Output
format is compatible with downstream processing when blocks are treated
as single "species" and "env" units.

Treats the whole `spec_df` as a single community matrix: all of its
columns together form ONE distance matrix
([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
with `spec_dist_method`). For each column of `env_df`, that single
column is converted into its own distance matrix
([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
with `env_dist_method`) and a Mantel test is run between the two
distance matrices.

## Usage

``` r
mantel_pairwise(
  spec_df,
  env_df,
  method = c("pearson", "spearman", "kendall"),
  alternative = c("two.sided", "less", "greater"),
  permutations = 999L,
  na_omit = TRUE,
  seed = NULL
)

mantel_between_blocks(
  spec,
  env,
  spec_select = NULL,
  env_select = NULL,
  test_type = c("mantel", "mantel.partial"),
  env_ctrl = NULL,
  method = c("pearson", "spearman", "kendall"),
  spec_dist_method = "euclidean",
  env_dist_method = "euclidean",
  na_omit = TRUE,
  permutations = 999L,
  seed = NULL
)

mantel_block_vs_col(
  spec_df,
  env_df,
  block_name = "block",
  method = c("pearson", "spearman", "kendall"),
  spec_dist_method = "bray",
  env_dist_method = "euclidean",
  permutations = 999L,
  na_omit = TRUE,
  seed = NULL
)
```

## Arguments

- spec_df:

  Data frame or matrix; rows = samples, columns = species (or any
  community variables). The full matrix is converted into ONE distance
  matrix.

- env_df:

  Data frame or matrix; rows = samples, columns = env variables. Each
  column is converted into its own distance matrix and tested
  separately.

- method:

  Correlation method for the Mantel test. One of `"pearson"`,
  `"spearman"`, or `"kendall"`.

- alternative:

  Alternative hypothesis for the test.

- permutations:

  Integer. Number of permutations for the test.

- na_omit:

  If `TRUE`, drop incomplete cases jointly across `spec_df` and `env_df`
  before computing distances.

- seed:

  Random seed for reproducibility.

- spec:

  Data frame of species abundances.

- env:

  Data frame of environmental variables.

- spec_select:

  Named list of column indices or names for species blocks. E.g.
  `list(block1 = 1:5, block2 = 6:10)`.

- env_select:

  Named list of column indices or names for env blocks.

- test_type:

  `"mantel"` or `"mantel.partial"`.

- env_ctrl:

  For `test_type = "mantel.partial"`, a data frame of controlling
  variables (same rows as spec/env).

- spec_dist_method:

  Distance method for the spec matrix
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)).
  Common ecological choices: `"bray"`, `"jaccard"`, `"euclidean"`.

- env_dist_method:

  Distance method for each env column
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)).
  Default `"euclidean"` is the standard choice for continuous env
  variables.

- block_name:

  Character (default `"block"`). Value placed in the `ID` column of the
  result, useful for tagging which block these rows came from when
  binding many results together.

## Value

A data frame with columns `ID` (species/block), `Type` (env/block),
`Correlation` (Mantel r), and `Pvalue`.

A data frame with one row per env column. Columns: `ID` (=
`block_name`), `Type` (env column name), `Correlation` (Mantel r),
`Pvalue` (Mantel p).

## Details

**NOTE (statistical caveat).** This is the **column-vs-column** Mantel
variant: each species column and each env column is reduced to a
single-variable distance matrix before the Mantel test. With one
variable per side,
[`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
is mathematically close to a (rank) correlation between the two columns
and does **not** carry the "community-vs-environment" interpretation
that ecology papers usually associate with a Mantel test. For the
standard ecological pattern, where a whole species block (community
matrix) is tested against each environmental gradient, use
`mantel_block_vs_col` instead.

This is the **ecologically meaningful** Mantel pattern (also used by
linkET / ggcor): "community structure of a spec block vs each
environmental gradient". Use this instead of `mantel_pairwise` when you
want a Mantel result that carries the standard
"community-vs-environment" interpretation.

## References

Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English
Edition. Elsevier.
