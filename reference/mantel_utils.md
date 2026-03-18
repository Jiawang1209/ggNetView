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

## Usage

``` r
mantel_pairwise(
  spec_df,
  env_df,
  method = c("pearson", "spearman", "kendall"),
  alternative = c("two.sided", "less", "greater"),
  permutations = 999L,
  na_omit = TRUE
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
```

## Arguments

- spec_df:

  Data frame or matrix of species abundances (samples as rows).

- env_df:

  Data frame or matrix of environmental variables (samples as rows).

- method:

  Correlation method: `"pearson"`, `"spearman"`, or `"kendall"`.

- alternative:

  Alternative hypothesis for the test.

- permutations:

  Number of permutations.

- na_omit:

  If `TRUE`, remove incomplete cases.

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

  Distance method for species matrix when using `mantel_between_blocks`.
  One of `"euclidean"`, `"bray"`, `"manhattan"`, etc. (see
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)).

- env_dist_method:

  Distance method for environmental matrix.

- seed:

  Random seed for reproducibility.

## Value

A data frame with columns `ID` (species/block), `Type` (env/block),
`Correlation` (Mantel r), and `Pvalue`.

## References

Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English
Edition. Elsevier.
