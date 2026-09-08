# Network + module-environment heatmaps with shared Mantel API

Render a `ggNetView` network in the centre, surrounded by up to four
environmental-correlation heatmap quadrants, with link segments
connecting each module's anchor to the corresponding env-variable points
on the diagonals. Each module is represented by a single per-sample
summary (eigengene or abundance) for downstream statistics. Supports
both Pearson/Spearman/Kendall correlation and
[`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
tests, and exposes the same Mantel API as
[`gglink_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md).

## Usage

``` r
ggnetview_modularity_heatmaps(
  graph_obj,
  env,
  otu_mat,
  env_select = NULL,
  module_index = c("eigengene", "abundance"),
  abundance_type = c("sum", "mean"),
  relation_method = c("correlation", "mantel"),
  cor.method = c("pearson", "kendall", "spearman"),
  cor.use = c("everything", "all", "complete", "pairwise", "na"),
  mantel.method2 = c("pearson", "kendall", "spearman"),
  mantel_kind = c("block_vs_col", "col_vs_col"),
  spec_dist_method = "bray",
  env_dist_method = "euclidean",
  permutations = 999L,
  mantel.seed = 1115,
  drop_nonsig = FALSE,
  layout = "gephi",
  layout_module = c("random", "adjacent", "order"),
  orientation = c("top_right", "bottom_right", "top_left", "bottom_left"),
  distance = 3,
  r = 6,
  HeatmapScale = 1,
  SigLineAlpha = 0.5,
  HeatmapLabelSize = 5,
  HeatmapSigSize = 5,
  HeatmapColorBar = NULL,
  HeatmapLabelOrient = 0,
  SigLineWidth = c(0.5, 2),
  SigLineColor = c("#fdbb84", "#d7301f"),
  HeatmapPointSize = 5,
  HeatmapPointFill = "#de77ae",
  HeatmapTileColor = NA,
  HeatmapTileSize = 0,
  ...,
  layout.module = deprecated()
)
```

## Arguments

- graph_obj:

  A `tbl_graph` (e.g. from `build_graph_from_mat` or
  `build_graph_from_df`). Must have a node `name` attribute and a module
  column (one of `"Modularity"`, `"modularity3"`, `"modularity2"`;
  auto-detected).

- env:

  Data frame or matrix of environmental variables. Rows are samples
  (rownames matched against `otu_mat` columns), columns are env factors.

- otu_mat:

  Numeric matrix. Rows = OTUs/ASVs (rownames matched against `graph_obj`
  node names), columns = samples (colnames matched against `env`
  rownames). Used to compute module eigengenes / abundances and, in
  block-vs-col Mantel mode, to assemble per-module community distance
  matrices.

- env_select:

  Named list (required). Column indices or names of `env` that form each
  env block, one block per heatmap quadrant. `length(env_select)` must
  equal `length(orientation)`. The list names (`names(env_select)`) are
  used by `comparisons_groups` and in the returned stats. Example:
  `list(Env01 = 1:5, Env02 = 6:10, Env03 = 11:15, Env04 = 16:20)`.

- module_index:

  Character (default `"eigengene"`). How each module is summarised into
  one per-sample value used downstream. `"eigengene"` = PC1 of the
  module's OTU sub-matrix (recommended); `"abundance"` = sum or mean of
  OTU abundances within the module (controlled by `abundance_type`).

- abundance_type:

  Character (default `"sum"`). Only used when
  `module_index = "abundance"`. Either `"sum"` or `"mean"`.

- relation_method:

  Character (default `"correlation"`). One of `"correlation"` or
  `"mantel"`.

- cor.method:

  Character (default `"pearson"`). Correlation method used by
  [`psych::corr.test`](https://rdrr.io/pkg/psych/man/corr.test.html) for
  env-env (heatmap tiles) and, when `relation_method = "correlation"`,
  for module-env links. One of `"pearson"`, `"kendall"`, `"spearman"`.

- cor.use:

  Character (default `"everything"`). Missing-value handling for
  [`psych::corr.test`](https://rdrr.io/pkg/psych/man/corr.test.html).
  One of `"everything"`, `"all"`, `"complete"`, `"pairwise"`, `"na"`.

- mantel.method2:

  Character (default `"pearson"`). Correlation coefficient passed to
  [`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
  as its `method` argument. One of `"pearson"`, `"kendall"`,
  `"spearman"`.

- mantel_kind:

  Character (default `"block_vs_col"`). Which Mantel algorithm to use;
  see **Details**. The same parameter is exposed in
  [`gglink_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md).

- spec_dist_method:

  Character (default `"bray"`). Dissimilarity method
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html))
  used to convert a module's OTU sub-matrix into ONE community distance
  matrix when `mantel_kind = "block_vs_col"`.

- env_dist_method:

  Character (default `"euclidean"`). Distance method
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html))
  used to convert each env column into its own distance matrix when
  `relation_method = "mantel"`.

- permutations:

  Integer (default `999L`). Number of permutations passed to
  [`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html).

- mantel.seed:

  Integer (default `1115`). Seed forwarded to the Mantel helpers so the
  permutation p-values are reproducible across runs.

- drop_nonsig:

  Logical (default `FALSE`). If `TRUE`, non-significant links (p \>
  0.05) are removed from the plots; the returned stats data frame is
  unaffected.

- layout:

  Character (default `"gephi"`). Layout passed to the underlying
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  call (e.g. `"gephi"`, `"square"`, `"WGCNA"`).

- layout_module:

  Character (default `"random"`). Module ordering strategy passed
  through to
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md).
  One of `"random"`, `"adjacent"`, `"order"`.

- orientation:

  Character vector (default
  `c("top_right","bottom_right","top_left","bottom_left")`). Which
  heatmap quadrants to draw, in the same order as `env_select`.

- distance:

  Numeric (default `3`). Offset between the central network's outer
  boundary and the env heatmaps. Positive pushes heatmaps outward; `0`
  places them flush; negative values pull them inward and may overlap
  the network.

- r:

  Numeric (default `6`). Effective radius for scaling the central
  network.

- HeatmapScale:

  Numeric (default `1`). Global scale for the overall heatmap layout.
  `>1` enlarges, `<1` shrinks.

- SigLineAlpha:

  Numeric in `[0, 1]` (default `0.5`). Transparency for module-env link
  segments.

- HeatmapLabelSize:

  Numeric (default `5`). Text size for the heatmap row/column labels.

- HeatmapSigSize:

  Numeric (default `5`). Text size for the significance marks (`*`,
  `**`, `***`) inside heatmap tiles.

- HeatmapColorBar:

  `NULL` or list (default `NULL`). Per-quadrant colour palettes. Three
  accepted forms:

  - `NULL`: built-in defaults.

  - Length-2 named list `list(low = ..., high = ...)`: applied to all
    quadrants.

  - List of length `length(orientation)`: each element is either
    `c(low, high)` or `list(low = ..., high = ...)`. Example:
    `list(c("#2166ac","#b2182b"), c("#1b7837","#762a83"), c("#4393c3","#d6604d"), c("#92c5de","#f4a582"))`.

- HeatmapLabelOrient:

  Numeric (default `0`). Rotation angle (in degrees) for heatmap
  row/column labels. Try 45 or 90 to avoid label overlap.

- SigLineWidth:

  Numeric vector of length 2 (default `c(0.5, 2)`). Min / max line width
  for module-env links; mapped from `-log10(p-value)` so smaller p -\>
  thicker line.

- SigLineColor:

  Character vector of length 2 (default `c("#fdbb84", "#d7301f")`).
  Colour gradient for module-env links, mapped from low / high
  correlation (or Mantel r).

- HeatmapPointSize:

  Numeric (default `5`). Point size for the central module anchor where
  the heatmap link lands.

- HeatmapPointFill:

  Character (default `"#de77ae"`). Fill colour for the central module
  anchor point.

- HeatmapTileColor:

  Character or `NA` (default `NA`). Border colour for heatmap tiles.

- HeatmapTileSize:

  Numeric (default `0`). Border line width for heatmap tiles.

- ...:

  Additional arguments forwarded to the underlying
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  network call. Commonly used: `shrink`, `inner_shrink` (intra-module
  compactness, only for `layout = "WGCNA"`), `node_jitter`,
  `module_outline`, `network_outline`, `module_label` (logical or
  character – module labels in ggNetView style), `module_label_size`,
  `module_label_segment_width`, `module_label_segment_alpha`,
  `node_fill_values`, `node_size_range`. Deprecated pre-0.2.0 names
  (`add_outer`, `label`, `fill`, ...) are still accepted with a
  lifecycle warning.

- layout.module:

  **\[deprecated\]** Use `layout_module`.

## Value

A list of length 3:

- `[[1]]`:

  ggplot object with straight link segments (`geom_segment`).

- `[[2]]`:

  ggplot object with curved link segments (`geom_curve`).

- `[[3]]`:

  Data frame of module-env stats (unfiltered, not affected by
  `drop_nonsig`). Columns: `ID` (module name), `Type` (env column name),
  `Correlation`, `Pvalue`, `p_signif`, `spec_block`, `env_block`,
  `method` (`"correlation"` or `"mantel"`). Schema is identical across
  all `relation_method` / `mantel_kind` combinations.

## Details

**Pipeline.**

1.  Read module membership from `graph_obj` (the node attribute selected
    by the package's module column, e.g. `"Modularity"`).

2.  Build a per-sample module summary matrix `spec_df` of shape (samples
    x modules) using either `module_eigengene` (PC1 of the module's OTU
    sub-matrix; recommended) or `module_abundance` (sum / mean of OTU
    abundances inside the module).

3.  For each `(env_block, modules)` pair, compute either a correlation
    or a Mantel test (see **Mantel API** below).

4.  Render the network via
    [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
    (parameters forwarded through `...`), draw four heatmap quadrants
    for env-env correlations, and overlay link segments from module
    anchors to env diagonals.

**Mantel API (shared with
[`gglink_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md)).**
Two algorithms exposed via `mantel_kind`; both go through the helpers in
[`mantel_block_vs_col`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md)
/
[`mantel_pairwise`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md)
so the two top-level functions stay numerically identical:

- `"block_vs_col"` (default, ecological standard): for each module, the
  OTUs that belong to it are pulled out of `otu_mat` (transposed to
  samples x OTUs) and turned into ONE community distance matrix with
  `spec_dist_method`; each env column is turned into its own distance
  matrix with `env_dist_method`; one Mantel test per (module, env_col).

- `"col_vs_col"` (legacy): the module's representative vector (eigengene
  or abundance) is treated as a single variable; its single-column
  distance matrix is tested against each env column's single-column
  distance matrix. Mathematically close to a rank correlation. Kept for
  backwards compatibility / sensitivity comparisons.

Prior versions of this function used the equivalent of `"col_vs_col"`
implicitly. The default has been switched to `"block_vs_col"` (with a
one-time [`message()`](https://rdrr.io/r/base/message.html) on the first
Mantel call) to match the standard ecological interpretation; pass
`mantel_kind = "col_vs_col"` to reproduce the old numbers.

**Output schema is stable across modes.** The returned stats data frame
always has columns
`ID, Type, Correlation, Pvalue, p_signif, spec_block, env_block, method`.
`ID` is the module name in all cases (`"M1"`, `"M2"`, ...), `Type` is
the env column name. This makes `drop_nonsig` and `comparisons_groups`
work identically across `relation_method` / `mantel_kind` combinations.

## Data inputs

The graph carrying module assignments, the env data table that defines
the heatmap quadrants, the OTU abundance matrix used to summarise each
module, and the named list that partitions env into blocks.

## Module representation

How each module's per-sample value is computed before being correlated
with env: either as the eigengene (PC1 of the OTU sub-matrix) or as a
summary (sum / mean) of within-module abundances.

## Statistics – correlation

Parameters that govern the env-env tile correlations and (when
`relation_method = "correlation"`) the module-env link correlations: the
correlation method and missing-value handling.

## Statistics – Mantel

Parameters used only when `relation_method = "mantel"`: the Mantel
variant, the dissimilarity / distance metrics, the Mantel correlation
method, and the permutation count.

## What gets analysed / drawn

Filters and selectors that decide what ends up on the plot: dropping
non-significant links, and which heatmap quadrants are rendered.

## Geometry

Spatial parameters that position the env heatmaps relative to the
central network and select the network's own layout: heatmap offset,
network radius, overall heatmap scale, and the layout/module-ordering
choices forwarded to
[`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md).

## Heatmap aesthetics

Visual styling of the env-env heatmap tiles: per-quadrant colour
palettes, label and significance-mark sizes, label rotation, the central
anchor point on each heatmap, and tile border styling.

## Link line aesthetics

Visual styling of the module-env link segments: line-width range (mapped
from p-value), colour gradient (mapped from correlation / Mantel r), and
overall transparency.

## Forwarded to ggNetView

Extra arguments captured via `...` and forwarded verbatim to the
underlying
[`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
call. Use them to customise the central network's appearance (labels,
fills, jitter, outer rings, point sizes) without leaving this wrapper.

## See also

[`gglink_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md)
for the spec-select counterpart that shares the same Mantel API;
[`mantel_block_vs_col`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md),
[`mantel_pairwise`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md)
for the underlying Mantel implementations;
[`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
for the central network rendering and the parameters forwarded via
`...`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: correlation, eigengene as module summary, four quadrants.
p <- ggnetview_modularity_heatmaps(
  graph_obj  = g,
  env        = env,
  otu_mat    = otu_mat,
  env_select = list(Env01 = 1:5, Env02 = 6:10,
                    Env03 = 11:15, Env04 = 16:20)
)
p[[1]]            # straight links
head(p[[3]])      # stats data frame

# Ecologically standard Mantel: one test per (module, env_col).
p2 <- ggnetview_modularity_heatmaps(
  graph_obj        = g,
  env              = env,
  otu_mat          = otu_mat,
  env_select       = list(Env01 = 1:5, Env02 = 6:10,
                          Env03 = 11:15, Env04 = 16:20),
  relation_method  = "mantel",
  mantel_kind      = "block_vs_col",
  spec_dist_method = "bray",
  env_dist_method  = "euclidean",
  permutations     = 999
)

# Reproduce legacy column-vs-column Mantel results.
p3 <- ggnetview_modularity_heatmaps(
  graph_obj       = g,
  env             = env,
  otu_mat         = otu_mat,
  env_select      = list(Env01 = 1:5, Env02 = 6:10,
                         Env03 = 11:15, Env04 = 16:20),
  relation_method = "mantel",
  mantel_kind     = "col_vs_col"
)
} # }
```
