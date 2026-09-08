# Visualize multi-orientation environmental-species correlation heatmaps (adaptive sizing)

An improved version of
[`gglink_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md)
that sizes each heatmap quadrant independently according to its own
number of variables. All tiles share the same size; larger env blocks
simply extend further. The central species network is always kept at the
centre of the canvas.

## Usage

``` r
gglink_heatmaps_2(
  env,
  spec,
  env_select = NULL,
  spec_select = NULL,
  spec_layout = "circle_outline",
  spec_orientation = c("up", "down", "left", "right"),
  spec_relation = TRUE,
  relation_method = c("correlation", "mantel"),
  cor.method = c("pearson", "kendall", "spearman"),
  cor.use = c("everything", "all", "complete", "pairwise", "na"),
  mantel.method = c("mantel", "mantel.partial", "mantelhaen.test", "mantel.correlog"),
  mantel.method2 = c("pearson", "kendall", "spearman"),
  mantel.alternative = c("two.sided", "less", "greater"),
  mantel.seed = 1115,
  drop_nonsig = FALSE,
  comparisons = TRUE,
  comparisons_groups = NULL,
  shape = 22,
  distance = 3,
  HeatmapLabelSize = 5,
  HeatmapSigSize = 5,
  HeatmapColorBar = NULL,
  HeatmapLabelOrient = 0,
  SigLineWidth = c(0.5, 2),
  SigLineColor = c("#fdbb84", "#d7301f"),
  HeatmapPointSize = 5,
  CorePointSize = 8.5,
  HeatmapPointFill = "#de77ae",
  CorePointFill = "#41b6c4",
  HeatmapTileColor = NA,
  HeatmapTileSize = 0,
  HeatmapScale = 1,
  SigLineAlpha = 0.5,
  fontsize = 5,
  orientation = c("top_right", "bottom_right", "top_left", "bottom_left"),
  r = 6,
  group_layout = c("circle", "row", "column", "square", "diamond", "triangle",
    "triangle_down", "snake"),
  anchor_dist = 6,
  scale_networks = TRUE,
  nrow = NULL,
  ncol = NULL
)
```

## Arguments

- env:

  Data frame or matrix. Environmental variables, one column per factor,
  one row per sample. Row order must match `spec`.

- spec:

  Data frame or matrix. Species abundance / trait data, one column per
  species (or taxonomic unit), one row per sample. Row order must match
  `env`.

- env_select:

  Named list (required). Each element gives the column indices or names
  of `env` that form one environmental block; each block becomes one
  heatmap quadrant. The list length must equal `length(orientation)`.
  Block names (used for `comparisons_groups`) come from
  `names(env_select)`, e.g.
  `list(Env01 = 1:14, Env02 = 15:28, Env03 = 29:42, Env04 = 43:56)`.

- spec_select:

  Named list (required). Each element gives the column indices or names
  of `spec` that form one species block. Each block is rendered as one
  central network (or one collapsed point if `spec_collapse = TRUE`).
  Block names come from `names(spec_select)`, e.g.
  `list(Spec01 = 1:15, Spec02 = 16:30)`.

- spec_layout:

  Character or character vector (default `"circle_outline"`). Shape of
  the per-block node layout. Length 1 applies to all blocks; a vector
  must have length equal to `length(spec_select)` and is matched
  element-wise. Valid values: `"circle_outline"`, `"diamond_outline"`,
  `"rectangle_outline"`, `"square_outline"`. Ignored when
  `spec_collapse = TRUE`.

- spec_orientation:

  Character (default `"up"`). Base orientation passed to the per-block
  layout function. One of `"up"`, `"down"`, `"left"`, `"right"`.

- spec_relation:

  Logical (default `TRUE`). Whether to compute within-block
  species-species correlations to drive the per-block layout (e.g.
  modularity of `"circle_outline"`). Set `FALSE` for a geometry-only
  layout. Ignored when `spec_collapse = TRUE`.

- relation_method:

  Character (default `"correlation"`). One of `"correlation"` or
  `"mantel"`.

- cor.method:

  Character (default `"pearson"`). Correlation method used by
  [`psych::corr.test`](https://rdrr.io/pkg/psych/man/corr.test.html) for
  env-env (the heatmap tiles) and, when
  `relation_method = "correlation"`, for spec-env links. One of
  `"pearson"`, `"kendall"`, `"spearman"`.

- cor.use:

  Character (default `"everything"`). Missing-value handling for
  [`psych::corr.test`](https://rdrr.io/pkg/psych/man/corr.test.html).
  One of `"everything"`, `"all"`, `"complete"`, `"pairwise"`, `"na"`.

- mantel.method:

  Character. Reserved for future use – currently accepted for backwards
  compatibility but **not** consumed by the active code path (only
  [`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
  via `mantel.method2` is used).

- mantel.method2:

  Character (default `"pearson"`). Correlation coefficient passed to
  [`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
  as its `method` argument. One of `"pearson"`, `"kendall"`,
  `"spearman"`.

- mantel.alternative:

  Character. Same status as `mantel.method` – accepted but not consumed.

- mantel.seed:

  Integer (default `1115`). Seed forwarded to the Mantel helpers so the
  permutation p-values (and hence the significance stars /
  solid-vs-dashed links) are reproducible across runs.

- drop_nonsig:

  Logical (default `FALSE`). If `TRUE`, non-significant links
  (`Pvalue > sig_threshold`) are removed from the plot but kept in the
  returned stats data frame.

- comparisons:

  Logical (default `TRUE`). Master switch for spec-env analysis. `FALSE`
  skips all spec-env stats and links (only the env-env heatmaps remain).

- comparisons_groups:

  List or NULL (default `NULL`). When `comparisons = TRUE`, restricts
  which (env_block, spec_block) pairs are computed and drawn. Each
  element is a length-2 character vector
  `c(env_block_name, spec_block_name)`; names must match
  `names(env_select)` / `names(spec_select)`. `NULL` means "all pairs".
  Example: `list(c("Env01","Spec01"), c("Env02","Spec02"))`.

- shape:

  Integer (default `22`). Reserved; the rendered point shapes are
  currently hard-coded to `21` internally.

- distance:

  Numeric (default `3`). Offset added between the central node group and
  the env heatmaps. Positive pushes heatmaps outward; `0` places them
  flush against the central group; **negative values are allowed** and
  will pull heatmaps inward (may visually overlap the central points). A
  [`message()`](https://rdrr.io/r/base/message.html) is emitted on
  negative values, and a
  [`warning()`](https://rdrr.io/r/base/warning.html) is emitted if the
  value is so negative that an anchor coordinate becomes \\\le 0\\ (the
  heatmap will then flip to the opposite side).

- HeatmapLabelSize:

  Numeric (default `5`). Text size for heatmap row/column labels
  (ID/Type).

- HeatmapSigSize:

  Numeric (default `5`). Text size for the significance marks (`*`,
  `**`, `***`) inside heatmap tiles.

- HeatmapColorBar:

  `NULL` or list (default `NULL`). Per-quadrant colour palettes. Three
  accepted forms:

  - `NULL`: use built-in defaults.

  - Length-2 named list `list(low = ..., high = ...)`: applied to all
    quadrants (each value can be a vector that is recycled).

  - List of length `length(orientation)`: each element is either
    `c(low, high)` or `list(low = ..., high = ...)` for that quadrant in
    order. Example:
    `list(c("#2166ac","#b2182b"), c("#1b7837","#762a83"), c("#4393c3","#d6604d"), c("#92c5de","#f4a582"))`.

- HeatmapLabelOrient:

  Numeric (default `0`). Rotation angle (in degrees) for heatmap
  row/column labels. Use 45 or 90 to avoid label overlap.

- SigLineWidth:

  Numeric vector of length 2 (default `c(0.5, 2)`). Min and max line
  width for the significant spec-env link segments (mapped through
  `link_width_by`). The minimum is also used as the fixed width for
  non-significant background lines.

- SigLineColor:

  Character vector of length 2 (default `c("#fdbb84", "#d7301f")`).
  Colour gradient endpoints (low, high) for the significant spec-env
  link segments (mapped through `link_color_by`).

- HeatmapPointSize:

  Numeric (default `5`). Point size for the diagonal anchor points on
  each heatmap (where the link lines land).

- CorePointSize:

  Numeric (default `8.5`). Point size for the central species nodes (or
  collapsed block points).

- HeatmapPointFill:

  Character vector (default `"#de77ae"`). Fill colour(s) for the heatmap
  diagonal points.

  - Length 1: same colour for all quadrants.

  - Length `length(orientation)`: one colour per quadrant, in the order
    given by `orientation`.

  - Other lengths: recycled (modulo) over quadrants.

- CorePointFill:

  Character vector (default `"#41b6c4"`). Fill colour(s) for the central
  species nodes.

  - Length 1: same colour for everyone.

  - Length `length(spec_select)`: one colour per spec block, in the
    order of `names(spec_select)`.

  - Other lengths: recycled (modulo) over spec blocks.

- HeatmapTileColor:

  Character or `NA` (default `NA`). Border colour for heatmap tiles
  (passed to `geom_tile(colour = ...)`).

- HeatmapTileSize:

  Numeric (default `0`). Border line width for heatmap tiles.

- HeatmapScale:

  Numeric (default `1`). Global scale for the overall heatmap layout
  (tile spacing). `>1` enlarges, `<1` shrinks.

- SigLineAlpha:

  Numeric in `[0, 1]` (default `0.5`). Transparency of spec-env link
  segments (applied to both the significant and non-significant layers).

- fontsize:

  Numeric (default `5`). Deprecated. Use `HeatmapLabelSize` (which now
  also drives the central species node label size).

- orientation:

  Character vector (default
  `c("top_right","bottom_right","top_left","bottom_left")`). Which
  heatmap quadrants to draw, in the same order as `env_select`'s
  elements (i.e. `env_select[[1]]` -\> `orientation[1]`).

- r:

  Numeric (default 6). Effective radius of a single central network (in
  plot units). When `spec_collapse = TRUE` this only affects how compact
  a single block looks before being collapsed and is essentially
  cosmetic.

- group_layout:

  Character (default `"circle"`). Arrangement of the per-block anchors
  when `spec_select` has multiple elements. One of `"circle"`, `"row"`,
  `"column"`, `"square"`, `"diamond"`, `"triangle"`, `"triangle_down"`,
  `"snake"`, `"arc"`. `"arc"` places anchors on a circular arc whose
  chord has the same row-like footprint; curvature is controlled by
  `group_arc_angle`.

- anchor_dist:

  Numeric (default `6`). Spacing of the anchor layout. For
  `row / column / snake` this is the centre-to-centre distance between
  adjacent anchors; for `circle / square / diamond / triangle` it is the
  radius from origin to each anchor; for `arc` it is the chord-projected
  spacing (chord length = `(n_blocks - 1) * anchor_dist`).

- scale_networks:

  Logical (default `TRUE`). If `TRUE`, normalise each per-block network
  to the same visual radius (`r`); if `FALSE`, `r` is the minimum
  network radius and larger networks scale proportionally to node count.
  Ignored when `spec_collapse = TRUE`.

- nrow, ncol:

  Integer or NULL (default `NULL`). Grid dimensions for
  `group_layout = "row" | "column" | "snake"`. If both are NULL,
  defaults are inferred from the layout choice.

## Value

A list of length 3: - \[\[1\]\]: ggplot object with straight link
segments. - \[\[2\]\]: ggplot object with curved link segments. -
\[\[3\]\]: data.frame of full species-environment correlation statistics
(unfiltered, not affected by `drop_nonsig`), with columns `ID`, `Type`,
`Correlation`, `Pvalue`, `spec_block`, `env_block`, and `method` (e.g.
`"correlation"` or `"mantel"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Adaptive-sized variant of `gglink_heatmaps()`.
p <- gglink_heatmaps_2(
  env  = env,
  spec = spec
)
} # }
```
