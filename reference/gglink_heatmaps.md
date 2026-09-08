# Multi-orientation species-environment correlation / Mantel heatmap

Render a "central species network(s)" surrounded by up to four
environmental-correlation heatmap quadrants, with curved or straight
link segments connecting each spec node (or each spec block) to each env
variable. Supports both Pearson/Spearman/Kendall correlation and
[`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html)
tests, and exposes a unified Mantel API
(`mantel_kind = "block_vs_col" | "col_vs_col"`) shared with
[`ggnetview_modularity_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/ggnetview_modularity_heatmaps.md).
Link line colour and width are user-configurable through expression
strings (`link_color_by`, `link_width_by`); non-significant links are
rendered as a separate flat-grey background layer for visual context.

## Usage

``` r
gglink_heatmaps(
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
  spec_dist_method = "bray",
  env_dist_method = "euclidean",
  mantel_kind = c("block_vs_col", "col_vs_col"),
  permutations = 999L,
  mantel.seed = 1115,
  spec_collapse = FALSE,
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
  SigLineMid = NULL,
  link_color_by = "Correlation",
  link_width_by = "-log10(Pvalue)",
  NonsigLineColor = "grey80",
  NonsigLineType = "dashed",
  sig_threshold = 0.05,
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
    "triangle_down", "snake", "arc"),
  group_angle = 0,
  group_arc_angle = pi/2,
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

- spec_dist_method:

  Character (default `"bray"`). Dissimilarity method (passed to
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html))
  used to convert the spec block into ONE community distance matrix when
  `mantel_kind = "block_vs_col"`.

- env_dist_method:

  Character (default `"euclidean"`). Distance method used to convert
  each env column into its own distance matrix under any `mantel_kind`.

- mantel_kind:

  Character (default `"block_vs_col"`). Which Mantel variant to use; see
  **Details**. Same parameter is exposed in
  [`ggnetview_modularity_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/ggnetview_modularity_heatmaps.md).

- permutations:

  Integer (default `999L`). Number of permutations passed to
  [`vegan::mantel`](https://vegandevs.github.io/vegan/reference/mantel.html).

- mantel.seed:

  Integer (default `1115`). Seed forwarded to the Mantel helpers so the
  permutation p-values (and hence the significance stars /
  solid-vs-dashed links) are reproducible across runs.

- spec_collapse:

  Logical (default `FALSE`). If `TRUE`, each block is rendered as ONE
  labelled point at its anchor position; the point's label is the block
  name from `names(spec_select)`. In this mode, `spec_layout`,
  `spec_relation`, `scale_networks` are all ignored, and link sources
  are always the collapsed point regardless of `relation_method`. Pairs
  naturally with `relation_method = "mantel"` +
  `mantel_kind = "block_vs_col"`.

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

- SigLineMid:

  Character or `NULL` (default `NULL`). If `NULL`, link colour uses
  `ggplot2::scale_colour_gradient(low, high)` – no centred midpoint,
  suitable when the mapped variable is one-sided (e.g. `"Pvalue"`,
  `"-log10(Pvalue)"`, `"abs(Correlation)"`). If a single colour string
  (e.g. `"white"`), link colour switches to
  `ggplot2::scale_colour_gradient2(low, mid, high, midpoint = 0)` –
  recommended when `link_color_by = "Correlation"` and the values span
  both signs.

- link_color_by:

  Character string (default `"Correlation"`). An R expression, written
  as a string, that is parsed via
  [`rlang::parse_expr`](https://rlang.r-lib.org/reference/parse_expr.html)
  and evaluated against the link data frame to produce the values mapped
  to link line colour. Accepts a bare column name (`"Correlation"`,
  `"Pvalue"`) or any numeric expression of the available columns (see
  **Details** for the full list). Common choices:

  - `"Correlation"` – signed effect size; pair with
    `SigLineMid = "white"` for a diverging palette around 0.

  - `"abs(Correlation)"` – unsigned effect size.

  - `"-log10(Pvalue)"` – significance strength (large = more
    significant).

  The expression must yield a numeric vector with one entry per link
  row; otherwise the function errors out and lists the available numeric
  columns. Only significant links (`Pvalue <= sig_threshold`) flow
  through this colour scale; non-significant links are drawn flat in
  `NonsigLineColor`.

- link_width_by:

  Character string (default `"-log10(Pvalue)"`). An R expression
  (string) evaluated against the link data frame to produce the values
  mapped to link line width. Same syntax and evaluation rules as
  `link_color_by`; same set of available columns. Common choices:
  `"-log10(Pvalue)"` (default; significance strength),
  `"abs(Correlation)"`, `"Correlation^2"`. Only significant links
  participate in this width scale; non-significant links are drawn at
  `min(SigLineWidth)` as a flat background.

- NonsigLineColor:

  Character (default `"grey80"`). Fixed colour for non-significant link
  segments (`Pvalue > sig_threshold`). Only used when
  `drop_nonsig = FALSE`.

- NonsigLineType:

  Character (default `"dashed"`). Line type for non-significant link
  segments. Any value accepted by ggplot2's `linetype` aesthetic (e.g.
  `"solid"`, `"dashed"`, `"dotted"`).

- sig_threshold:

  Numeric in `(0, 1)` (default `0.05`). P-value threshold separating
  significant from non-significant links. Used to (a) decide which links
  go through the colour/width scales, (b) drive `drop_nonsig` filtering.

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

- group_angle:

  Numeric (default `0`). Extra rotation applied to the entire anchor
  set, on top of `group_layout`'s default orientation. Accepts radians
  (`|x| <= 2*pi`) or degrees (`|x| > 2*pi`); use
  [`deg`](https://jiawang1209.github.io/ggNetView/reference/deg.md)`(x)`
  to force a small value to mean degrees. Examples: `group_angle = 45`
  tilts a row by 45 degrees; `group_angle = pi/2` turns a row into a
  column.

- group_arc_angle:

  Numeric (default `pi/2`). Only used when `group_layout = "arc"`.
  Central angle subtended by the arc. `0` degenerates to a flat row;
  `pi/2` (= 90 degrees) is a quarter-circle (default); `pi` is a
  half-circle. Negative values flip the arc to the opposite side. Same
  unit auto-detection as `group_angle`.

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

A list of length 3:

- `[[1]]`:

  A `ggplot` object with straight link segments (rendered with
  [`ggplot2::geom_segment`](https://ggplot2.tidyverse.org/reference/geom_segment.html)).

- `[[2]]`:

  A `ggplot` object with curved link segments (rendered with
  [`ggplot2::geom_curve`](https://ggplot2.tidyverse.org/reference/geom_segment.html),
  curvature `0.25`).

- `[[3]]`:

  A data frame of the full spec-env statistics, **unfiltered** (not
  affected by `drop_nonsig` or `sig_threshold`). Columns: `ID`, `Type`,
  `Correlation`, `Pvalue`, `p_signif` (one of `""`, `"*"`, `"**"`,
  `"***"` at the fixed 0.05 / 0.01 / 0.001 cutoffs), `spec_block`,
  `env_block`, `method` (`"correlation"` or `"mantel"`). For
  `mantel_kind = "block_vs_col"`, `ID` is the spec_block name; otherwise
  `ID` is the spec column name.

Both ggplot objects can be post-processed with the usual ggplot2 idioms
(`+ ggplot2::labs(...)`, `+ ggplot2::theme(...)`, etc.).

## Details

**Layered geometry.** The plot is built in three nested layers:

1.  **Single network geometry** (per spec block) controlled by
    `spec_layout` / `spec_orientation` – how nodes inside one block are
    arranged.

2.  **Group geometry** (across spec blocks) controlled by `group_layout`
    / `group_angle` / `group_arc_angle` / `anchor_dist` – how the
    per-block anchors are scattered on the canvas.

3.  **Heatmap geometry** (around the centre) controlled by `orientation`
    / `distance` / `HeatmapScale` – where the four env heatmap quadrants
    sit relative to the central group.

**Two collapse modes for `spec_select`.**

- `spec_collapse = FALSE` (default): each block is drawn as a small
  network and each spec column has its own node. Link sources are
  per-column nodes (correlation / col_vs_col mantel) or block centroids
  (block_vs_col mantel).

- `spec_collapse = TRUE`: each block becomes ONE labelled point at its
  anchor position. Pairs naturally with `relation_method = "mantel"` +
  `mantel_kind = "block_vs_col"`, where one block = one statistical
  unit.

**Mantel API.** Two algorithms exposed via `mantel_kind` and shared with
[`ggnetview_modularity_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/ggnetview_modularity_heatmaps.md):

- `"block_vs_col"` (default, ecological standard, linkET / ggcor style):
  the whole spec block becomes ONE community-distance matrix
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  with `spec_dist_method`); each env column becomes its own distance
  matrix (`env_dist_method`); one Mantel test per (spec_block, env_col).
  See
  [`mantel_block_vs_col`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md).

- `"col_vs_col"` (legacy, kept for sensitivity / backwards
  compatibility): each spec column and each env column is reduced to a
  single-column distance matrix; the resulting "Mantel" is
  mathematically close to a rank correlation. See
  [`mantel_pairwise`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md).

**Angle inputs.** `group_angle` and `group_arc_angle` auto-detect their
unit from the magnitude: `|x| <= 2*pi` is taken as radians, `|x| > 2*pi`
as degrees. Use
[`deg`](https://jiawang1209.github.io/ggNetView/reference/deg.md) to
force a small value to mean degrees (e.g. `deg(5)` = 5 degrees).

**Link line styling.** Spec-env links are rendered in two stacked layers
controlled by `sig_threshold`:

- **Significant layer** (`Pvalue <= sig_threshold`): colour and width
  are mapped from user-supplied expressions `link_color_by` and
  `link_width_by` (each a string like `"Correlation"`,
  `"-log10(Pvalue)"`, or `"abs(Correlation)"`). The colour gradient
  endpoints come from `SigLineColor`; pass `SigLineMid` (e.g. `"white"`)
  to switch to a diverging `scale_colour_gradient2(midpoint = 0)`,
  recommended whenever the mapped variable can be negative (signed
  correlations). The width scale range comes from `SigLineWidth`.

- **Non-significant layer** (`Pvalue > sig_threshold`): drawn flat in
  `NonsigLineColor` (default light grey) with line type `NonsigLineType`
  (default dashed) at the minimum of `SigLineWidth`, purely as visual
  context. This layer is omitted when `drop_nonsig = TRUE`.

Both layers share `SigLineAlpha`. Because non-significant links do not
participate in the colour/width scales, the legends only reflect the
significant layer.

**Available columns for link expressions.** `link_color_by` and
`link_width_by` are evaluated against the link data frame, which
contains the spec-env stat columns plus link-source / link-target
coordinates: `ID`, `Type`, `Correlation`, `Pvalue`, `p_signif`,
`spec_block`, `env_block`, `method`, `is_sig`, plus `x`, `y`, `x_to`,
`y_to`. Most users will only build expressions from `Correlation` and
`Pvalue`.

## Data inputs

Matrices and block selectors that drive the whole computation: the env
and spec data tables plus the named lists that partition each into
visual blocks (heatmap quadrants and central networks).

## Per-block (spec) geometry

Controls how nodes are arranged inside one spec block: the layout shape
(e.g. circle vs diamond outline), its base orientation, whether the
block is collapsed to a single labelled point, and how its visual radius
scales.

## Across-block (group) geometry

Controls how the per-block anchors are scattered across the figure: the
macro arrangement (circle / row / arc / ...), its rotation, the
anchor-to-anchor spacing, and any grid dimensions when applicable.

## Heatmap geometry

Where the env-env heatmaps are placed and scaled relative to the central
spec blocks: which quadrants are drawn, how far the heatmaps sit from
the centre, and an overall size multiplier.

## Statistics – correlation

Parameters that govern the env-env and (when
`relation_method = "correlation"`) spec-env correlation pipeline: the
correlation method, missing-value handling, and the choice between
correlation and Mantel as the spec-env relation.

## Statistics – Mantel

Parameters used only when `relation_method = "mantel"`: the Mantel
variant, the dissimilarity / distance metrics applied to the spec block
and env columns, the Mantel correlation method, and the permutation
count.

## What gets analysed

Switches that decide which (env_block, spec_block) pairs are actually
computed and drawn – the master on/off, the optional restriction list,
and the rule for hiding non-significant links.

## Heatmap aesthetics

Visual styling of the env-env heatmap tiles and labels: per-quadrant
colour palettes, label and significance-mark sizes, label rotation, tile
borders, and the diagonal anchor points where link lines land.

## Central node aesthetics

Visual styling of the central spec nodes (or the collapsed block points
when `spec_collapse = TRUE`): point size and per-block fill colour.

## Link line aesthetics

Visual styling of the spec-env link segments: which numeric expression
drives the colour and width scales for significant links, the colour
endpoints and optional centred midpoint, line widths, transparency, and
the separate styling for non-significant links.

## Deprecated / unused parameters

These arguments are kept in the signature for backwards compatibility
with old call sites but are NOT consumed by the active code path. They
may be removed in a future release.

## See also

[`ggnetview_modularity_heatmaps`](https://jiawang1209.github.io/ggNetView/reference/ggnetview_modularity_heatmaps.md)
for the modularity-based counterpart with the same Mantel API;
[`mantel_block_vs_col`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md),
[`mantel_pairwise`](https://jiawang1209.github.io/ggNetView/reference/mantel_utils.md)
for the underlying Mantel implementations;
[`deg`](https://jiawang1209.github.io/ggNetView/reference/deg.md) to
make a numeric explicitly mean degrees in `group_angle` /
`group_arc_angle`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call (defaults: correlation, no collapse).
p <- gglink_heatmaps(
  env  = env,
  spec = spec,
  env_select  = list(Env01 = 1:14, Env02 = 15:28,
                     Env03 = 29:42, Env04 = 43:56),
  spec_select = list(Spec01 = 1:15, Spec02 = 16:30)
)
p[[1]]   # straight links
p[[2]]   # curved links
head(p[[3]])

# Ecologically standard Mantel + collapse each block to one point.
p2 <- gglink_heatmaps(
  env  = env,
  spec = spec,
  env_select  = list(Env01 = 1:14, Env02 = 15:28,
                     Env03 = 29:42, Env04 = 43:56),
  spec_select = list(Spec01 = 1:15, Spec02 = 16:30),
  relation_method  = "mantel",
  mantel_kind      = "block_vs_col",
  spec_dist_method = "bray",
  env_dist_method  = "euclidean",
  spec_collapse    = TRUE,
  group_layout     = "row",
  group_angle      = 45,            # degrees, auto-detected
  anchor_dist      = 4,
  distance         = -1,            # heatmaps pulled slightly inward
  CorePointFill    = c("#2166ac", "#b2182b"),
  CorePointSize    = 10
)
p2[[1]]

# Customise link aesthetics with R expressions passed as strings: colour
# by signed correlation (with a centred diverging palette), width by the
# absolute effect size, drop non-significant links, and tighten the
# threshold. Any string parseable as an expression of the link data
# frame's columns is accepted.
p3 <- gglink_heatmaps(
  env  = env,
  spec = spec,
  env_select  = list(Env01 = 1:14, Env02 = 15:28,
                     Env03 = 29:42, Env04 = 43:56),
  spec_select = list(Spec01 = 1:15, Spec02 = 16:30),
  link_color_by = "Correlation",
  link_width_by = "abs(Correlation)",
  SigLineColor  = c("#2166ac", "#b2182b"),
  SigLineMid    = "white",
  sig_threshold = 0.01,
  drop_nonsig   = TRUE
)
p3[[2]]
} # }
```
