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

  Data Frame A data frame or matrix containing environmental variables.
  Each column represents an environmental factor.

- spec:

  Data Frame A data frame or matrix containing species abundance or
  trait data. Each column represents a species or taxonomic unit.

- env_select:

  Optional list specifying the column indices (or names) of
  environmental variables to include in each environmental block. Each
  list element corresponds to one quadrant of the heatmap layout.

- spec_select:

  Optional list specifying column indices (or names) of species to
  include. If multiple elements are provided, they define separate
  species clusters in the visualization.

- spec_layout:

  Character or character vector (default = `"circle_outline"`). Spatial
  arrangement of species nodes for each central network. If a single
  value, applied to all spec blocks. If a vector, length must match the
  number of spec blocks; each element specifies the layout for that
  block in order. Valid options: `"circle_outline"`,
  `"diamond_outline"`, `"rectangle_outline"`, `"square_outline"`.

- spec_orientation:

  Character spec_oritation. Options include: "up","down","left","right"

- spec_relation:

  Logical (defalt = TRUE) Whether to compulate the ralationship of spec

- relation_method:

  Character Method for computing relationships between species and
  environmental factors. Options are \`"correlation"\` or \`"mantel"\`.

- cor.method:

  Character Correlation method to use when \`relation_method =
  "correlation"\`. One of \`"pearson"\`, \`"kendall"\`, or
  \`"spearman"\`.

- cor.use:

  Character Method for handling missing values in correlation
  computation. One of \`"everything"\`, \`"all"\`, \`"complete"\`,
  \`"pairwise"\`, or \`"na"\`.

- mantel.method:

  Character Type of Mantel test to use when \`relation_method =
  "mantel"\`. Options include \`"mantel"\`, \`"mantel.partial"\`,
  \`"mantelhaen.test"\`, and \`"mantel.correlog"\`.

- mantel.method2:

  Character Correlation coefficient used in the Mantel test. One of
  \`"pearson"\`, \`"kendall"\`, or \`"spearman"\`.

- mantel.alternative:

  Character Alternative hypothesis for Mantel test. One of
  \`"two.sided"\`, \`"less"\`, or \`"greater"\`.

- drop_nonsig:

  Logical if \`TRUE\`, non-significant correlations are dropped from the
  final visualization.

- comparisons:

  Logical (default = TRUE). Whether to perform species–environment
  correlation or Mantel analysis. If `FALSE`, no spec–env links are
  computed or drawn.

- comparisons_groups:

  List or NULL (default = NULL). When `comparisons = TRUE`, constrains
  which (env_block, spec_block) pairs are analyzed. Each element must be
  a length-2 character vector: `c(env_block_name, spec_block_name)`,
  e.g. `list(c("Env01", "Spec01"), c("Env02", "Spec01"))`. Block names
  must match `names(env_select)` and `names(spec_select)`. If `NULL`,
  all env–spec block pairs are analyzed (default).

- shape:

  Intrger Integer or numeric specifying the shape of species nodes in
  the plot (passed to \`geom_point()\`).

- distance:

  Numeric the offset distance between central nodes and the
  environmental heatmaps.

- HeatmapLabelSize:

  Numeric (default = 5) Text size for heatmap axis labels (ID/Type).

- HeatmapSigSize:

  Numeric (default = 5) Text size for significance symbols (e.g. \`\*\`,
  \`\*\*\`, \`\*\*\*\`) inside heatmap tiles.

- HeatmapColorBar:

  NULL or list Controls the colorbar palettes used by each heatmap
  quadrant.

  \- If \`NULL\`, the built-in default palettes are used (same as
  current behavior). - If a list of length 2 with names \`low\` and
  \`high\`, each should be a character vector of colors (recycled if
  shorter) used across quadrants. - If a list of length equal to the
  number of quadrants, each element should be either \`c(low, high)\` or
  \`list(low=..., high=...)\` for that quadrant (in order). - example
  \`HeatmapColorBar = list( c("#2166ac", "#b2182b"), \# 第1个象限
  low/high c("#1b7837", "#762a83"), \# 第2个 c("#4393c3", "#d6604d"), \#
  第3个 c("#92c5de", "#f4a582") \# 第4个 )\`

- HeatmapLabelOrient:

  Numeric (default = 0) Rotation angle (in degrees) for heatmap axis
  labels (ID/Type). Use this to avoid overlap of the top/bottom labels;
  e.g. 45 or 90.

- SigLineWidth:

  Numeric vector of length 2 (default = c(0.5, 2)) Controls the minimum
  and maximum line width for species–environment links, scaled by
  significance (p-value). Smaller p leads to thicker lines.

- SigLineColor:

  Character vector (length 2, default = c("#fdbb84", "#d7301f")) Colors
  used for the species–environment link color gradient, corresponding to
  low and high correlation values respectively.

- HeatmapPointSize:

  Numeric (default = 5) Point size for the heatmap diagonal nodes.

- CorePointSize:

  Numeric (default = 8.5) Point size for the central species nodes.

- HeatmapPointFill:

  Character (default = "#de77ae") Fill color for heatmap diagonal
  points.

- CorePointFill:

  Character (default = "#41b6c4") Fill color for central species nodes.

- HeatmapTileColor:

  Character or NA (default = NA) Border color for heatmap tiles (passed
  to \`geom_tile(colour=...)\`).

- HeatmapTileSize:

  Numeric (default = 0) Border line width for heatmap tiles (passed to
  \`geom_tile(size=...)\`).

- HeatmapScale:

  Numeric (default = 1) Global scale factor for the overall heatmap
  size. Values \> 1 enlarge the whole heatmap layout, while values \< 1
  shrink it.

- SigLineAlpha:

  Numeric (default = 0.5) Transparency for species-network to heatmap
  link lines. Must be between 0 and 1.

- fontsize:

  Numeric (default = 5) (Deprecated) Use \`HeatmapLabelSize\` instead.

- orientation:

  Character Character vector defining which heatmap quadrants to
  display. Can include any combination of \`"top_right"\`,
  \`"bottom_right"\`, \`"top_left"\`, and \`"bottom_left"\`.

- r:

  Numeric radius of the central species layout (in plot units).

- group_layout:

  Character (default = "circle"). Arrangement of multiple species
  networks when `spec_select` has multiple elements. Options:
  `"circle"`, `"row"`, `"column"`, `"square"`, `"diamond"`,
  `"triangle"`, `"triangle_down"`, `"snake"`.

- anchor_dist:

  Numeric (default = 6). Distance between species networks when multiple
  `spec_select` blocks are used.

- scale_networks:

  Logical (default = TRUE). If `TRUE`, normalize each species network to
  the same scale (radius `r`) before placing on anchors. If `FALSE`, `r`
  is the minimum network radius; larger networks scale proportionally by
  node count.

- nrow:

  Integer (default = NULL). Number of rows for `group_layout = "row"`,
  `"column"`, or `"snake"`.

- ncol:

  Integer (default = NULL). Number of columns for
  `group_layout = "row"`, `"column"`, or `"snake"`.

## Value

A list of length 3: - \[\[1\]\]: ggplot object with straight link
segments. - \[\[2\]\]: ggplot object with curved link segments. -
\[\[3\]\]: data.frame of full species-environment correlation statistics
(unfiltered, not affected by `drop_nonsig`), with columns `ID`, `Type`,
`Correlation`, `Pvalue`, `spec_block`, `env_block`, and `method` (e.g.
`"correlation"` or `"mantel"`).

## Examples

``` r
NULL
#> NULL
```
