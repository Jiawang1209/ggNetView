# Visualize network with module-environment heatmaps

Combines ggNetView network visualization with four environmental
correlation heatmaps. Module-environment relationships are computed
using module eigengenes (recommended) or module relative abundance, then
correlated with environmental factors. Links connect module centroids to
heatmap diagonals.

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
  drop_nonsig = FALSE,
  layout = "gephi",
  layout.module = c("random", "adjacent", "order"),
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
  ...
)
```

## Arguments

- graph_obj:

  A graph object from `build_graph_from_mat` or `build_graph_from_df`.

- env:

  Data frame or matrix of environmental variables (samples as rows).

- otu_mat:

  Numeric matrix. Rows = OTUs/ASVs, columns = samples. Used to compute
  module eigengenes or module abundance. Must align with `graph_obj`
  node names (rownames) and `env` sample IDs.

- env_select:

  Named list. Column indices or names for each env block. Each block
  corresponds to one heatmap quadrant. E.g.
  `list(Env01 = 1:5, Env02 = 6:10, Env03 = 11:15, Env04 = 16:20)`.

- module_index:

  Character. How to represent each module for env correlation:
  `"eigengene"` (default, PC1 of module OTUs) or `"abundance"` (sum/mean
  of module OTU abundances).

- abundance_type:

  Character. When `module_index = "abundance"`, use `"sum"` or `"mean"`.

- relation_method:

  `"correlation"` or `"mantel"`.

- cor.method:

  Correlation method: `"pearson"`, `"kendall"`, `"spearman"`.

- cor.use:

  Handling of missing values in correlation.

- mantel.method2:

  Correlation for Mantel test.

- drop_nonsig:

  Logical. Drop non-significant links from plot.

- layout:

  Character. Layout for ggNetView (e.g. `"gephi"`, `"square"`).

- orientation:

  Character vector. Heatmap quadrants: `"top_right"`, `"bottom_right"`,
  `"top_left"`, `"bottom_left"`.

- distance:

  Numeric. Gap between the scaled network boundary and the environmental
  heatmaps.

- r:

  Numeric. Effective radius for scaling the central network.

- HeatmapScale:

  Numeric (default = 1). Global scale factor for the overall heatmap
  size. Values \> 1 enlarge the whole heatmap layout.

- SigLineAlpha:

  Numeric (default = 0.5). Transparency for module–heatmap link lines.
  Must be between 0 and 1.

- SigLineWidth:

  Numeric vector of length 2 (default = c(0.5, 2)). Min and max line
  width for module–heatmap links. Line width is mapped from
  `-log10(p-value)`: smaller p (more significant) → thicker line. E.g.
  p=0.05→1.3, p=0.01→2, p=0.001→3; values are scaled to this range.

- SigLineColor:

  Character vector of length 2. Colors for link gradient (low and high
  correlation).

- ...:

  Additional arguments passed to layout and network rendering,
  including: `layout`, `layout.module`, `shrink`, `jitter`, `add_outer`,
  `add_group_outer`, `label` (logical or character: module labels in
  ggNetView style), `labelsize`, `labelsegmentsize`,
  `labelsegmentalpha`, `fill`, `color`, `pointsize`.

## Value

A list: `[[1]]` ggplot with straight links, `[[2]]` ggplot with curved
links, `[[3]]` data frame of module-env correlation stats.
