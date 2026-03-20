# Visualize multiple networks and show the connections between them to highlight key members within the networks

Provide either `(mat, group_info)` or `graph_obj_list` as input. When
`graph_obj_list` is provided, `mat` and `group_info` are ignored.

## Usage

``` r
ggNetView_multi_link(
  mat = NULL,
  group_info = NULL,
  graph_obj_list = NULL,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  r.threshold = 0.7,
  p.threshold = 0.05,
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
  cor.method = c("pearson", "kendall", "spearman"),
  proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  SpiecEasi.method = c("mb", "glasso"),
  sparcc_R = 20,
  node_annotation = NULL,
  top_modules = 15,
  layout = NULL,
  node_add = 7,
  ring_n = NULL,
  r = 1,
  center = TRUE,
  idx = NULL,
  shrink = 1,
  k_nn = 12,
  push_others_delta = 0,
  layout.module = c("random", "adjacent", "order"),
  group.by = "Modularity",
  fill.by = "Modularity",
  fill = NULL,
  color = NULL,
  pointsize = c(1, 5),
  jitter = FALSE,
  jitter_sd = 0.01,
  mapping_line = FALSE,
  linealpha = 0.25,
  linecolor = "grey70",
  inner_curve = FALSE,
  inner_curvature = 0.12,
  inner_curve_adaptive = TRUE,
  inner_curve_adaptive_range = c(0.7, 1.3),
  inner_curve_adaptive_bins = 7,
  add_outer = "circle",
  q_outer = 0.88,
  expand_outer = 1.02,
  outerwidth = 1.25,
  outerlinetype = 2,
  outeralpha = 0.5,
  link_level = "Module",
  link_curve = FALSE,
  link_curvature = 0.2,
  link_curve_mode = "outward",
  link_curve_adaptive = TRUE,
  link_curve_adaptive_range = c(0.7, 1.3),
  link_curve_adaptive_bins = 7,
  link_color_node = NULL,
  link_color_module = NULL,
  link_linewidth_node = 1,
  link_linewidth_module = 1,
  link_linetype_node = 2,
  link_linetype_module = 1,
  link_linealpha_node = 0.25,
  link_linealpha_module = 0.5,
  dropOthers = FALSE,
  calculate_topology = FALSE,
  comparisons = TRUE,
  comparisons_groups = NULL,
  order = NULL,
  group_layout = "circle",
  scale_groups = TRUE,
  orientation = "up",
  angle = 0,
  anchor_dist = 6,
  layout_anchor_dist = NULL,
  nrow = NULL,
  ncol = NULL,
  sine_period = 4,
  label_offset = 0.2,
  label_size = 4,
  add_group_outer = FALSE,
  add_group_outer_expand = 2,
  add_group_outer_color = "grey50",
  add_group_outer_fill = NULL,
  add_group_outer_fill_alpha = 0.2,
  add_group_outer_linetype = 1,
  add_group_outer_linewidth = 0.5,
  seed = 1115
)
```

## Arguments

- mat:

  Numeric matrix. Required when `graph_obj_list` is `NULL`. A numeric
  matrix with variables (e.g. genes, taxa) in rows and samples in
  columns.

- group_info:

  DataFrame. Required when `graph_obj_list` is `NULL`. The group
  information contains: Sample and Group. If Sample or Group contain
  underscores (`_`), they are automatically replaced with hyphens (`-`)
  to avoid parsing issues; `mat` column names, `order`, and
  `comparisons_groups` are updated accordingly.

- graph_obj_list:

  Named list of `tbl_graph` objects. Alternative to `mat` and
  `group_info`. Each element is a graph object (e.g. from
  `build_graph_from_mat` or `build_graph_from_df`). List names define
  group names (e.g. `list(WT = g1, KO = g2)`). When provided, `mat` and
  `group_info` are ignored. Enables custom pre-built networks.

- transfrom.method:

  Character. Data transformation methods applied before correlation
  analysis. Options include: "none" (raw data), "scale" (z-score
  standardization), "center" (mean centering only), "log2" (log2
  transfrom), "log10" (log10 transfrom), "ln" (natural transfrom ),
  "rrarefy" (random rarefaction using
  [`vegan::rrarefy`](https://vegandevs.github.io/vegan/reference/rarefy.html)),
  "rrarefy_relative" (rarefy then convert to relative abundance).

- r.threshold:

  Numeric. Correlation coefficient threshold; edges are kept only if
  \|r\| \>= r.threshold.

- p.threshold:

  Numeric. \#' Significance threshold for correlations; edges are kept
  only if p \< p.threshold.

- method:

  Character. Relationship analysis methods. Options include: "WGCNA",
  "SpiecEasi", "SPARCC" and "cor".

- cor.method:

  Character. Correlation analysis method. Options include "pearson",
  "kendall", and "spearman".

- proc:

  Character. Correlation p-value adjustment methods. Options include:
  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and
  "none".

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- SpiecEasi.method:

  Character. Method used in `SpiecEasi` network inference; options
  include "mb" and "glasso".

- sparcc_R:

  Integer. Number of bootstrap/permutation replicates for SparCC
  p-values (when `method = "SPARCC"`). Default 20.

- node_annotation:

  Data frame. Optional node annotation table, containing metadata such
  as taxonomy or functional categories.

- top_modules:

  Integer. Number of top-ranked modules to retain for downstream
  visualization or analysis.

- layout:

  Character string. Custom layouts; one of "gephi", "square", "square2",
  "petal", "petal2", "heart_centered","diamond", "star",
  "star_concentric","rectangle, "rightiso_layers" etc.

- node_add:

  Integer (default = 7). Number of nodes to add in each layer of the
  layout.

- ring_n:

  Numeric (default = 7) Numbers of ring in rings layout.

- r:

  Numeric (default = 1). Radius increment for concentric or layered
  layouts.

- center:

  Logical (default = TRUE). Whether to place a node at the center of the
  layout.

- idx:

  Optional. Index of nodes to be emphasized or centered in the layout

- shrink:

  Numeric (default = 1). Shrinkage factor applied to the center points.

- k_nn:

  Numeric (default = 8). Number of nearest neighbors used to build the
  local adjacency graph.

- push_others_delta:

  Numeric (default = 0). Radial offset applied to the "Others" module to
  slightly

- layout.module:

  Character (default = "random") - random : modules are distributed more
  randomly and independently. - adjacent : modules are positioned close
  to each other, minimizing inter-module gaps. - order : modules are
  distributed by order, applicable to \`Bipartite, Tripartite,
  Quadripartite, Multipartite, Pentapartite Layout\`

- group.by:

  Character (default = "Modularity"). Change group for nodes

- fill.by:

  Character (default = "Modularity"). Change fill for nodes

- fill:

  Named vector of colors for node/module fill (e.g.
  `c("M1" = "red", "M2" = "blue")`). If `NULL` (default), uses viridis
  discrete fill scale (`scale_fill_viridis_d`); if provided, uses
  `scale_fill_manual(values = fill)`.

- color:

  Color setting for node/module border. Supports either a single color
  string (fixed border color) or a named vector (module-to-color
  mapping, similar to `fill`). If `NULL`, mapped borders use viridis
  discrete color scale (`scale_color_viridis_d`).

- pointsize:

  Numeric vector of length 2 (default = `c(1, 5)`). The range of point
  size when mapping `Degree` to size. First value is minimum size,
  second is maximum size.

- jitter:

  Logical (default = FALSE). Whether to apply jitter to points.

- jitter_sd:

  Integer (default = 0.1). The standard deviation of the jitter applied
  when \`jitter = TRUE\`.

- mapping_line:

  Logical or Character (default = FALSE). Whether to map line color in
  ggNetView. If a character string is provided, it must be a variable
  name in edge data.

- linealpha:

  Integer (default = 0.25). Change line alpha.

- linecolor:

  Character (default = "grey70"). Change line color.

- inner_curve:

  Logical (default = FALSE). Whether to draw within-group edges as
  curves.

- inner_curvature:

  Numeric (default = 0.12). Curvature for within-group edges when
  `inner_curve = TRUE`.

- inner_curve_adaptive:

  Logical (default = TRUE). Whether to adapt within-group edge curvature
  by edge length when `inner_curve = TRUE`.

- inner_curve_adaptive_range:

  Numeric vector of length 2 (default = c(0.7, 1.3)). Multipliers
  applied to `inner_curvature` for shortest and longest within-group
  edges.

- inner_curve_adaptive_bins:

  Integer (default = 7). Number of bins used to approximate per-edge
  adaptive curvature for within-group edges.

- add_outer:

  Logical or Character (default = "circle"). Add outer boundaries for
  matched modules. Supported values: `"circle"` (use
  [`ggforce::geom_mark_circle`](https://ggforce.data-imaginist.com/reference/geom_mark_circle.html)),
  `"manual"` (use smoothed polygon boundary like `ggNetView`), and
  `"none"` (disable). Logical `TRUE`/`FALSE` are accepted and mapped to
  `"circle"`/`"none"`.

- q_outer:

  Numeric (default = 0.88). Quantile of radial distance used to
  construct the smooth outer boundary when `add_outer = "manual"`.

- expand_outer:

  Numeric (default = 1.02). Global scaling factor applied to the smooth
  outer boundary when `add_outer = "manual"`.

- outerwidth:

  Numeric (default = 1.25). Line width for module outer boundaries.

- outerlinetype:

  Integer or character (default = 2). Linetype for module outer
  boundaries (e.g. 1 = solid, 2 = dashed).

- outeralpha:

  Numeric (default = 0.5). Alpha for module outer boundaries.

- link_level:

  Character (default = "Module"). Cross-group link granularity. One of
  `"None"`, `"Module"`, `"Node"`, `"NodeinModule"`, `"Module&Node"`, or
  `"Module&Node2"`. `"None"` draws no cross-group links. `"Module"`
  links significant module matches (module centroids). `"Node"` links
  shared node names across groups. `"NodeinModule"` links both module
  centroids and nodes within significantly overlapping modules (module
  links show which modules overlap; node links show shared nodes).
  `"Module&Node"` draws both module-level links and node-level links.
  `"Module&Node2"` same as `"Module&Node"` but adds outer boundaries for
  all modules (not only matched ones).

- link_curve:

  Logical (default = FALSE). Whether to draw cross-group links as curves
  (`geom_curve`) instead of straight segments.

- link_curvature:

  Numeric (default = 0.2). Curvature used when `link_curve = TRUE`.

- link_curve_mode:

  Character (default = "outward"). Curve direction strategy used when
  `link_curve = TRUE`. `"outward"` bends links away from the global
  center. `"inward"` bends links toward the global center. `"cross"`
  follows a cross-axis rule: left links bend left, right links bend
  right, upper links bend up, and lower links bend down.

- link_curve_adaptive:

  Logical (default = TRUE). Whether to adapt link curvature by link
  length when `link_curve = TRUE`. Longer links get larger curvature and
  shorter links get smaller curvature.

- link_curve_adaptive_range:

  Numeric vector of length 2 (default = c(0.7, 1.3)). Multipliers
  applied to `link_curvature` for shortest and longest links. The first
  value is the minimum multiplier; the second is the maximum multiplier.

- link_curve_adaptive_bins:

  Integer (default = 7). Number of bins used to approximate per-link
  adaptive curvature.

- link_color_node:

  Character or NULL (default = NULL). Colors for node-to-node
  cross-group links. `NULL` = use default palette. A single value, named
  vector (e.g. `c("WT|KO" = "red")`), or unnamed vector (by pair index).

- link_color_module:

  Character or NULL (default = NULL). Colors for module-to-module
  cross-group links. Same rules as `link_color_node`.

- link_linewidth_node:

  Numeric (default = 1). Line width for node-to-node cross-group links.
  Single value or vector (by pair index/named).

- link_linewidth_module:

  Numeric (default = 1). Line width for module-to-module cross-group
  links. Single value or vector.

- link_linetype_node:

  Integer or character (default = 2). Linetype for node-to-node
  cross-group links (e.g. 2 = dashed, 1 = solid). Single value or
  vector.

- link_linetype_module:

  Integer or character (default = 1). Linetype for module-to-module
  cross-group links. Single value or vector.

- link_linealpha_node:

  Numeric (default = 0.25). Alpha (transparency) for node-to-node
  cross-group links. Single value or vector.

- link_linealpha_module:

  Numeric (default = 0.5). Alpha (transparency) for module-to-module
  cross-group links. Single value or vector.

- dropOthers:

  Logical (default = FALSE). If TRUE, remove nodes in the `"Others"`
  module from each group's `graph_obj` before layout, plotting, and
  module-overlap comparison.

- calculate_topology:

  Logical (default = FALSE). Whether to compute topology for each group
  using
  [`get_network_topology_parallel()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology_parallel.md)
  and
  [`get_sample_subgraph_topology_parallel()`](https://jiawang1209.github.io/ggNetView/reference/get_sample_subgraph_topology_parallel.md).

- comparisons:

  Logical (default = TRUE). Whether to perform cross-group comparisons
  and draw links. If `FALSE`, no module overlap or node comparison is
  done, and no cross-group links are drawn.

- comparisons_groups:

  List or NULL (default = NULL). When `comparisons = TRUE`, constrains
  which group pairs are compared. Each element must be a length-2
  character vector, e.g. `list(c("WT", "OE"), c("WT", "KO"))`. Group
  names must exist in `group_info$Group` or `names(graph_obj_list)`. If
  `NULL`, all pairwise group comparisons are performed.

- order:

  Character vector or NULL (default = NULL). Order of groups for layout
  positions. Groups are placed evenly on a circle: 1st = top, 2nd = next
  clockwise, etc. (e.g. 3 groups = triangle, 4 = square). Must contain
  all unique groups from `group_info$Group` or `names(graph_obj_list)`
  exactly once. If `NULL`, uses `unique(group_info$Group)` or
  `names(graph_obj_list)` order (first occurrence).

- group_layout:

  Character (default = "circle"). Arrangement of groups in the
  multi-group plot. `"circle"`: groups placed evenly on a circle
  (default). `"row"`: groups in a grid, filled row-by-row; uses `nrow`
  and `ncol`. `"column"`: groups in a grid, filled column-by-column;
  uses `nrow` and `ncol`. `"square"`: 4 groups at corners of a square
  (top-left, top-right, bottom-right, bottom-left). `"diamond"`: 4
  groups at top, right, bottom, left (like a rotated square).
  `"triangle"`: 3 groups at vertices of an upright triangle (point up).
  `"triangle_down"`: 3 groups at vertices of an inverted triangle (point
  down). `"snake"`: groups in a grid with snake-like ordering; first row
  left-to-right, second row right-to-left, third row left-to-right, etc.
  Uses `nrow` and `ncol`. `"snake_vertical"`: groups along a smooth sine
  curve; uses `sine_period`. `"snake_vertical_sin"`: same as
  `"snake_vertical"`. `"snake_vertical_cos"`: groups along a smooth
  cosine curve; uses `sine_period`. `"snake_vertical_neg_sin"`: smooth
  `-sin` curve. `"snake_vertical_neg_cos"`: smooth `-cos` curve.
  `"sin"`: vertex-only; group 1 at center; groups 2+ alternate peaks and
  troughs. `"cos"`: vertex-only; group 1 at peak; groups 2+ alternate
  peak and trough. `"-sin"`: vertex-only; group 1 at center; groups 2+
  alternate down, up. `"-cos"`: vertex-only; group 1 at trough; groups
  2+ alternate trough and peak. `"center_pairs"`: group 1 at center
  (y=0); remaining groups form pairs on peak (y=1) and trough (y=-1),
  with one `anchor_dist` gap between center and first pair. Pairs
  alternate right/left of center. If a pair has only one group, it goes
  on top (peak). Empty pair slots keep their positions.

- scale_groups:

  Logical (default = TRUE). Whether to normalize each group to a
  comparable coordinate scale before placing groups on anchors for
  cross-group visual comparison.

- orientation:

  Character string. Custom orientation; one of
  "up","down","left","right".

- angle:

  Integer (default = 0). Change orientation angle.

- anchor_dist:

  Numeric (default = 6). Distance between groups when placing each group
  network on the outer circle in `ggNetView_multi_link`.

- layout_anchor_dist:

  Numeric (default = NULL). Anchor distance passed to the single-group
  layout function (module spacing within each group). If `NULL`, it
  falls back to `anchor_dist` for backward compatibility.

- nrow:

  Integer (default = NULL). Number of rows. Used by: (1)
  `group_layout = "row"`, `"column"`, or `"snake"` for group grid; (2)
  layout functions like `"consensus_module_equal_gephi"` for module
  grid.

- ncol:

  Integer (default = NULL). Number of columns. Used by: (1)
  `group_layout = "row"`, `"column"`, or `"snake"` for group grid; (2)
  layout functions like `"consensus_module_equal_gephi"` for module
  grid.

- sine_period:

  Numeric (default = 4). Groups per wavelength for `snake_vertical*`;
  ignored for `sin`, `cos`, `-sin`, `-cos`.

- label_offset:

  Numeric (default = 0.2). Vertical offset of group labels above each
  group's layout (added to max y).

- label_size:

  Numeric (default = 4). Font size for group labels (Group, Node, Edge,
  etc.).

- add_group_outer:

  Logical (default = FALSE). Whether to add a circle boundary around
  each group (mimics
  [`ggforce::geom_mark_circle`](https://ggforce.data-imaginist.com/reference/geom_mark_circle.html)).

- add_group_outer_expand:

  Numeric (default = 2). Expansion in mm for the group circle to account
  for point size; passed to `geom_mark_circle(expand = ...)`.

- add_group_outer_color:

  Character (default = "grey50"). Color of the group outer circle
  border. A single value applies to all groups; a named vector maps
  group names to colors (e.g. `c("WT" = "blue", "KO" = "red")`); an
  unnamed vector is used by index (recycled if needed).

- add_group_outer_fill:

  Character or NULL (default = NULL). Fill color of the group outer
  circle. `NULL` = no fill (transparent). A single value, named vector,
  or unnamed vector works like `add_group_outer_color`.

- add_group_outer_fill_alpha:

  Numeric (default = 0.2). Alpha (transparency) of the group outer
  circle fill; 0 = fully transparent, 1 = opaque.

- add_group_outer_linetype:

  Integer or character (default = 1). Linetype of the group outer circle
  (e.g. 1 = solid, 2 = dashed).

- add_group_outer_linewidth:

  Numeric (default = 0.5). Line width of the group outer circle.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

- pointstroke:

  Integer (default = 0.3).

## Value

A list containing plot, module-overlap info, link info, group graphs,
and optional topology results.

## Examples

``` r
NULL
#> NULL
```
