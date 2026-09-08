# ggNetView 0.2.0

## Argument renaming in `ggNetView()` (lifecycle-managed)

The plotting arguments of `ggNetView()` now follow a tidyverse-style
`element_property` naming scheme. **Every old name keeps working** and emits a
`lifecycle` deprecation warning pointing to the new name; old names will be
removed in a future major release.

| Old | New |
|---|---|
| `fill.by` / `color.by` / `shape` | `node_fill` / `node_color` / `node_shape` |
| `pointsize` / `pointalpha` / `pointstroke` | `node_size_range` / `node_alpha` / `node_stroke` |
| `fill` / `color` | `node_fill_values` / `node_color_values` (+ `edge_color_values`) |
| `jitter` / `jitter_sd` | `node_jitter` / `node_jitter_sd` |
| `pointlabel` / `pointlabelsize` | `node_label` / `node_label_size` |
| `plot_line` | `show_edges` |
| `linecolor` / `mapping_line` | `edge_color` (`mapping_line = TRUE` -> `edge_color = "corr_direction"`) |
| `linealpha` / `curve` / `curvature` | `edge_alpha` / `edge_curve` / `edge_curvature` |
| `label` / `labelsize` / `labelsegmentsize` / `labelsegmentalpha` | `module_label` / `module_label_size` / `module_label_segment_width` / `module_label_segment_alpha` |
| `label_layout` / `label_wrap_width` / `label_outer_pad` | `module_label_layout` / `module_label_wrap` / `module_label_pad` |
| `add_outer` / `q_outer` / `expand_outer` / `bandwidth_scale` | `module_outline` / `module_outline_q` / `module_outline_expand` / `module_outline_bandwidth` |
| `outerwidth` / `outerlinetype` / `outeralpha` | `module_outline_width` / `module_outline_linetype` / `module_outline_alpha` |
| `add_group_outer*` | `network_outline*` (`add_group_outer_linewidth` -> `network_outline_width`) |
| `layout.module` / `group.by` | `layout_module` / `group_by` |
| `remove` / `dropOthers` | `hide_others` / `drop_others` |
| `ring_n` / `nodelabsize` | removed (they had no effect) |

`ggNetView_multi_link()` follows the same scheme (old names kept as
deprecated aliases), plus: `inner_curve*` -> `edge_curve*`,
`link_linewidth_*` -> `link_width_*`, `link_linealpha_*` -> `link_alpha_*`,
`label_offset` / `label_size` -> `group_label_offset` / `group_label_size`,
`mapping_line` / `linecolor` -> `edge_color` (+ new `edge_color_values`),
`fill` / `color` -> `node_fill_values` / `node_color_values`.

The same renaming applies to arguments forwarded through `...` / `full_args` /
`sub_args` in `ggNetView_multi()`, `ggnetview_modularity_heatmaps()` and
`ggnetview_subgraph()` (`ggnetview_subgraph(sub_pointsize)` ->
`sub_node_size_range`; `ggnetview_modularity_heatmaps(layout.module)` ->
`layout_module`). `ggNetView_multi()` now forwards plotting arguments via
`...` instead of mirroring the full `ggNetView()` signature.

## New features

* `gglink_heatmap_triple()` gains `cor.method`, `cor.use`, `env_p_adjust`,
  `link_p_adjust` and `sig_breaks`. The correlation step was previously a fixed
  `psych::corr.test()` call with no way to switch to Spearman/Kendall, control
  missing-value handling, apply a multiple-testing correction, or move the
  significance cut points. All defaults reproduce the previous output exactly.
  `cor.method` / `cor.use` use the same vocabulary as `gglink_heatmaps()`, but
  `cor.use` defaults to `"pairwise"` (what `psych::corr.test()` itself uses and
  what this plot has always done) rather than that function's `"everything"`.
  Note that `psych::corr.test()` reports **raw** p-values in `$p` for the
  two-matrix call whatever its `adjust` argument says, so its `"holm"` default
  never reached the plot; `link_p_adjust` reads `$p.adj` instead and therefore
  actually takes effect.

* `node_fill`, `node_color`, `node_shape`, `node_size`, `edge_color`,
  `edge_width` and `edge_linetype` accept **either a column name (mapping) or a
  literal value (constant)**.
* `edge_width` can map a numeric edge attribute (e.g. `"weight"`) to line
  width (`edge_width_range` controls the range); `edge_linetype` can be
  mapped or fixed.
* `node_size` can map any numeric node column (default `"Degree"`) or be a
  fixed size.
* Solid point shapes (0-20) now work with a mapped `node_fill`: the mapping
  is routed to the colour aesthetic and the legend uses the real shape
  (previously the nodes and the legend rendered grey).
* Node and edge colour palettes are separate (`node_color_values` vs
  `edge_color_values`); the old shared `color` is split automatically.
* `ggNetView()` no longer emits "Coordinate system already present" when
  module labels / outlines are drawn.

## Bug fixes

* `ggNetView_multi_link(drop_others = TRUE)` no longer changes which modules are
  linked across networks. Previously the `"Others"` nodes were removed from each
  graph *before* `compare_modules_by_overlap()` ran, so the hypergeometric test
  was re-run on a drastically smaller universe (e.g. 537 -> 81 shared nodes for
  one group pair). Because module membership is computed per group, a node that
  is a core member of a module in group A but falls into `"Others"` in group B
  was dropped from both, shrinking every real module and flipping module pairs in
  and out of significance -- one example dataset went from 16 module links to 8,
  keeping only 5 of the original pairs. The module-overlap comparison now always
  runs on the complete network, making `drop_others` a display-only switch: the
  cross-group module links are identical for `TRUE` and `FALSE`, and only the
  plotted nodes differ.

* `ggNetView_multi_link()` now excludes the `"Others"` bucket from cross-group
  module links unconditionally. The previous filter (`modA != "Others" | modB !=
  "Others"`) only removed the `Others`-to-`Others` pair and would have drawn a
  link into `"Others"` -- which is a display bucket for every module ranked below
  `top_modules`, not a community, and is not given a module outline -- had such a
  pair ever reached significance.

* `get_sample_subgraph_topology_parallel()` no longer registers a global
  `progressr` handler (which errored with "should not be called with handlers
  on the stack" whenever the function was invoked inside `tryCatch()` /
  `withCallingHandlers()`, e.g. from testthat or user pipelines) and no longer
  permanently overwrites the caller's `future::plan()` in the
  `parallel = FALSE` branch. Progress display is unchanged
  (`progressr::with_progress()`).
* `gglink_heatmap_triple()`: hub anchor coordinates became `NA` whenever
  `ncol(Experiment) > ceiling(ncol(Environment)/2)`, silently dropping every
  hub node and hub link segment from the figure. Anchors are now generated
  for any number of Experiment variables (values in the previously working
  range are unchanged). The function also gained input hardening: `edge$weight`
  defaults to 1, `node$annotation` is derived automatically
  (Experiment/Environment), hub nodes default to the Experiment variables with
  automatic ordering, and mismatched hub counts give an informative error.

## Dependencies

* New import: `lifecycle`.

## Other breaking changes (accumulated since 0.1.0)

* The module outer boundary drawn by `ggNetView(add_outer = TRUE)` is now
  computed via 2D kernel density estimation followed by a Highest-Density-Region
  (HDR) contour, replacing the previous polar-quantile / radial-spline
  algorithm. The user-facing parameters `q_outer` and `expand_outer` are kept,
  but their meaning has been reinterpreted in HDR terms (see `?ggNetView`).
  Visual results from `add_outer = TRUE` therefore differ from prior versions;
  in particular, sparse satellite / outlier nodes will now generally fall
  outside the contour rather than dragging the boundary toward themselves.
* `generateMask_ggnetview()` now returns an additional `polygon_id` column to
  support clusters whose HDR contour has multiple disconnected components.
  Downstream `geom_polygon()` aesthetics in `ggNetView()`,
  `ggnetview_modularity_heatmaps()` and `ggNetView_multi_link()` were updated
  accordingly. Custom callers reusing the mask table should switch
  `group = cluster` to `group = interaction(cluster, polygon_id)`.
* The multipartite layout family (`create_layout_tripartite_*`,
  `create_layout_quadripartite_*`, `create_layout_cross_quadripartite_*`,
  `create_layout_pentapartite_*`) now requires the graph to have exactly
  3 / 4 / 5 modules and raises an error otherwise. The previous behaviour
  silently truncated the layout to the first N modules while leaving the
  remaining nodes in the graph, which produced row-count mismatches in
  downstream `bind_cols()` calls. Filter `graph_obj` to the expected number
  of modules before calling these layouts.
* `build_graph_from_igraph()` now raises an error when an explicit
  `module_attr` is supplied but is not present on the graph, instead of
  silently falling back to community detection.
* `trans_adjacency_matrix_to_df()` now returns a `from`, `to`, `weight`
  data frame (previously the `weight` column was silently dropped even
  though the graph was built with `weighted = TRUE`).

## Bug fixes

* `create_layout_bipartite_layout()` and `create_layout_bipartite_gephi_layout()`
  no longer crash on default invocation: the broken `scale = scale` default
  (which resolved to `base::scale`, a function) is now `scale = TRUE`.
* `build_graph_from_double_mat()` and
  `build_graph_from_double_mat_with_module()` now respect the user's
  `directed` argument. A duplicate `igraph::graph_from_data_frame()` call
  had been silently hard-coding `directed = FALSE`.
* In `build_graph_from_module()`, `build_graph_from_double_mat_with_module()`
  and `build_graph_from_adj_mat_module()`, nodes assigned to `"Others"`
  are no longer turned into `NA` after the `factor()` step. The internal
  `factor_levels` vector now includes `"Others"`.
* The significance-bin `case_when()` chains in `gglink_heatmaps.R`,
  `gglink_heatmaps_2.R`, `build_graph_from_double_mat.R`,
  `build_graph_from_double_mat_with_module.R`,
  `build_graph_from_multi_mat.R` and `use_function.R` no longer return
  `NA` for `Pvalue` values that fall exactly on a boundary (`0.05`,
  `0.01`, `0.001`).
* Functions across the package now handle empty / single-element /
  no-bootstrap inputs gracefully where they previously crashed or
  produced silently wrong output. Affected files include
  `get_location.R`, `get_network_topology.R`,
  `get_network_topology_parallel.R`, `get_geo_neighbors.R`,
  `create_layout_multirings.R`, `use_function.R`, the
  `create_layout_petal*`, `create_layout_square*`,
  `create_layout_rectangle`, `create_layout_gephi`,
  `create_layout_circular_modules_grid_layout`, and the
  `build_graph_*` family (`1:max_model` / `1:top_modules` -> `seq_len(...)`).
* `gglink_heatmap_triple()` no longer hard-codes the layout split
  position at 31; `create_layout2()` now attaches the correct
  non-hub / hub split point as an attribute on the returned layout.
* Added numerical guards across `get_network_topology*`,
  `build_graph_from_consensus.R`, `ggNetView_RMT.R`, `sparcc_matrix.R`
  and `compare_modules_info.R` for `mean(numeric(0))`, `min(numeric(0))`,
  `log(<= 0)` and divisions by zero / negative `phyper` parameters.
* `generateMask_ggnetview()` now returns `NULL` (so the convex-hull
  fallback is used) when the KDE-derived HDR threshold is NA.

# ggNetView 0.1.0

## Initial CRAN release

* Initial submission of `ggNetView` to CRAN.
* Provides a unified, reproducible framework for analyzing and
  visualizing complex biological, ecological, and microbial association
  networks.
* Includes tools for building correlation and co-occurrence networks
  via `WGCNA`, `SpiecEasi`, `SparCC`, and standard correlation methods
  (Pearson, Spearman, Kendall).
* Computes node-level and network-level topological metrics, including
  robustness analyses.
* Supports module-level analyses (modularity detection, Zi-Pi
  classification, sample-level subgraph topology).
* Offers a large family of deterministic layout generators built on
  `ggraph` and `ggplot2` (bipartite, tripartite, quadripartite,
  pentapartite, circular-modules, petal, diamond, heart, star, and
  more), all with reproducible seeds.
* Ships 18 example datasets covering OTU tables, taxonomy tables,
  environmental metadata, PPI networks, and modularity examples.
