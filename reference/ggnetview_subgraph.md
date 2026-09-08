# Visualize a network with a magnified (zoomed) module subgraph

Reproduces the common "local-magnification" network figure: the full
network is drawn on the left (optionally with the selected module(s)
outlined), an arrow points to the right, and the extracted module
subgraph is redrawn as its own panel, with an optional node / edge /
component summary in its subtitle.

## Usage

``` r
ggnetview_subgraph(
  graph_obj,
  select_module,
  full_layout = "gephi",
  sub_layout = "same",
  full_args = list(),
  sub_args = list(),
  sub_fill = NULL,
  sub_node_size_range = c(4, 10),
  arrow = TRUE,
  show_stats = TRUE,
  full_title = "Full Network",
  sub_title = NULL,
  widths = c(1, 0.12, 0.62),
  seed = 1115,
  sub_pointsize = deprecated()
)
```

## Arguments

- graph_obj:

  A `tbl_graph` object from
  [`build_graph_from_mat`](https://jiawang1209.github.io/ggNetView/reference/build_graph_from_mat.md)
  /
  [`build_graph_from_df`](https://jiawang1209.github.io/ggNetView/reference/build_graph_from_df.md)
  (or any other `build_graph_from_*`). Its node table must contain a
  `Modularity` column, exactly as produced by those builders.

- select_module:

  Character or numeric vector. The module name(s) (from
  `levels(Modularity)`) to extract into the magnified panel. Multiple
  modules are allowed; they are combined into a single induced subgraph.

- full_layout:

  Character (default `"gephi"`). Layout used for the full network panel;
  any layout accepted by
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md).

- sub_layout:

  Character (default `"same"`). Layout of the magnified subgraph panel.
  `"same"` inherits the full network's node coordinates (a true zoom,
  identical relative positions). Any other value is treated as a
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  layout name (e.g. `"circle"`, `"gephi"`, `"fr"`) and re-lays-out the
  subgraph independently. When `select_module` names several modules,
  the subgraph contains exactly those modules, so the multipartite
  layouts become applicable: `"bipartite_gephi_layout"` for 2 modules,
  `"tripartite_gephi_layout"` for 3, `"quadripartite_gephi_layout"` for
  4 and `"pentapartite_gephi_layout"` for 5 (each requires the module
  count to match). Unused module levels are dropped automatically so the
  count is exact.

- full_args, sub_args:

  Named lists of extra arguments passed through to
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  for the full-network and subgraph panels respectively (e.g.
  `full_args = list(module_label = TRUE, node_size_range = c(2, 6))`).
  Deprecated pre-0.2.0 names are still accepted with a lifecycle
  warning. Anything you can pass to
  [`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  can go here, including `module_outline = TRUE` to outline modules
  (with `module_outline_q`, `module_outline_alpha`,
  `module_outline_width`, etc.). On the full network this is drawn by
  `ggNetView` and outlines every module. On the subgraph the outline is
  drawn by this function (for *any* `sub_layout`, including `"same"`)
  using the full-network palette, so the outer colour always matches the
  node / full-network colours; it outlines the selected module(s). For
  `sub_layout = "same"` the subgraph is drawn directly, so only
  `edge_color`, `edge_alpha` and the `module_outline` styling keys in
  `sub_args` take effect there.

- sub_fill:

  Optional single colour used to recolour every node of the magnified
  subgraph (as in the classic teal "extracted module" figure). When
  `NULL` (default) the subgraph keeps its original module colour, so it
  visually matches the module in the full network.

- sub_node_size_range:

  Numeric length-2 vector (default `c(4, 10)`). Point size range for the
  enlarged subgraph.

- arrow:

  Logical (default `TRUE`). Whether to draw the connecting arrow panel
  between the full network and the magnified subgraph.

- show_stats:

  Logical (default `TRUE`). Whether to show a nodes / edges / components
  summary in the subgraph panel's subtitle.

- full_title:

  Character (default `"Full Network"`). Title of the left panel.

- sub_title:

  Character or `NULL`. Title of the right panel. When `NULL` (default) a
  title of the form `"Extracted Subgraph (Module X)"` is generated
  automatically.

- widths:

  Numeric length-3 vector (default `c(1, 0.12, 0.62)`) giving the
  relative widths of the full-network, arrow and subgraph panels; the
  subgraph panel is deliberately smaller than the full network. When
  `arrow = FALSE` the middle entry is ignored.

- seed:

  Integer (default `1115`). Seed forwarded to
  [`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  for both panels to keep the figure reproducible.

- sub_pointsize:

  **\[deprecated\]** Use `sub_node_size_range`.

## Value

A `patchwork` object (a composed `ggplot`) that can be printed, further
modified with `+`, or saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The function is a thin, deterministic wrapper that reuses the package's
own machinery end to end: the full network is rendered by
[`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
(so network construction, module colouring, layouts and seeds are
identical to the rest of the package), the subgraph is extracted by
[`get_subgraph`](https://jiawang1209.github.io/ggNetView/reference/get_subgraph.md),
and the panels are composed with `patchwork`. By default
(`sub_layout = "same"`) the magnified panel inherits the exact node
coordinates of the full network, so it is a true zoom of that module
rather than a re-layout. No non-CRAN dependency (e.g. `ggmagnify`) is
required.

## See also

[`get_subgraph`](https://jiawang1209.github.io/ggNetView/reference/get_subgraph.md),
[`ggNetView`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build a microbial co-occurrence network exactly as in a normal
# ggNetView workflow.
data(otu_rare_relative)
data(tax_tab)
obj <- build_graph_from_mat(
  mat              = otu_rare_relative,  # variables x samples abundance matrix
  transfrom.method = "none",             # pre-correlation transform
  method           = "WGCNA",            # WGCNA::corAndPvalue backend
  cor.method       = "pearson",          # Pearson correlation
  proc             = "BH",               # Benjamini-Hochberg correction
  r.threshold      = 0.7,                # |r| edge cutoff
  p.threshold      = 0.05,               # adjusted p-value cutoff
  node_annotation  = tax_tab             # taxonomy joined onto nodes
)

# The module is magnified keeping the SAME layout as the full network (a
# true zoom). Full-panel styling is passed through via `full_args`.
full_style <- list(
  layout_module   = "adjacent",  # neighbouring modules close together
  node_size_range = c(1, 5),     # node size range
  center        = FALSE,       # do not pull nodes to module centre
  shrink        = 0.9,         # compact layout
  edge_alpha    = 0.2,         # edge transparency
  edge_color    = "#d9d9d9"    # edge colour
)
ggnetview_subgraph(obj, select_module = "1", full_args = full_style)

# Outline modules via ggNetView's own module_outline, on BOTH panels.
ggnetview_subgraph(
  obj, select_module = "1",
  full_args = c(full_style, list(module_outline = TRUE)),
  sub_args  = list(module_outline = TRUE)
)

# Re-lay-out the subgraph as a clean circle, recoloured teal.
ggnetview_subgraph(
  obj,
  select_module = "1",
  sub_layout    = "circle",
  sub_fill      = "#2C7C82",
  full_args     = full_style
)

# Select THREE modules and lay the subgraph out as a tripartite network.
ggnetview_subgraph(
  obj,
  select_module = c("1", "2", "3"),
  sub_layout    = "tripartite_gephi_layout",
  full_args     = full_style
)
} # }
```
