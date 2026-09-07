#' Visualize a network with custom layouts
#'
#' @description
#' `ggNetView()` renders a `tbl_graph` produced by one of the
#' `build_graph_from_*()` builders with a deterministic, module-aware layout.
#'
#' Since ggNetView 0.2.0 the plotting arguments follow a tidyverse-style
#' `element_property` naming scheme: `node_*` for nodes, `edge_*` for edges,
#' `module_label_*` for module text labels, `module_outline_*` for the per-module
#' boundary and `network_outline_*` for the whole-network circle. Every old
#' argument name is still accepted with a deprecation warning
#' (see *Deprecated arguments* below).
#'
#' @section Column name or literal value:
#' `node_fill`, `node_color`, `node_shape`, `node_size`, `edge_color`,
#' `edge_width` and `edge_linetype` accept **either** the name of a column in
#' the node / edge table of `graph_obj` (an aesthetic *mapping*) **or** a
#' literal value (a fixed colour, shape code, size, ...). A single string that
#' matches a column name is treated as a mapping; anything else is used as a
#' constant.
#'
#' @section Deprecated arguments:
#' \lifecycle{deprecated} The following arguments were renamed in 0.2.0 and
#' will be removed in a future release. Old names keep working but emit a
#' [lifecycle::deprecate_warn()] message.
#'
#' | Old | New |
#' |---|---|
#' | `fill.by` | `node_fill` |
#' | `color.by` | `node_color` |
#' | `shape` | `node_shape` |
#' | `pointsize` | `node_size_range` |
#' | `pointalpha` | `node_alpha` |
#' | `pointstroke` | `node_stroke` |
#' | `fill` | `node_fill_values` |
#' | `color` | `node_color_values` + `edge_color_values` |
#' | `jitter`, `jitter_sd` | `node_jitter`, `node_jitter_sd` |
#' | `pointlabel`, `pointlabelsize` | `node_label`, `node_label_size` |
#' | `plot_line` | `show_edges` |
#' | `linecolor`, `mapping_line` | `edge_color` (`mapping_line = TRUE` is `edge_color = "corr_direction"`) |
#' | `linealpha` | `edge_alpha` |
#' | `curve`, `curvature` | `edge_curve`, `edge_curvature` |
#' | `label`, `labelsize` | `module_label`, `module_label_size` |
#' | `labelsegmentsize`, `labelsegmentalpha` | `module_label_segment_width`, `module_label_segment_alpha` |
#' | `label_layout`, `label_wrap_width`, `label_outer_pad` | `module_label_layout`, `module_label_wrap`, `module_label_pad` |
#' | `add_outer`, `q_outer`, `expand_outer`, `bandwidth_scale` | `module_outline`, `module_outline_q`, `module_outline_expand`, `module_outline_bandwidth` |
#' | `outerwidth`, `outerlinetype`, `outeralpha` | `module_outline_width`, `module_outline_linetype`, `module_outline_alpha` |
#' | `add_group_outer*` | `network_outline*` (`..._linewidth` is `network_outline_width`) |
#' | `layout.module`, `group.by` | `layout_module`, `group_by` |
#' | `remove`, `dropOthers` | `hide_others`, `drop_others` |
#' | `ring_n`, `nodelabsize` | removed (they had no effect) |
#'
#' @param graph_obj A `tbl_graph` from `build_graph_from_mat()` or any other
#'   `build_graph_from_*()` builder. The network object to be visualized.
#' @param layout Character string naming the layout, dispatched to the
#'   internal `create_layout_<layout>()` function. Popular choices: `"gephi"`,
#'   `"fr"`, `"kk"`, `"stress"`, `"circle"`, `"square"`, `"square2"`,
#'   `"petal"`, `"petal2"`, `"heart_centered"`, `"diamond"`, `"star"`,
#'   `"star_concentric"`, `"rectangle"`, `"rightiso_layers"`, `"WGCNA"`,
#'   `"circlepack"`, `"bipartite_gephi_layout"`, `"tripartite_gephi_layout"`,
#'   `"circular_modules_*"`. An unknown name errors with the full list of
#'   available layouts.
#' @param node_add Integer (default = 7).
#'   Number of nodes to add in each layer of the layout.
#' @param r Numeric (default = 1).
#'   Radius increment for concentric or layered layouts.
#' @param center Logical (default = TRUE).
#'   Whether to place a node at the center of the layout.
#' @param idx Optional.
#'   Index of nodes to be emphasized or centered in the layout.
#' @param shrink Numeric (default = 1).
#'   Shrinkage factor applied to the center points.
#' @param inner_shrink Numeric (default = 1).
#'   Intra-module compactness factor for `layout = "WGCNA"` only.
#'   Controls how tightly nodes fill each module's allocated disc during
#'   the WGCNA bubble-pack layout: at `1` (default) nodes spread to
#'   fill 95 percent of the disc (original behaviour); smaller values
#'   (e.g. `0.65`) contract the FR/uniform fill toward each module
#'   centre, exposing hub/periphery structure and producing visible
#'   inter-module whitespace.  Module disc centres and radii are
#'   invariant under `inner_shrink`; this parameter only affects
#'   the point cloud inside each module.  Ignored by all other layouts.
#' @param k_nn Numeric (default = 12).
#'   Number of nearest neighbors used to build the local adjacency graph
#'   for `layout_module = "adjacent"` / `"order"`.
#' @param push_others_delta Numeric (default = 0).
#'   Radial offset applied to the "Others" module to push it slightly
#'   outward from the rest of the network.
#' @param layout_module Character (default = `"random"`). How modules are
#'   arranged relative to each other:
#'   - `"random"`: modules are distributed more randomly and independently.
#'   - `"adjacent"`: modules are positioned close to each other, minimizing inter-module gaps.
#'   - `"order"`: modules are distributed by order; required by the
#'     `*partite*` and `multirings` layouts.
#' @param group_by Character (default = `"Modularity"`).
#'   Node column used to group nodes into modules (drives module labels and
#'   top-N node labels). `"pie"` switches to the scatterpie rendering path.
#'
#' @param node_fill Character (default = `"Modularity"`).
#'   Node fill: a node column name (mapping) or a single colour (constant).
#'   Categorical columns use `scale_fill_ggnetview()` (or `node_fill_values`);
#'   numeric columns use a blue-red gradient. For shapes without a fill slot
#'   (0-20) the mapping is applied to the point colour instead.
#' @param node_color Character or NULL (default = NULL).
#'   Node border colour: a node column name (mapping) or a single colour
#'   (constant). Numeric columns use `scale_color_gradient()`, otherwise
#'   `scale_color_ggnetview()` / `node_color_values`.
#' @param node_shape Integer or character (default = 21).
#'   ggplot2 point shape (constant) or a node column name (mapping; shapes
#'   21-25 are cycled).
#' @param node_size Character or numeric (default = `"Degree"`).
#'   Node size: a numeric node column name (mapping, scaled to
#'   `node_size_range`) or a single number (constant size).
#' @param node_size_range Numeric length-2 (default = `c(1, 10)`).
#'   Size range used when `node_size` is a mapping.
#' @param node_alpha Numeric (default = 1). Node alpha.
#' @param node_stroke Numeric (default = 0.3). Node border width.
#' @param node_fill_values Named colour vector or NULL (default = NULL).
#'   Manual palette for a categorical `node_fill` mapping, e.g.
#'   `c("1" = "red", "2" = "blue")`. `NULL` uses `scale_fill_ggnetview()`.
#'   When `node_fill` is a module field this palette is also used for module
#'   labels and module outlines.
#' @param node_color_values Named colour vector or NULL (default = NULL).
#'   Manual palette for a categorical `node_color` mapping.
#' @param node_jitter Logical (default = FALSE). Whether to jitter node
#'   positions.
#' @param node_jitter_sd Numeric (default = 0.1). Standard deviation of the
#'   jitter applied when `node_jitter = TRUE`.
#' @param node_label Character or NULL (default = NULL).
#'   Text labels for the top-degree nodes within each module:
#'   `"topN"` (e.g. `"top1"`, `"top7"`) or `"ALL"`.
#' @param node_label_size Numeric (default = 5). Node label size.
#'
#' @param show_edges Logical (default = TRUE). Whether to draw edges.
#' @param edge_color Character (default = `"grey70"`).
#'   Edge colour: an edge column name (mapping) or a single colour
#'   (constant). `"corr_direction"` (added by the correlation builders)
#'   colours positive edges red and negative edges blue and adds the
#'   positive/negative counts to the title. Numeric columns use a blue-red
#'   gradient.
#' @param edge_color_values Named colour vector or NULL (default = NULL).
#'   Manual palette for a categorical `edge_color` mapping.
#' @param edge_width Character or numeric (default = 0.5).
#'   Edge line width: an edge column name (mapping, e.g. `"weight"`, scaled
#'   to `edge_width_range`) or a single number (constant).
#' @param edge_width_range Numeric length-2 (default = `c(0.2, 1.5)`).
#'   Line-width range used when `edge_width` is a mapping.
#' @param edge_linetype Character, integer or NULL (default = 1).
#'   Edge linetype: an edge column name (mapping) or a single linetype
#'   (constant).
#' @param edge_alpha Numeric (default = 0.25). Edge alpha.
#' @param edge_curve Logical (default = FALSE). Draw curved edges.
#' @param edge_curvature Numeric (default = 0.25). Curvature of curved edges.
#'
#' @param module_label Logical or character (default = FALSE).
#'   Controls module text labels and the module legend prefix.
#'   `FALSE`: no module text labels, legend prefix `"Modularity"`.
#'   `TRUE`: draw labels, legend prefix `"Modularity"`.
#'   A string: draw labels and use that string as prefix for module text and
#'   legend labels.
#' @param module_label_size Numeric (default = 10). Module label size.
#' @param module_label_segment_width Numeric (default = 1). Width of the
#'   module label leader segment.
#' @param module_label_segment_alpha Numeric (default = 1). Alpha of the
#'   module label leader segment.
#' @param module_label_layout Character (default = `"two_column"`).
#'   Strategy used to place module text labels when `module_label` is not
#'   `FALSE`. One of:
#'   \itemize{
#'     \item `"two_column"` (default): modules whose centroid `x` is left of
#'       the network median go to a fixed left column, the rest to a fixed
#'       right column, with labels evenly distributed along `y`. Best for
#'       very crowded networks or layouts that read left-to-right (bipartite,
#'       grid).
#'     \item `"two_column_follow"`: 360-degree ring layout with two-segment
#'       "L-shape" leaders. Modules are sorted by their actual angle from the
#'       network centroid and assigned angularly equispaced target angles
#'       around the full circle, so labels are evenly distributed around the
#'       network while preserving the angular neighbour order. Each label is
#'       projected onto an outer ellipse whose semi-axes are
#'       `(1 + module_label_pad)` times the network's half-width /
#'       half-height. The first leader leg runs from the module centroid to
#'       an elbow on the network boundary; the second leg is drawn by
#'       `ggrepel` from the label back to that elbow.
#'     \item `"label_circle"`: like `"two_column_follow"` but every label
#'       sits at its module's true angle on the outer ellipse and `ggrepel`
#'       draws an L-shaped leader (`segment.square = TRUE`).
#'   }
#' @param module_label_wrap Integer or NULL (default = NULL).
#'   If a positive integer, module text labels are wrapped to roughly that
#'   many characters per line via `stringr::str_wrap()`.
#' @param module_label_pad Numeric (default = 0.25).
#'   Fractional outward push of the label anchors: a multiple of the
#'   network's `x`-range (`"two_column"`) or the fractional enlargement of
#'   the outer label-anchor ellipse (`"two_column_follow"`,
#'   `"label_circle"`). Try `0.20` for a tight layout, `0.40` for more
#'   breathing room, `0.55+` when modules sit inside a thick
#'   `network_outline` ring or labels are long.
#'
#' @param module_outline Logical (default = FALSE).
#'   Draw a smooth outer boundary around each module. The boundary is
#'   computed by 2D kernel density estimation followed by a
#'   Highest-Density-Region (HDR) contour: it encloses the densest portion of
#'   each module rather than all of its points. A module whose points fall
#'   into two well-separated clouds may produce two disconnected polygons.
#' @param module_outline_q Numeric (default = 0.88).
#'   HDR coverage of the module outline: the contour is drawn at the density
#'   level whose iso-density region contains a fraction `module_outline_q`
#'   of the module's probability mass. Higher values are more inclusive
#'   (closer to the convex hull). Modules with fewer than 10 nodes bypass the
#'   KDE path and use a convex hull.
#' @param module_outline_expand Numeric (default = 1.02).
#'   Multiplicative scaling applied to each outline polygon from its own
#'   centroid (> 1 expands, < 1 shrinks).
#' @param module_outline_bandwidth Numeric (default = 2).
#'   Multiplier on the robust normal-reference 2D KDE bandwidth used to build
#'   the module outline (> 1 smoother/wider, < 1 tighter).
#' @param module_outline_width Numeric (default = 1). Module outline line width.
#' @param module_outline_linetype Integer or character (default = 1).
#'   Module outline linetype.
#' @param module_outline_alpha Numeric (default = 0.5). Module outline alpha.
#'
#' @param network_outline Logical (default = FALSE).
#'   Draw a circle around the entire network (via
#'   `ggforce::geom_mark_circle()`).
#' @param network_outline_expand Numeric (default = 2).
#'   Expansion in mm of the network circle; passed to
#'   `geom_mark_circle(expand = )`.
#' @param network_outline_color Character (default = `"grey50"`).
#'   Border colour of the network circle.
#' @param network_outline_fill Character or NULL (default = NULL).
#'   Fill colour of the network circle; `NULL` = transparent.
#' @param network_outline_fill_alpha Numeric (default = 0.2).
#'   Fill alpha of the network circle.
#' @param network_outline_linetype Integer or character (default = 1).
#'   Linetype of the network circle.
#' @param network_outline_width Numeric (default = 0.5).
#'   Line width of the network circle.
#'
#' @param hide_others Logical (default = FALSE).
#'   Hide the `"Others"` module at the visualization stage only (post-layout),
#'   so the layout of the remaining modules is unchanged.
#' @param drop_others Logical (default = FALSE).
#'   Remove `"Others"` nodes from `graph_obj` *before* layout, then recompute
#'   layout and plot from the reduced graph.
#' @param orientation Character string (default = `"up"`).
#'   Orientation for directional layouts: `"up"`, `"down"`, `"left"`, `"right"`.
#' @param angle Numeric (default = 0). Rotation angle of the layout.
#' @param scale Logical (default = TRUE).
#'   Whether the `*partite*` layouts scale the module radius.
#' @param anchor_dist Numeric (default = 6).
#'   Distance between modules in the `*partite*` layouts.
#' @param nrow,ncol Integer (default = NULL).
#'   Grid dimensions for `layout = "consensus_module_equal_gephi"` /
#'   `"consensus_module_gephi"`.
#' @param seed Integer (default = 1115). Random seed for reproducibility.
#' @param scale_radius Numeric or NULL (default = NULL).
#'   When non-NULL, scale the layout so the network fits within this radius.
#'   Used by `ggnetview_modularity_heatmaps()` for coordinate alignment.
#' @param return_layout Logical (default = FALSE).
#'   When TRUE, return a list with `$plot` (ggplot) and `$layout_data`
#'   (`graph_ly_final`, `graph_obj`, `ggplot_data`, `module_centroids`) for
#'   downstream use (e.g. adding heatmaps and links).
#'
#' @param fill.by,color.by,shape,pointsize,pointalpha,pointstroke,fill,color,jitter,jitter_sd,pointlabel,pointlabelsize,nodelabsize,ring_n,plot_line,linecolor,mapping_line,linealpha,curve,curvature,label,labelsize,labelsegmentsize,labelsegmentalpha,label_layout,label_wrap_width,label_outer_pad,add_outer,q_outer,expand_outer,bandwidth_scale,outerwidth,outerlinetype,outeralpha,add_group_outer,add_group_outer_expand,add_group_outer_color,add_group_outer_fill,add_group_outer_fill_alpha,add_group_outer_linetype,add_group_outer_linewidth,layout.module,group.by,remove,dropOthers
#'   \lifecycle{deprecated} Old argument names kept for backward
#'   compatibility; see the *Deprecated arguments* section for the new name
#'   of each one. Using them emits a deprecation warning and the value is
#'   forwarded to the new argument.
#'
#' @returns A ggplot object, or when `return_layout = TRUE`, a list with
#'   `$plot` and `$layout_data`.
#' @md
#' @export
#'
#' @examples
#' \donttest{
#' library(ggNetView)
#' data("otu_rare_relative")
#' data("tax_tab")
#'
#' graph_obj <- build_graph_from_mat(
#'   mat = otu_rare_relative, transfrom.method = "none",
#'   r.threshold = 0.7, p.threshold = 0.05, method = "WGCNA",
#'   cor.method = "pearson", proc = "bonferroni",
#'   module.method = "Fast_greedy", node_annotation = tax_tab,
#'   top_modules = 15, seed = 1115
#' )
#'
#' # tidyverse-style arguments (>= 0.2.0)
#' ggNetView(graph_obj, layout = "gephi",
#'           node_fill = "Modularity", node_shape = 21,
#'           edge_color = "corr_direction", edge_width = "weight",
#'           module_label = TRUE)
#'
#' # solid shapes (no fill slot): the fill mapping is applied to colour
#' ggNetView(graph_obj, layout = "gephi", node_shape = 16, node_fill = "Phylum")
#' }
ggNetView <- function(graph_obj,
                      layout = NULL,
                      # ---- layout geometry -------------------------------
                      node_add = 7,
                      r = 1,
                      center = TRUE,
                      idx = NULL,
                      shrink = 1,
                      inner_shrink = 1,
                      k_nn = 12,
                      push_others_delta = 0,
                      layout_module = c("random", "adjacent", "order"),
                      group_by = "Modularity",
                      # ---- nodes ------------------------------------------
                      node_fill = "Modularity",
                      node_color = NULL,
                      node_shape = 21,
                      node_size = "Degree",
                      node_size_range = c(1, 10),
                      node_alpha = 1,
                      node_stroke = 0.3,
                      node_fill_values = NULL,
                      node_color_values = NULL,
                      node_jitter = FALSE,
                      node_jitter_sd = 0.1,
                      node_label = NULL,
                      node_label_size = 5,
                      # ---- edges ------------------------------------------
                      show_edges = TRUE,
                      edge_color = "grey70",
                      edge_color_values = NULL,
                      edge_width = 0.5,
                      edge_width_range = c(0.2, 1.5),
                      edge_linetype = 1,
                      edge_alpha = 0.25,
                      edge_curve = FALSE,
                      edge_curvature = 0.25,
                      # ---- module labels ----------------------------------
                      module_label = FALSE,
                      module_label_size = 10,
                      module_label_segment_width = 1,
                      module_label_segment_alpha = 1,
                      module_label_layout = c("two_column", "two_column_follow", "label_circle"),
                      module_label_wrap = NULL,
                      module_label_pad = 0.25,
                      # ---- module outline ---------------------------------
                      module_outline = FALSE,
                      module_outline_q = 0.88,
                      module_outline_expand = 1.02,
                      module_outline_bandwidth = 2,
                      module_outline_width = 1,
                      module_outline_linetype = 1,
                      module_outline_alpha = 0.5,
                      # ---- network outline --------------------------------
                      network_outline = FALSE,
                      network_outline_expand = 2,
                      network_outline_color = "grey50",
                      network_outline_fill = NULL,
                      network_outline_fill_alpha = 0.2,
                      network_outline_linetype = 1,
                      network_outline_width = 0.5,
                      # ---- misc -------------------------------------------
                      hide_others = FALSE,
                      drop_others = FALSE,
                      orientation = "up",
                      angle = 0,
                      scale = TRUE,
                      anchor_dist = 6,
                      nrow = NULL,
                      ncol = NULL,
                      seed = 1115,
                      scale_radius = NULL,
                      return_layout = FALSE,
                      # ---- deprecated (< 0.2.0) names ---------------------
                      fill.by = deprecated(),
                      color.by = deprecated(),
                      shape = deprecated(),
                      pointsize = deprecated(),
                      pointalpha = deprecated(),
                      pointstroke = deprecated(),
                      fill = deprecated(),
                      color = deprecated(),
                      jitter = deprecated(),
                      jitter_sd = deprecated(),
                      pointlabel = deprecated(),
                      pointlabelsize = deprecated(),
                      nodelabsize = deprecated(),
                      ring_n = deprecated(),
                      plot_line = deprecated(),
                      linecolor = deprecated(),
                      mapping_line = deprecated(),
                      linealpha = deprecated(),
                      curve = deprecated(),
                      curvature = deprecated(),
                      label = deprecated(),
                      labelsize = deprecated(),
                      labelsegmentsize = deprecated(),
                      labelsegmentalpha = deprecated(),
                      label_layout = deprecated(),
                      label_wrap_width = deprecated(),
                      label_outer_pad = deprecated(),
                      add_outer = deprecated(),
                      q_outer = deprecated(),
                      expand_outer = deprecated(),
                      bandwidth_scale = deprecated(),
                      outerwidth = deprecated(),
                      outerlinetype = deprecated(),
                      outeralpha = deprecated(),
                      add_group_outer = deprecated(),
                      add_group_outer_expand = deprecated(),
                      add_group_outer_color = deprecated(),
                      add_group_outer_fill = deprecated(),
                      add_group_outer_fill_alpha = deprecated(),
                      add_group_outer_linetype = deprecated(),
                      add_group_outer_linewidth = deprecated(),
                      layout.module = deprecated(),
                      group.by = deprecated(),
                      remove = deprecated(),
                      dropOthers = deprecated()
                      ){

  # ---- lifecycle: translate deprecated (< 0.2.0) argument names ----------
  .fn_env <- environment()
  .old_args <- .ggnv_collect_deprecated(.fn_env)
  if (length(.old_args) > 0L) {
    .explicit_new <- setdiff(names(match.call(expand.dots = FALSE))[-1L],
                             names(.ggnv_deprecated_arg_map))
    .renamed <- .ggnv_rename_args(.old_args, fn = "ggNetView",
                                  env = .fn_env, user_env = parent.frame())
    for (.nm in names(.renamed)) {
      if (.nm %in% .explicit_new) next          # explicit new name wins
      assign(.nm, .renamed[[.nm]], envir = .fn_env)
    }
  }

  layout_module <- match.arg(layout_module)
  module_label_layout <- match.arg(module_label_layout)

  if (!is.null(module_label_wrap)) {
    if (!is.numeric(module_label_wrap) || length(module_label_wrap) != 1 ||
        is.na(module_label_wrap) || module_label_wrap < 1) {
      stop("`module_label_wrap` must be NULL or a single positive integer.")
    }
    module_label_wrap <- as.integer(module_label_wrap)
  }

  if (!is.numeric(module_label_pad) || length(module_label_pad) != 1 ||
      is.na(module_label_pad) || module_label_pad < 0) {
    stop("`module_label_pad` must be a single non-negative numeric value.")
  }

  if (is.null(node_size_range)) node_size_range <- c(1, 10)
  if (is.null(edge_width_range)) edge_width_range <- c(0.2, 1.5)
  if (!is.numeric(node_size_range) || length(node_size_range) != 2L ||
      anyNA(node_size_range)) {
    stop("`node_size_range` must be a numeric vector of length 2.")
  }
  if (!is.numeric(edge_width_range) || length(edge_width_range) != 2L ||
      anyNA(edge_width_range)) {
    stop("`edge_width_range` must be a numeric vector of length 2.")
  }

  .ggnv_local_seed(seed)

  # drop_others acts on the source graph_obj BEFORE layout:
  # it removes "Others" nodes first, then downstream layout/plot are rebuilt.
  if (isTRUE(drop_others)) {
    node_tbl <- graph_obj %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()

    module_candidates <- c("Modularity", "modularity3", "modularity2")
    module_col <- module_candidates[module_candidates %in% colnames(node_tbl)]
    module_col <- if (length(module_col) > 0) module_col[[1]] else NULL

    if (!is.null(module_col)) {
      if ("name" %in% colnames(node_tbl)) {
        keep_names <- node_tbl %>%
          dplyr::filter(as.character(.data[[module_col]]) != "Others") %>%
          dplyr::pull(name)

        graph_obj <- graph_obj %>%
          tidygraph::activate(nodes) %>%
          tidygraph::filter(name %in% keep_names)
      } else {
        graph_obj <- graph_obj %>%
          tidygraph::activate(nodes) %>%
          tidygraph::filter(as.character(.data[[module_col]]) != "Others")
      }
    } else {
      warning("`drop_others = TRUE` but no module column found in `graph_obj` nodes.")
    }
  }

  if (is.logical(module_label)) {
    if (length(module_label) != 1 || is.na(module_label)) {
      stop("`module_label` must be a single logical or character string.")
    }
    show_module_label <- isTRUE(module_label)
    module_label_prefix <- "Modularity"
  } else if (is.character(module_label)) {
    if (length(module_label) != 1 || is.na(module_label)) {
      stop("`module_label` must be a single logical or character string.")
    }
    module_label_prefix <- trimws(module_label)
    if (identical(module_label_prefix, "")) {
      module_label_prefix <- "Modularity"
    }
    show_module_label <- TRUE
  } else {
    stop("`module_label` must be a single logical or character string.")
  }

  module_label_fun <- function(x) {
    x_chr <- as.character(x)
    ifelse(x_chr == "Others", "Others", paste0(module_label_prefix, x_chr))
  }

  is_module_field <- function(var_name) {
    is.character(var_name) &&
      length(var_name) == 1 &&
      !is.na(var_name) &&
      tolower(var_name) %in% c("modularity", "modularity2", "modularity3")
  }

  # Point-level aesthetics (`node_shape`, `node_fill`) should reflect raw node
  # attributes by default, but if mapping a modularity field, use the same
  # prefix style as module labels regardless of `module_label` visibility.
  point_legend_label_fun <- function(x, var_name) {
    if (is_module_field(var_name)) {
      module_label_fun(x)
    } else {
      as.character(x)
    }
  }

  # Positive / negative edge counts go into the title when edges are coloured
  # by correlation sign.
  edge_tbl_names <- graph_obj %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    colnames()
  show_sign_stats <- identical(edge_color, "corr_direction") &&
    "corr_direction" %in% edge_tbl_names
  stat_graph <- stat_graph(graph_obj, show_sign_stats)



  # validate `layout` against the available create_layout_* functions before
  # dispatching, so an unknown/misspelled layout gives a friendly error listing
  # the valid options instead of an obscure getFromNamespace() failure.
  if (is.null(layout) || !is.character(layout) || length(layout) != 1L) {
    stop("`layout` must be a single character string naming a layout.",
         call. = FALSE)
  }
  available_layouts <- sub("^create_layout_", "",
                           grep("^create_layout_",
                                ls(getNamespace("ggNetView")),
                                value = TRUE))
  if (!(layout %in% available_layouts)) {
    stop("Unknown `layout = \"", layout, "\"`.\n",
         "  Available layouts: ",
         paste(sort(available_layouts), collapse = ", "),
         call. = FALSE)
  }

  # find layout function
  func_name <- paste0("create_layout_", layout)

  # find layout functions from ggNetView package
  lay_func <- utils::getFromNamespace(func_name, "ggNetView")

  # get ly1
  if (layout == "consensus_module_equal_gephi" | layout == "consensus_module_gephi") {
    ly1 = lay_func(graph_obj = graph_obj,
                   node_add = node_add,
                   r = r,
                   scale = scale,
                   anchor_dist = anchor_dist,
                   orientation = orientation,
                   angle = angle,
                   nrow = nrow,
                   ncol = ncol)
  }else if (layout == "WGCNA") {
    # WGCNA layout has an extra `inner_shrink` parameter that controls
    # intra-module compactness independently of inter-module spacing.
    # Other layouts do not accept this argument, so we dispatch it here
    # only.  Default `inner_shrink = 1` reproduces historical behaviour.
    ly1 = lay_func(graph_obj = graph_obj,
                   node_add = node_add,
                   r = r,
                   inner_shrink = inner_shrink,
                   scale = scale,
                   anchor_dist = anchor_dist,
                   orientation = orientation,
                   angle = angle)
  }else{
    ly1 = lay_func(graph_obj = graph_obj,
                   node_add = node_add,
                   r = r,
                   scale = scale,
                   anchor_dist = anchor_dist,
                   orientation = orientation,
                   angle = angle)
  }



  # get ly1_1

  ly1_1 <- NULL

  # `circlepack` is a self-arranging layout: it already packs modules and places
  # nodes, so it bypasses the module_layout*() arrangers (which would re-arrange
  # the modules) and goes through a no-op passthrough that just wraps the
  # coordinates. Guarding the branches below with `is.null(ly1_1)` leaves every
  # other layout's dispatch byte-for-byte unchanged.
  if (func_name == "create_layout_circlepack") {
    ly1_1 <- module_layout_passthrough(graph_obj,
                                       layout = ly1,
                                       jitter = node_jitter,
                                       jitter_sd = node_jitter_sd)
  }

  if (is.null(ly1_1) && layout_module == "random") {
    ly1_1 <- module_layout(graph_obj,
                           layout = ly1,
                           center = center,
                           idx = idx,
                           shrink = shrink,
                           jitter = node_jitter,
                           jitter_sd = node_jitter_sd# ,
                           # seed = seed
    )
  }

  if (is.null(ly1_1) && layout_module == "adjacent") {
    k_nn_try <- k_nn
    k_nn_cap <- max(1, nrow(ly1) - 1)
    ly1_1 <- NULL
    while (is.null(ly1_1)) {
      ly_try <- tryCatch(
        module_layout3(graph_obj,
                       layout = ly1,
                       center = center,
                       k_nn = k_nn_try,
                       push_others_delta = push_others_delta,
                       shrink = shrink,
                       jitter = node_jitter,
                       jitter_sd = node_jitter_sd
                       # seed = seed
        ),
        error = function(e) e
      )

      if (!inherits(ly_try, "error")) {
        ly1_1 <- ly_try
        break
      }

      err_msg <- conditionMessage(ly_try)
      is_slot_error <- grepl("consecutive slots", err_msg, ignore.case = TRUE)
      if (!is_slot_error || k_nn_try >= k_nn_cap) {
        stop(ly_try)
      }

      next_k <- min(k_nn_cap, max(k_nn_try + 20, ceiling(k_nn_try * 1.25)))
      if (next_k <= k_nn_try) {
        stop(ly_try)
      }
      warning(sprintf(
        "`layout_module = 'adjacent'` failed at k_nn = %d; retrying with k_nn = %d.",
        k_nn_try, next_k
      ))
      k_nn_try <- next_k
    }
  }

  if (is.null(ly1_1) && layout_module == "order" & func_name != "create_layout_multirings") {
    ly1_1 <- module_layout4(graph_obj,
                            layout = ly1,
                            center = center,
                            k_nn = k_nn,
                            push_others_delta = push_others_delta,
                            shrink = shrink,
                            jitter = node_jitter,
                            jitter_sd = node_jitter_sd
                            # seed = seed
    )
  }

  if (is.null(ly1_1) && layout_module == "order" & func_name == "create_layout_multirings") {
    ly1_1 <- module_layout5(graph_obj,
                            layout = ly1,
                            center = center,
                            k_nn = k_nn,
                            push_others_delta = push_others_delta,
                            shrink = shrink,
                            jitter = node_jitter,
                            jitter_sd = node_jitter_sd
                            # seed = seed
    )
  }

  # Normal layout


  if (group_by != "pie") {

    # Optional: scale layout to fit in radius (for use with ggnetview_modularity_heatmaps)
    if (!is.null(scale_radius) && is.finite(scale_radius) && scale_radius > 0) {
      xr_net <- range(ly1_1$graph_ly_final$x, na.rm = TRUE)
      yr_net <- range(ly1_1$graph_ly_final$y, na.rm = TRUE)
      scale_net <- max(diff(xr_net), diff(yr_net), 1e-8)
      cx <- mean(xr_net)
      cy <- mean(yr_net)
      ly1_1[["graph_ly_final"]] <- ly1_1[["graph_ly_final"]] %>%
        dplyr::mutate(
          x = (x - cx) / scale_net * scale_radius,
          y = (y - cy) / scale_net * scale_radius
        )
      ly1_1[["layout"]] <- dplyr::select(ly1_1[["graph_ly_final"]], x, y)
      ly1_1[["ggplot_data"]] <- get_location(ly1_1[["graph_ly_final"]], ly1_1[["graph_obj"]])
    }



    module_number <- ly1_1$graph_ly_final$Modularity %>% as.character() %>% unique() %>% length()

    # module info
    if (module_number == 1) {
      module_info <-  ly1_1$graph_ly_final$Modularity %>% as.character() %>% unique()

      ly1_1[["graph_ly_final"]] <- ly1_1[["graph_ly_final"]] %>%
        dplyr::mutate(Modularity = as.character(Modularity)) %>%
        dplyr::mutate(Modularity = factor(Modularity))

      ly1_1[["graph_obj"]] <- ly1_1[["graph_obj"]] %>%
        tidygraph::mutate(Modularity = as.character(Modularity)) %>%
        tidygraph::mutate(Modularity = factor(Modularity, ordered = TRUE))

    }else{
      module_info <- levels(ly1_1$graph_ly_final$Modularity)
      module_info <- module_info[module_info!="Others"]
    }

    # hide_others acts only at the visualization stage (post-layout):
    # it drops "Others" from rendered data while keeping the existing layout.
    if (isFALSE(hide_others)) {
      ly1_1 <- ly1_1
    }else{
      ly1_1[["graph_ly_final"]] <- ly1_1[["graph_ly_final"]] %>%
        dplyr::filter(as.character(Modularity) %in% module_info)

      ly1_1[["graph_obj"]] <- ly1_1[["graph_obj"]] %>%
        tidygraph::filter(name %in% ly1_1[["graph_ly_final"]]$name)


      ly1_1[["layout"]] <- dplyr::select(ly1_1[["graph_ly_final"]], x, y)

      ly1_1[["ggplot_data"]] <- get_location(ly1_1[["graph_ly_final"]],
                                             ly1_1[["graph_obj"]])
    }

    xr <- NULL; yr <- NULL; x_mid <- NULL
    dx <- NULL; pad <- NULL; lab_df <- NULL
    plot_xlim <- NULL; plot_ylim <- NULL; label_force <- NULL
    # Non-NULL -> render code adds a geom_segment first-leg from module
    # centroid to the elbow on the network boundary (used by
    # two_column_follow). NULL -> no manual first leg.
    lab_leader_df <- NULL
    # TRUE -> pass segment.square = TRUE to geom_text_repel so ggrepel
    # draws an L-shape leader. Requires ggrepel >= 0.9.4. Used by
    # label_circle. Combined with segment.squareShape = 0 to force the
    # elbow at the true L corner instead of collapsing to the label end.
    label_segment_square       <- FALSE
    label_segment_square_shape <- 1
    # point.padding passed to geom_text_repel. 0 closes any gap at the
    # aes point -- needed by two_column_follow whose aes point IS the
    # elbow shared with the manual first leg (any padding leaves a
    # visible disconnect at the elbow).
    label_point_padding        <- 0.15

    .build_label_location <- function(){
      # ---- basic geometry shared by both strategies ----
      xr <<- range(ly1_1[["layout"]]$x)
      yr <<- range(ly1_1[["layout"]]$y)
      x_mid <<- stats::median(ly1_1[["layout"]]$x)
      dx <<- diff(xr) * module_label_pad
      pad <<- dx * 1.2

      # text-wrap helper (no-op when module_label_wrap is NULL)
      .wrap_label <- function(txt) {
        if (is.null(module_label_wrap)) {
          as.character(txt)
        } else {
          stringr::str_wrap(as.character(txt), width = module_label_wrap)
        }
      }

      # Partial equispacing: enforce a minimum angular gap between
      # neighbouring labels (analog of p2's rank-rescaling within a
      # side, but only "spreading where needed"). Sparse modules keep
      # theta_target = theta_actual (zero displacement); clustered
      # modules get pushed apart just enough to clear min_gap. The
      # circle is "broken" at the largest natural gap so naturally
      # close-but-wrapped modules (near +/- pi) see each other.
      # Re-centring after the forward pass minimises overall drift so
      # the cluster's mean angle is preserved.
      #   min_gap = 2*pi / (2*N) = half the average per-label slice.
      .partial_equispace <- function(theta_in) {
        n <- length(theta_in)
        if (n < 2) return(theta_in)
        sort_order  <- order(theta_in)
        sorted_th   <- theta_in[sort_order]
        gaps        <- c(diff(sorted_th),
                         sorted_th[1] + 2 * pi - sorted_th[n])
        max_gap_idx <- which.max(gaps)
        if (max_gap_idx < n) {
          rotation     <- c((max_gap_idx + 1):n, seq_len(max_gap_idx))
          rotated_th   <- sorted_th[rotation]
          wrap_n       <- max_gap_idx
          rotated_th[(n - wrap_n + 1):n] <-
            rotated_th[(n - wrap_n + 1):n] + 2 * pi
          rotated_orig <- sort_order[rotation]
        } else {
          rotated_th   <- sorted_th
          rotated_orig <- sort_order
        }
        min_gap  <- 2 * pi / (2 * n)
        adjusted <- rotated_th
        for (i in 2:n) {
          if (adjusted[i] - adjusted[i - 1] < min_gap) {
            adjusted[i] <- adjusted[i - 1] + min_gap
          }
        }
        # re-centre so the spread sequence's mean matches the originals'
        adjusted <- adjusted + (mean(rotated_th) - mean(adjusted))
        # normalise into [-pi, pi]
        adjusted <- atan2(sin(adjusted), cos(adjusted))
        out <- numeric(n)
        out[rotated_orig] <- adjusted
        out
      }

      # shared base: one row per module (excluding "Others"), with side / y_rank
      base_df <- ly1_1[["graph_ly_final"]] %>%
        dplyr::distinct(modularity3, .keep_all = TRUE) %>%
        dplyr::filter(modularity3 != "Others") %>%
        dplyr::mutate(side = ifelse(x < x_mid, "left", "right")) %>%
        dplyr::group_by(side) %>%
        dplyr::arrange(y, .by_group = TRUE) %>%
        dplyr::mutate(
          y_rank   = dplyr::row_number(),
          y_target = scales::rescale(y_rank, to = yr)
        ) %>%
        dplyr::ungroup()

      if (identical(module_label_layout, "two_column")) {
        # ===== two fixed columns: x_anchor pinned to xr[1] - dx / xr[2] + dx
        lab_df <<- base_df %>%
          dplyr::mutate(
            x_anchor = dplyr::if_else(side == "left", xr[1] - dx, xr[2] + dx),
            y_anchor = y_target,
            nudge_x  = x_anchor - x,
            nudge_y  = y_target - y,
            hjust    = dplyr::if_else(side == "left", 1, 0),
            vjust    = 0.5,
            .label_text = .wrap_label(module_label_fun(modularity3))
          )

        plot_xlim                  <<- c(xr[1] - pad, xr[2] + pad)
        plot_ylim                  <<- yr
        label_force                <<- 0.05
        lab_leader_df              <<- NULL    # no manual first leg
        label_segment_square       <<- FALSE   # straight ggrepel segment
        label_segment_square_shape <<- 1
        label_point_padding        <<- 0.15    # default gap from node
      } else if (identical(module_label_layout, "two_column_follow")) {
        # ===== two_column_follow: label at module's actual angle ==========
        # No equispacing -- each label sits at its module's REAL angular
        # position on the outer ellipse, so left-side modules get
        # left-side labels and right-side modules get right-side labels.
        # Leader has two segments:
        #   - first leg (manual geom_segment in the render code) from
        #     module centroid to an elbow on the network's bounding
        #     ellipse at the same theta_actual -- this is the radial
        #     "exit the network" stub
        #   - second leg drawn by ggrepel from label back to the elbow.
        # ggrepel force = 1 lets it slide overlapping labels (along the
        # tangent direction in practice) without breaking the leader,
        # because ggrepel always retargets its segment to the elbow.
        centroids <- ly1_1[["graph_ly_final"]] %>%
          dplyr::filter(modularity3 != "Others") %>%
          dplyr::group_by(modularity3) %>%
          dplyr::summarise(mx = mean(x), my = mean(y), .groups = "drop")

        cx <- mean(xr)
        cy <- mean(yr)
        R_x_net <- (xr[2] - xr[1]) / 2
        R_y_net <- (yr[2] - yr[1]) / 2
        if (!is.finite(R_x_net) || R_x_net <= 0) R_x_net <- 1
        if (!is.finite(R_y_net) || R_y_net <= 0) R_y_net <- 1
        R_x_outer <- R_x_net * (1 + module_label_pad)
        R_y_outer <- R_y_net * (1 + module_label_pad)

        base_follow <- ly1_1[["graph_ly_final"]] %>%
          dplyr::distinct(modularity3, .keep_all = TRUE) %>%
          dplyr::filter(modularity3 != "Others") %>%
          dplyr::left_join(centroids, by = "modularity3") %>%
          dplyr::mutate(theta_actual = atan2(my - cy, mx - cx))

        # Partial equispacing of theta_target (closer to p2's intent):
        # sparse modules keep theta_target = theta_actual so their
        # labels sit exactly above/right/etc. of their modules;
        # clustered modules get pushed apart just enough to clear the
        # min angular gap so labels don't overlap. The L bend will be
        # visible only for clustered modules where theta_target
        # actually differs from theta_actual -- which is exactly where
        # an L is needed to distinguish neighbouring labels.
        theta_target_ord <- .partial_equispace(base_follow$theta_actual)

        lab_df <<- base_follow %>%
          dplyr::mutate(
            theta_target = theta_target_ord,
            # label sits on the outer ellipse at the equispaced
            # theta_target
            x_anchor = cx + R_x_outer * cos(theta_target),
            y_anchor = cy + R_y_outer * sin(theta_target),
            # elbow on the network's bounding ellipse at the module's
            # ACTUAL angle: first leg of the L goes straight radially
            # from the module centroid out to the network boundary.
            # Because theta_target != theta_actual (thanks to
            # equispacing), the second leg ggrepel draws -- from the
            # label back to this elbow -- is slanted, giving a clear L.
            elbow_x  = cx + R_x_net * cos(theta_actual),
            elbow_y  = cy + R_y_net * sin(theta_actual),
            # ggrepel's "point" = the elbow. ggrepel will draw the
            # second leg from the label back to this elbow.
            x        = elbow_x,
            y        = elbow_y,
            nudge_x  = x_anchor - elbow_x,
            nudge_y  = y_anchor - elbow_y,
            # hjust / vjust fan the text outward along theta_target
            hjust    = (1 - cos(theta_target)) / 2,
            vjust    = (1 - sin(theta_target)) / 2,
            .label_text = .wrap_label(module_label_fun(modularity3))
          )

        # lab_leader_df is just a flag (re-using lab_df) so the render
        # code knows to draw the manual first leg via geom_segment from
        # the module centroid (mx, my) to the elbow (elbow_x, elbow_y).
        lab_leader_df <<- lab_df

        # plot_xlim matches two_column (left/right labels fit in the
        # 0.6 R_x_net side pad). plot_ylim must extend to the outer
        # ellipse so top/bottom label anchors don't sit outside the
        # plot panel -- two_column doesn't need this because its
        # labels never leave yr in the y direction. The theme's 20pt
        # plot.margin covers the remaining text-height overflow above
        # / below each anchor.
        plot_xlim   <<- c(xr[1] - pad, xr[2] + pad)
        plot_ylim   <<- c(cy - R_y_outer, cy + R_y_outer)
        # force = 0: labels stay exactly at the outer ellipse position
        # we computed. Using force = 1 here can push labels back INTO
        # the network when ggrepel resolves label-label overlaps, which
        # the user explicitly forbids. If labels overlap because too
        # many modules sit at similar angles, increase module_label_pad
        # so the outer ring has more tangential room.
        label_force                <<- 0
        label_segment_square       <<- FALSE   # ggrepel draws single line
        label_segment_square_shape <<- 1
        # ★ no padding around aes point: the aes IS the elbow shared
        # with the manual first leg, so ggrepel's segment must reach
        # it exactly or you see a visible gap at the L corner
        label_point_padding        <<- 0
      } else {
        # ===== label_circle: pure ggrepel + module's actual angle =========
        # Same "label at module's actual angle on the outer ellipse" as
        # two_column_follow, but the leader is drawn entirely by ggrepel
        # (no manual first leg). ggrepel's aes(x, y) is the module
        # centroid, nudge places the label initially at outer ellipse at
        # theta_actual, force = 1 lets ggrepel spread overlapping labels.
        # segment.square = TRUE (with squareShape = 0) makes ggrepel draw
        # an L-shape leader. Simpler code than two_column_follow but the
        # L bend can visually collapse when label and module sit on the
        # same radial. Requires ggrepel >= 0.9.4.
        centroids <- ly1_1[["graph_ly_final"]] %>%
          dplyr::filter(modularity3 != "Others") %>%
          dplyr::group_by(modularity3) %>%
          dplyr::summarise(mx = mean(x), my = mean(y), .groups = "drop")

        cx <- mean(xr)
        cy <- mean(yr)
        R_x_net <- (xr[2] - xr[1]) / 2
        R_y_net <- (yr[2] - yr[1]) / 2
        if (!is.finite(R_x_net) || R_x_net <= 0) R_x_net <- 1
        if (!is.finite(R_y_net) || R_y_net <= 0) R_y_net <- 1
        R_x_outer <- R_x_net * (1 + module_label_pad)
        R_y_outer <- R_y_net * (1 + module_label_pad)

        base_follow <- ly1_1[["graph_ly_final"]] %>%
          dplyr::distinct(modularity3, .keep_all = TRUE) %>%
          dplyr::filter(modularity3 != "Others") %>%
          dplyr::left_join(centroids, by = "modularity3") %>%
          dplyr::mutate(theta_actual = atan2(my - cy, mx - cx))

        # Partial equispacing: sparse modules keep theta_target =
        # theta_actual, clustered modules get pushed apart just
        # enough to clear the min angular gap. See two_column_follow
        # for the algorithm rationale.
        theta_target_ord <- .partial_equispace(base_follow$theta_actual)

        lab_df <<- base_follow %>%
          dplyr::mutate(
            theta_target = theta_target_ord,
            x_anchor     = cx + R_x_outer * cos(theta_target),
            y_anchor     = cy + R_y_outer * sin(theta_target),
            # ggrepel's "point" = module centroid; nudge places the
            # label on the outer ellipse at the equispaced theta_target
            x        = mx,
            y        = my,
            nudge_x  = x_anchor - mx,
            nudge_y  = y_anchor - my,
            hjust    = (1 - cos(theta_target)) / 2,
            vjust    = (1 - sin(theta_target)) / 2,
            .label_text = .wrap_label(module_label_fun(modularity3))
          )

        # plot_xlim matches two_column; plot_ylim must extend to the
        # outer ellipse so top/bottom label anchors fit inside the
        # panel (theme's 20pt margin covers the remaining text height).
        plot_xlim   <<- c(xr[1] - pad, xr[2] + pad)
        plot_ylim   <<- c(cy - R_y_outer, cy + R_y_outer)
        lab_leader_df              <<- NULL    # ggrepel draws everything
        # force = 0: same reason as two_column_follow -- ggrepel's
        # collision avoidance would otherwise pull labels back toward
        # their aes point (the module centroid, inside the network).
        # Keep labels exactly at the outer ellipse; if they overlap
        # because of clustered angles, bump label_outer_pad.
        label_force                <<- 0
        # segment.square = FALSE so ggrepel draws a single straight
        # slanted line from label to module. With segment.square = TRUE,
        # the L would visually collapse here because the label sits on
        # the module's radial line (label, network centre, module are
        # collinear -- the bend has no perpendicular component to show).
        # Using a single straight line makes label_circle visually
        # distinct from two_column_follow's right-angle L.
        label_segment_square       <<- FALSE
        label_segment_square_shape <<- 1       # unused when square=FALSE
        label_point_padding        <<- 0.15    # default gap from module
      }
    }



    # outlier df location
    .build_mask_table <- function(){
      maskTable <- generateMask_ggnetview(
        dims = ly1_1[["layout"]],
        clusters = ly1_1[["graph_obj"]] %>%
          tidygraph::activate(nodes) %>%
          tidygraph::as_tibble() %>%
          dplyr::pull(modularity3),
        q = module_outline_q,
        expand = module_outline_expand,
        bandwidth_scale = module_outline_bandwidth
      )

      return(maskTable)
    }

    ####----Plot----####
    # base plot
    p1_1 <- ggplot2::ggplot()

    node_df <- ly1_1[["ggplot_data"]][[1]]
    edge_df <- ly1_1[["ggplot_data"]][[2]]

    # "column name or literal value" resolver shared by node_* / edge_* args:
    # a single string matching a column of `df` is an aesthetic mapping,
    # anything else is a constant.
    .is_mapped <- function(x, df) {
      is.character(x) && length(x) == 1L && !is.na(x) && x %in% colnames(df)
    }

    # ---- network outline (whole-network circle) -----------------------------
    if (isTRUE(network_outline) && nrow(node_df) > 0) {
      group_circle_df <- node_df %>%
        dplyr::mutate(.group_outer = 1L)
      circle_n_grp <- max(40, min(300, as.integer(round(8 * sqrt(nrow(group_circle_df))))))
      fill_grp <- if (is.null(network_outline_fill) || length(network_outline_fill) == 0L) NA else network_outline_fill[1L]
      alpha_grp <- if (is.na(fill_grp)) 1 else network_outline_fill_alpha
      p1_1 <- p1_1 +
        ggforce::geom_mark_circle(
          data = group_circle_df,
          mapping = ggplot2::aes(x = x, y = y, group = .group_outer),
          fill = fill_grp,
          alpha = alpha_grp,
          color = network_outline_color,
          linetype = network_outline_linetype,
          linewidth = network_outline_width,
          n = circle_n_grp,
          expand = grid::unit(network_outline_expand, "mm")
        )
    }

    # ---- edges ----------------------------------------------------------------
    edge_color_mapped    <- .is_mapped(edge_color, edge_df)
    edge_width_mapped    <- .is_mapped(edge_width, edge_df)
    edge_linetype_mapped <- .is_mapped(edge_linetype, edge_df)

    if (is.character(edge_color) && length(edge_color) != 1L) {
      stop("`edge_color` must be a single colour or a single edge column name.")
    }
    if (is.character(edge_width) && !edge_width_mapped) {
      stop("`edge_width` must be a single number or the name of a numeric edge column.")
    }
    if (edge_width_mapped && !is.numeric(edge_df[[edge_width]])) {
      stop("`edge_width = \"", edge_width, "\"` must refer to a numeric edge column.")
    }

    if (isTRUE(show_edges) && nrow(edge_df) > 0) {
      edge_aes <- list(x = quote(from_x), xend = quote(to_x),
                       y = quote(from_y), yend = quote(to_y))
      edge_params <- list(data = edge_df, alpha = edge_alpha)

      # colour
      edge_color_scale <- NULL
      if (edge_color_mapped) {
        edge_aes$colour <- rlang::sym(edge_color)
        edge_vals <- edge_df[[edge_color]]
        if (is.numeric(edge_vals)) {
          edge_color_scale <- ggplot2::scale_color_gradient(
            low = "#4393c3", high = "#d6604d", name = edge_color)
        } else if (!is.null(edge_color_values)) {
          edge_color_scale <- ggplot2::scale_color_manual(
            values = edge_color_values, name = edge_color)
        } else if (identical(edge_color, "corr_direction")) {
          edge_color_scale <- ggplot2::scale_color_manual(
            values = c("Positive" = "#d6604d", "Negative" = "#4393c3"),
            name = edge_color)
        } else {
          edge_color_scale <- scale_color_ggnetview(
            .ggnv_class_order(edge_vals), name = edge_color)
        }
      } else {
        edge_params$colour <- edge_color
      }

      # width
      edge_width_scale <- NULL
      if (edge_width_mapped) {
        edge_aes$linewidth <- rlang::sym(edge_width)
        edge_width_scale <- ggplot2::scale_linewidth(
          range = edge_width_range, name = edge_width,
          guide = ggplot2::guide_legend(ncol = 1, order = 5))
      } else {
        edge_params$linewidth <- edge_width
      }

      # linetype
      edge_linetype_scale <- NULL
      if (edge_linetype_mapped) {
        edge_aes$linetype <- rlang::sym(edge_linetype)
        edge_linetype_scale <- ggplot2::scale_linetype(
          name = edge_linetype,
          guide = ggplot2::guide_legend(ncol = 1, order = 6))
      } else if (!is.null(edge_linetype)) {
        edge_params$linetype <- edge_linetype
      }

      edge_geom <- if (isTRUE(edge_curve)) ggplot2::geom_curve else ggplot2::geom_segment
      if (isTRUE(edge_curve)) edge_params$curvature <- edge_curvature

      edge_layer <- do.call(edge_geom,
                            c(list(mapping = ggplot2::aes(!!!edge_aes)), edge_params))

      p1_1 <- p1_1 +
        edge_layer +
        edge_color_scale +
        edge_width_scale +
        edge_linetype_scale +
        theme_ggnetview()
      # Edge colour and node colour are independent scales.
      if (edge_color_mapped) {
        p1_1 <- p1_1 + ggnewscale::new_scale_color()
      }
    }

    # ---- nodes ----------------------------------------------------------------
    node_fill_mapped  <- .is_mapped(node_fill, node_df)
    node_color_mapped <- .is_mapped(node_color, node_df)
    node_shape_mapped <- .is_mapped(node_shape, node_df)
    node_size_mapped  <- .is_mapped(node_size, node_df)

    if (is.character(node_shape) && !node_shape_mapped) {
      stop("`node_shape` must be a shape code or the name of a node column.")
    }
    if (is.character(node_size) && !node_size_mapped) {
      stop("`node_size` must be a single number or the name of a numeric node column.")
    }
    if (node_size_mapped && !is.numeric(node_df[[node_size]])) {
      stop("`node_size = \"", node_size, "\"` must refer to a numeric node column.")
    }
    if (is.character(node_fill) && length(node_fill) != 1L) {
      stop("`node_fill` must be a single colour or a single node column name.")
    }
    if (!is.null(node_color) && is.character(node_color) && length(node_color) != 1L) {
      stop("`node_color` must be NULL, a single colour or a single node column name.")
    }

    # Shapes 0-20 have no fill slot: route the fill mapping to `colour`.
    shape_has_fill <- node_shape_mapped || (is.numeric(node_shape) && all(node_shape %in% 21:25))
    fill_aes <- "fill"
    if (!shape_has_fill && node_fill_mapped) {
      if (node_color_mapped) {
        warning("`node_shape = ", node_shape, "` has no fill slot; `node_fill = \"",
                node_fill, "\"` is ignored because `node_color` is also mapped. ",
                "Use a fillable shape (21-25) to map both.", call. = FALSE)
        node_fill_mapped <- FALSE
      } else {
        fill_aes <- "colour"
        if (!is.null(node_color)) node_color <- NULL   # mapping wins over constant
      }
    }

    point_label_df <- NULL
    point_label_col <- NULL
    if (!is.null(node_label)) {
      if (!is.character(node_label) || length(node_label) != 1) {
        stop("`node_label` must be NULL, 'ALL', or 'topN' (N is a positive integer, e.g. 'top7').")
      }
      node_label_clean <- toupper(trimws(node_label))
      is_all <- identical(node_label_clean, "ALL")
      is_top_n <- grepl("^TOP[1-9][0-9]*$", node_label_clean)
      if (!is_all && !is_top_n) {
        stop("`node_label` must be NULL, 'ALL', or 'topN' (N is a positive integer, e.g. 'top7').")
      }

      if (!"Degree" %in% colnames(node_df)) {
        stop("`Degree` column is required in the node table for `node_label`.")
      }

      group_col <- if (group_by %in% colnames(node_df)) {
        group_by
      } else if ("Modularity" %in% colnames(node_df)) {
        "Modularity"
      } else {
        stop("No valid module column found for `node_label` grouping.")
      }

      point_label_col <- if ("ID" %in% colnames(node_df)) {
        "ID"
      } else if ("name" %in% colnames(node_df)) {
        "name"
      } else {
        stop("`node_label` requires an `ID` or `name` column in the node table.")
      }

      if (is_all) {
        point_label_df <- node_df
      } else {
        top_n <- as.integer(sub("^TOP", "", node_label_clean))

        point_label_df <- node_df %>%
          dplyr::group_by(.data[[group_col]]) %>%
          dplyr::slice_max(order_by = Degree, n = top_n, with_ties = FALSE) %>%
          dplyr::ungroup()
      }
    }

    merge_point_legends <- node_shape_mapped && node_fill_mapped &&
      identical(node_shape, node_fill)
    same_fill_color_mapping <- node_color_mapped && node_fill_mapped &&
      identical(node_color, node_fill)

    # -- fill scale (or colour scale when routed) --
    fill_scale_points <- NULL
    if (node_fill_mapped) {
      fill_vals <- ly1_1[["graph_ly_final"]][[node_fill]]
      if (is.numeric(fill_vals)) {
        fill_scale_points <- if (fill_aes == "fill") {
          ggplot2::scale_fill_gradient(low = "#4393c3", high = "#d6604d", name = node_fill)
        } else {
          ggplot2::scale_color_gradient(low = "#4393c3", high = "#d6604d", name = node_fill)
        }
      } else {
        fill_classes <- .ggnv_class_order(fill_vals)
        fill_guide <- ggplot2::guide_legend(ncol = 1, order = 1)
        if (fill_aes == "fill") {
          fill_scale_points <- if (is.null(node_fill_values)) {
            scale_fill_ggnetview(fill_classes,
                                 breaks = fill_classes,
                                 labels = function(x) point_legend_label_fun(x, node_fill),
                                 guide = fill_guide)
          } else {
            ggplot2::scale_fill_manual(values = node_fill_values,
                                       breaks = fill_classes,
                                       labels = function(x) point_legend_label_fun(x, node_fill),
                                       guide = fill_guide)
          }
        } else {
          fill_scale_points <- if (is.null(node_fill_values)) {
            scale_color_ggnetview(fill_classes,
                                  breaks = fill_classes,
                                  labels = function(x) point_legend_label_fun(x, node_fill),
                                  guide = fill_guide)
          } else {
            ggplot2::scale_color_manual(values = node_fill_values,
                                        breaks = fill_classes,
                                        labels = function(x) point_legend_label_fun(x, node_fill),
                                        guide = fill_guide)
          }
        }
      }
    }

    # -- colour scale (node border) --
    color_scale_points <- NULL
    if (node_color_mapped) {
      color_values <- node_df[[node_color]]
      if (is.numeric(color_values)) {
        color_scale_points <- ggplot2::scale_color_gradient(
          low = "#4393c3",
          high = "#d6604d",
          name = node_color,
          guide = if (same_fill_color_mapping) "none" else "legend"
        )
      } else if (is.null(node_color_values)) {
        color_scale_points <- scale_color_ggnetview(
          .ggnv_class_order(color_values),
          labels = function(x) point_legend_label_fun(x, node_color),
          guide = if (same_fill_color_mapping) "none" else ggplot2::guide_legend(ncol = 1, order = 2)
        )
      } else {
        color_scale_points <- ggplot2::scale_color_manual(
          values = node_color_values,
          labels = function(x) point_legend_label_fun(x, node_color),
          guide = if (same_fill_color_mapping) "none" else ggplot2::guide_legend(ncol = 1, order = 2)
        )
      }
    }

    # -- shape scale --
    shape_scale_points <- NULL
    if (node_shape_mapped) {
      shape_classes <- .ggnv_class_order(ly1_1[["graph_ly_final"]][[node_shape]])
      shape_values <- rep(21:25, length.out = length(shape_classes))
      shape_scale_points <- ggplot2::scale_shape_manual(
        values = shape_values,
        breaks = shape_classes,
        labels = function(x) point_legend_label_fun(x, node_shape),
        guide = ggplot2::guide_legend(ncol = 1, order = if (merge_point_legends) 1 else 2)
      )
    }

    # -- size scale --
    size_scale_points <- NULL
    if (node_size_mapped) {
      size_scale_points <- ggplot2::scale_size(
        range = node_size_range, name = node_size,
        guide = ggplot2::guide_legend(ncol = 1, order = 3))
    }

    # -- legend key overrides --
    legend_shape <- if (node_shape_mapped) 21 else node_shape[1L]
    size_guide_points <- NULL
    if (node_size_mapped && isTRUE(node_stroke == 0)) {
      size_guide_points <- ggplot2::guides(
        size = ggplot2::guide_legend(
          ncol = 1,
          order = 3,
          override.aes = list(
            shape = legend_shape,
            fill = "grey70",
            colour = "grey70",
            stroke = 0.3
          )
        )
      )
    }
    fill_guide_points <- NULL
    if (node_fill_mapped && fill_aes == "fill" && !node_color_mapped && !merge_point_legends) {
      fill_guide_points <- ggplot2::guides(
        fill = ggplot2::guide_legend(
          ncol = 1,
          order = 1,
          override.aes = list(
            shape = legend_shape,
            colour = if (!is.null(node_color)) node_color else "grey40",
            stroke = node_stroke
          )
        )
      )
    }

    # -- point layer --
    pt_aes <- list(x = quote(x), y = quote(y))
    if (node_size_mapped)  pt_aes$size  <- rlang::sym(node_size)
    if (node_fill_mapped)  pt_aes[[fill_aes]] <- rlang::sym(node_fill)
    if (node_shape_mapped) pt_aes$shape <- rlang::sym(node_shape)
    if (node_color_mapped) pt_aes$colour <- rlang::sym(node_color)

    pt_params <- list(data = node_df, alpha = node_alpha, stroke = node_stroke)
    if (!node_shape_mapped) pt_params$shape <- node_shape
    if (!node_size_mapped)  pt_params$size  <- node_size
    if (!node_fill_mapped && !is.null(node_fill) && shape_has_fill) pt_params$fill <- node_fill
    if (!node_fill_mapped && !is.null(node_fill) && !shape_has_fill && is.null(node_color)) {
      pt_params$colour <- node_fill
    }
    if (!node_color_mapped && !is.null(node_color)) pt_params$colour <- node_color

    point_layer <- do.call(ggplot2::geom_point,
                           c(list(mapping = ggplot2::aes(!!!pt_aes)), pt_params))

    p1_1 <- p1_1 +
      point_layer +
      size_scale_points +
      # the module label / outline blocks below add their own coord_equal()
      (if (isFALSE(show_module_label) && isFALSE(module_outline)) ggplot2::coord_fixed() else NULL) +
      theme_ggnetview() +
      fill_scale_points +
      color_scale_points +
      shape_scale_points +
      size_guide_points +
      fill_guide_points

    # Module labels / outlines are coloured by module with the node-fill
    # palette when `node_fill` is a module field; otherwise the default palette.
    module_palette <- if (node_fill_mapped && is_module_field(node_fill)) node_fill_values else NULL

    # label = FALSE module_outline = FALSE
    if (isFALSE(show_module_label) & isFALSE(module_outline)) {
      p1_1 <- p1_1

    }

    # label = TRUE module_outline = FALSE
    if (isTRUE(show_module_label) & isFALSE(module_outline)) {

      .build_label_location()

      lab_classes <- .ggnv_class_order(lab_df$Modularity)
      color_scale_lab <- if (is.null(module_palette)) scale_color_ggnetview(lab_classes, labels = module_label_fun) else ggplot2::scale_color_manual(values = module_palette, labels = module_label_fun)

      p1_1 <- p1_1 +
        ggnewscale::new_scale_color() +
        # First leg of the two-segment leader (two_column_follow only):
        # from module centroid (mx, my) to the elbow on the network
        # boundary (elbow_x, elbow_y). The second leg is drawn by
        # ggrepel (from label to elbow), so even if ggrepel pushes the
        # label to avoid overlap, the leader stays connected.
        (if (!is.null(lab_leader_df))
           ggplot2::geom_segment(
             data = lab_leader_df,
             mapping = ggplot2::aes(x = mx, y = my,
                                    xend = elbow_x, yend = elbow_y,
                                    color = .data[[group_by]]),
             linewidth = module_label_segment_width,
             alpha     = module_label_segment_alpha,
             lineend   = "round",
             show.legend = FALSE
           )
         else NULL) +
        ggrepel::geom_text_repel(data = lab_df,
                                 mapping = ggplot2::aes(x = x,
                                               y = y,
                                               label = .label_text,
                                               color = .data[[group_by]]),
                                 size = module_label_size,
                                 nudge_x = lab_df$nudge_x,
                                 nudge_y = lab_df$nudge_y,
                                 hjust   = lab_df$hjust,
                                 vjust   = lab_df$vjust,
                                 # ggrepel draws the second leg of the
                                 # two_column_follow leader (and the full
                                 # leader for two_column / label_circle).
                                 # min.segment.length = 0 forces it on.
                                 min.segment.length = 0,
                                 segment.size  = module_label_segment_width,
                                 segment.alpha = module_label_segment_alpha,
                                 # segment.square (>= ggrepel 0.9.4) is
                                 # TRUE only for label_circle so ggrepel
                                 # draws an L-shape there; FALSE elsewhere
                                 segment.square      = label_segment_square,
                                 segment.squareShape = label_segment_square_shape,
                                 max.overlaps = Inf,
                                 box.padding = 0.15,
                                 point.padding = label_point_padding,
                                 force = label_force,
                                 show.legend = FALSE
        ) +
        color_scale_lab +
        ggplot2::coord_equal(clip = "off",
                             xlim = plot_xlim,
                             ylim = plot_ylim) +
        theme_ggnetview()
    }

    # label = FALSE module_outline = TRUE
    if (isFALSE(show_module_label) & isTRUE(module_outline)) {

      maskTable <- .build_mask_table()

      maskTable <- maskTable %>% dplyr::mutate(cluster = factor(cluster, levels = levels(ly1_1[["graph_ly_final"]]$Modularity), ordered = TRUE))

      mask_classes <- .ggnv_class_order(maskTable$cluster)
      fill_scale_mask <- if (is.null(module_palette)) scale_fill_ggnetview(mask_classes, labels = module_label_fun) else ggplot2::scale_fill_manual(values = module_palette, labels = module_label_fun)
      color_scale_mask <- if (is.null(module_palette)) scale_color_ggnetview(mask_classes, labels = module_label_fun) else ggplot2::scale_color_manual(values = module_palette, labels = module_label_fun)

      p1_1 <- p1_1 +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggplot2::geom_polygon(data=maskTable %>%
                                dplyr::filter(cluster != "Others"),
                              mapping = ggplot2::aes(x = x, y = y,
                                                     group = interaction(cluster, polygon_id),
                                                     fill = cluster, color = cluster),
                              linewidth = module_outline_width,
                              linetype = module_outline_linetype,
                              alpha = module_outline_alpha,
                              show.legend = FALSE) +
        fill_scale_mask +
        color_scale_mask +
        ggplot2::coord_equal(clip = "off") +
        theme_ggnetview()
    }

    # label = TRUE module_outline = TRUE
    if (isTRUE(show_module_label) & isTRUE(module_outline)) {

      .build_label_location()
      maskTable <- .build_mask_table()

      maskTable <- maskTable %>% dplyr::mutate(cluster = factor(cluster, levels = levels(ly1_1[["graph_ly_final"]]$Modularity), ordered = TRUE))

      lab_classes_outer <- .ggnv_class_order(lab_df$Modularity)
      mask_classes_outer <- .ggnv_class_order(maskTable$cluster)
      color_scale_lab_outer <- if (is.null(module_palette)) scale_color_ggnetview(lab_classes_outer, labels = module_label_fun) else ggplot2::scale_color_manual(values = module_palette, labels = module_label_fun)
      fill_scale_mask_outer <- if (is.null(module_palette)) scale_fill_ggnetview(mask_classes_outer, na_value = NA, labels = module_label_fun) else ggplot2::scale_fill_manual(values = module_palette, labels = module_label_fun)
      color_scale_mask_outer <- if (is.null(module_palette)) scale_color_ggnetview(mask_classes_outer, na_value = NA, labels = module_label_fun) else ggplot2::scale_color_manual(values = module_palette, labels = module_label_fun)

      p1_1 <- p1_1 +
        ggnewscale::new_scale_color() +
        # First leg of the two-segment leader (two_column_follow only):
        # from module centroid (mx, my) to the elbow on the network
        # boundary (elbow_x, elbow_y). The second leg is drawn by
        # ggrepel (from label to elbow), so even if ggrepel pushes the
        # label to avoid overlap, the leader stays connected.
        (if (!is.null(lab_leader_df))
           ggplot2::geom_segment(
             data = lab_leader_df,
             mapping = ggplot2::aes(x = mx, y = my,
                                    xend = elbow_x, yend = elbow_y,
                                    color = modularity2),
             linewidth = module_label_segment_width,
             alpha     = module_label_segment_alpha,
             lineend   = "round",
             show.legend = FALSE
           )
         else NULL) +
        ggrepel::geom_text_repel(data = lab_df,
                                 mapping = ggplot2::aes(x = x,
                                               y = y,
                                               label = .label_text,
                                               color = modularity2),
                                 size = module_label_size,
                                 nudge_x = lab_df$nudge_x,
                                 nudge_y = lab_df$nudge_y,
                                 hjust   = lab_df$hjust,
                                 vjust   = lab_df$vjust,
                                 # ggrepel draws the second leg of the
                                 # two_column_follow leader (and the full
                                 # leader for two_column / label_circle).
                                 # min.segment.length = 0 forces it on.
                                 min.segment.length = 0,
                                 segment.size  = module_label_segment_width,
                                 segment.alpha = module_label_segment_alpha,
                                 # segment.square (>= ggrepel 0.9.4) is
                                 # TRUE only for label_circle so ggrepel
                                 # draws an L-shape there; FALSE elsewhere
                                 segment.square      = label_segment_square,
                                 segment.squareShape = label_segment_square_shape,
                                 max.overlaps = Inf,
                                 box.padding = 0.15,
                                 point.padding = label_point_padding,
                                 force = label_force,
                                 show.legend = FALSE
        ) +
        color_scale_lab_outer +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggplot2::geom_polygon(data= maskTable %>% dplyr::filter(cluster != "Others"),
                              mapping = ggplot2::aes(x = x, y = y,
                                                     group = interaction(cluster, polygon_id),
                                                     fill = cluster, color = cluster),
                              linewidth = module_outline_width,
                              linetype = module_outline_linetype,
                              alpha = module_outline_alpha,
                              show.legend = FALSE) +
        fill_scale_mask_outer +
        color_scale_mask_outer +
        ggplot2::coord_equal(clip = "off",
                             xlim = plot_xlim,
                             ylim = plot_ylim) +
        theme_ggnetview()
    }

    # add point labels after module boundary/text layers
    if (!is.null(point_label_df) && nrow(point_label_df) > 0) {
      p1_1 <- p1_1 +
        ggplot2::geom_text(
          data = point_label_df,
          mapping = ggplot2::aes(x = x, y = y, label = .data[[point_label_col]]),
          size = node_label_size,
          show.legend = FALSE
        )
    }


    if (isTRUE(show_sign_stats)) {
      gglabel = paste0("Node = ", stat_graph$node, "\n",
                       "Edge = ", stat_graph$edge, "\n",
                       "Positive = ", stat_graph$position_edge, "\n",
                       "Negative = ", stat_graph$negative_edge)
    }else{
      gglabel = paste0("Node = ", stat_graph$node, "\n",
                       "Edge = ", stat_graph$edge, "\n")
    }

    p1_1 <- p1_1 +
      ggplot2::ggtitle(label = gglabel) +
      ggplot2::theme(
        legend.box = "horizontal",
        legend.box.just = "left"
      )

  }

  # specific layout dendrogram
  if (layout == "dendrogram") {
    color_default_dendro <- c('#66c2a5','#fc8d62','#a6d854','#e78ac3')
    color_scale_dendro <- if (is.null(node_color_values)) {
      color_default_dendro
    } else {
      node_color_values
    }
    p1_1 <- ggraph::ggraph(graph_obj,layout = layout, circular = TRUE) +
      ggraph::geom_node_point(ggplot2::aes(size=node_size, color=type),alpha=node_alpha) +
      ggraph::geom_edge_diagonal(ggplot2::aes(color = node1.node), alpha=edge_alpha) +
      ggraph::scale_edge_color_manual(values = color_scale_dendro) +
      ggplot2::scale_color_manual(values = color_scale_dendro) +
      ggplot2::scale_size(range = c(3,15)) +
      ggraph::geom_node_text(
        ggplot2::aes(
          x = 1.0175 * x,
          y = 1.0175 * y,
          label = node,
          angle = -((-ggraph::node_angle(x, y) + 90) %% 180) + 90,
          filter = leaf,
          color = type
        ),
        size = 2, hjust = 'outward'
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label=node,
            filter = !leaf,
            color = type),
        fontface="bold",
        size=3,
        family="sans"
      ) +
      ggplot2::coord_fixed(clip = "off") +
      theme_ggnetview()

    return(p1_1)
  }

  # specific layout pie
  if (group_by == "pie") {

    ly <- ggraph::create_layout(graph_obj, layout = layout)

    col_index_start = which(colnames(ly) == "name")
    col_index_end = which(colnames(ly) == ".ggraph.orig_index")
    col_index = colnames(ly)[(col_index_start+1) : (col_index_end -1)]


    fill_default_pie <- c('#66c2a5','#fc8d62','#a6d854','#e78ac3')
    fill_scale_pie <- if (is.null(node_fill_values)) {
      ggplot2::scale_fill_manual(values = fill_default_pie)
    } else {
      ggplot2::scale_fill_manual(values = node_fill_values)
    }
    p1_1 <- ggraph::ggraph(ly, layout = "manual", x = ly[["x"]], y = ly[["y"]]) +
      ggraph::geom_edge_link(color = "#6baed6") +
      scatterpie::geom_scatterpie(
        data = ly,
        cols = col_index,
        colour = "#000000",
        pie_scale = 2
      ) +
      fill_scale_pie +
      ggplot2::coord_fixed() +
      theme_ggnetview()

    return(p1_1)
  }

  if (isTRUE(return_layout) && exists("ly1_1", inherits = FALSE) &&
      !is.null(ly1_1$graph_ly_final) && "Modularity" %in% colnames(ly1_1$graph_ly_final)) {
    module_centroids <- ly1_1$graph_ly_final %>%
      dplyr::filter(as.character(.data$Modularity) != "Others") %>%
      dplyr::group_by(.data$Modularity) %>%
      dplyr::summarise(x = mean(.data$x, na.rm = TRUE), y = mean(.data$y, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(ID = as.character(.data$Modularity)) %>%
      dplyr::select("ID", "x", "y")
    layout_data <- list(
      graph_ly_final = ly1_1$graph_ly_final,
      graph_obj = ly1_1$graph_obj,
      ggplot_data = ly1_1$ggplot_data,
      module_centroids = module_centroids
    )
    return(list(plot = p1_1, layout_data = layout_data))
  }

  return(p1_1)
}
