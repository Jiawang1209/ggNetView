#' Visualize a network with a magnified (zoomed) module subgraph
#'
#' Reproduces the common "local-magnification" network figure: the full
#' network is drawn on the left (optionally with the selected module(s)
#' outlined), an arrow points to the right, and the extracted module subgraph
#' is redrawn as its own panel, with an optional node / edge / component
#' summary in its subtitle.
#'
#' The function is a thin, deterministic wrapper that reuses the package's
#' own machinery end to end: the full network is rendered by
#' \code{\link{ggNetView}} (so network construction, module colouring, layouts
#' and seeds are identical to the rest of the package), the subgraph is
#' extracted by \code{\link{get_subgraph}}, and the panels are composed with
#' \code{patchwork}. By default (\code{sub_layout = "same"}) the magnified
#' panel inherits the exact node coordinates of the full network, so it is a
#' true zoom of that module rather than a re-layout. No non-CRAN dependency
#' (e.g. \code{ggmagnify}) is required.
#'
#' @param graph_obj A \code{tbl_graph} object from
#'   \code{\link{build_graph_from_mat}} / \code{\link{build_graph_from_df}}
#'   (or any other \code{build_graph_from_*}). Its node table must contain a
#'   \code{Modularity} column, exactly as produced by those builders.
#' @param select_module Character or numeric vector. The module name(s) (from
#'   \code{levels(Modularity)}) to extract into the magnified panel. Multiple
#'   modules are allowed; they are combined into a single induced subgraph.
#' @param full_layout Character (default \code{"gephi"}). Layout used for the
#'   full network panel; any layout accepted by \code{\link{ggNetView}}.
#' @param sub_layout Character (default \code{"same"}). Layout of the magnified
#'   subgraph panel. \code{"same"} inherits the full network's node coordinates
#'   (a true zoom, identical relative positions). Any other value is treated as
#'   a \code{\link{ggNetView}} layout name (e.g. \code{"circle"}, \code{"gephi"},
#'   \code{"fr"}) and re-lays-out the subgraph independently. When
#'   \code{select_module} names several modules, the subgraph contains exactly
#'   those modules, so the multipartite layouts become applicable:
#'   \code{"bipartite_gephi_layout"} for 2 modules,
#'   \code{"tripartite_gephi_layout"} for 3,
#'   \code{"quadripartite_gephi_layout"} for 4 and
#'   \code{"pentapartite_gephi_layout"} for 5 (each requires the module count to
#'   match). Unused module levels are dropped automatically so the count is
#'   exact.
#' @param full_args,sub_args Named lists of extra arguments passed through to
#'   \code{\link{ggNetView}} for the full-network and subgraph panels
#'   respectively (e.g. \code{full_args = list(module_label = TRUE, node_size_range = c(2, 6))}).
#'   Deprecated pre-0.2.0 names are still accepted with a lifecycle warning.
#'   Anything you can pass to \code{ggNetView()} can go here, including
#'   \code{module_outline = TRUE} to outline modules (with \code{module_outline_q},
#'   \code{module_outline_alpha}, \code{module_outline_width}, etc.). On the full network this is
#'   drawn by \code{ggNetView} and outlines every module. On the subgraph the
#'   outline is drawn by this function (for \emph{any} \code{sub_layout},
#'   including \code{"same"}) using the full-network palette, so the outer
#'   colour always matches the node / full-network colours; it outlines the
#'   selected module(s). For \code{sub_layout = "same"} the subgraph is drawn
#'   directly, so only \code{edge_color}, \code{edge_alpha} and the
#'   \code{module_outline} styling keys in \code{sub_args} take effect there.
#' @param sub_fill Optional single colour used to recolour every node of the
#'   magnified subgraph (as in the classic teal "extracted module" figure).
#'   When \code{NULL} (default) the subgraph keeps its original module colour,
#'   so it visually matches the module in the full network.
#' @param sub_node_size_range Numeric length-2 vector (default \code{c(4, 10)}). Point
#'   size range for the enlarged subgraph.
#' @param sub_pointsize \lifecycle{deprecated} Use \code{sub_node_size_range}.
#' @param arrow Logical (default \code{TRUE}). Whether to draw the connecting
#'   arrow panel between the full network and the magnified subgraph.
#' @param show_stats Logical (default \code{TRUE}). Whether to show a
#'   nodes / edges / components summary in the subgraph panel's subtitle.
#' @param full_title Character (default \code{"Full Network"}). Title of the
#'   left panel.
#' @param sub_title Character or \code{NULL}. Title of the right panel. When
#'   \code{NULL} (default) a title of the form
#'   \code{"Extracted Subgraph (Module X)"} is generated automatically.
#' @param widths Numeric length-3 vector (default \code{c(1, 0.12, 0.62)})
#'   giving the relative widths of the full-network, arrow and subgraph panels;
#'   the subgraph panel is deliberately smaller than the full network. When
#'   \code{arrow = FALSE} the middle entry is ignored.
#' @param seed Integer (default \code{1115}). Seed forwarded to
#'   \code{\link{ggNetView}} for both panels to keep the figure reproducible.
#'
#' @returns A \code{patchwork} object (a composed \code{ggplot}) that can be
#'   printed, further modified with \code{+}, or saved with
#'   \code{ggplot2::ggsave()}.
#'
#' @seealso \code{\link{get_subgraph}}, \code{\link{ggNetView}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a microbial co-occurrence network exactly as in a normal
#' # ggNetView workflow.
#' data(otu_rare_relative)
#' data(tax_tab)
#' obj <- build_graph_from_mat(
#'   mat              = otu_rare_relative,  # variables x samples abundance matrix
#'   transfrom.method = "none",             # pre-correlation transform
#'   method           = "WGCNA",            # WGCNA::corAndPvalue backend
#'   cor.method       = "pearson",          # Pearson correlation
#'   proc             = "BH",               # Benjamini-Hochberg correction
#'   r.threshold      = 0.7,                # |r| edge cutoff
#'   p.threshold      = 0.05,               # adjusted p-value cutoff
#'   node_annotation  = tax_tab             # taxonomy joined onto nodes
#' )
#'
#' # The module is magnified keeping the SAME layout as the full network (a
#' # true zoom). Full-panel styling is passed through via `full_args`.
#' full_style <- list(
#'   layout_module   = "adjacent",  # neighbouring modules close together
#'   node_size_range = c(1, 5),     # node size range
#'   center        = FALSE,       # do not pull nodes to module centre
#'   shrink        = 0.9,         # compact layout
#'   edge_alpha    = 0.2,         # edge transparency
#'   edge_color    = "#d9d9d9"    # edge colour
#' )
#' ggnetview_subgraph(obj, select_module = "1", full_args = full_style)
#'
#' # Outline modules via ggNetView's own module_outline, on BOTH panels.
#' ggnetview_subgraph(
#'   obj, select_module = "1",
#'   full_args = c(full_style, list(module_outline = TRUE)),
#'   sub_args  = list(module_outline = TRUE)
#' )
#'
#' # Re-lay-out the subgraph as a clean circle, recoloured teal.
#' ggnetview_subgraph(
#'   obj,
#'   select_module = "1",
#'   sub_layout    = "circle",
#'   sub_fill      = "#2C7C82",
#'   full_args     = full_style
#' )
#'
#' # Select THREE modules and lay the subgraph out as a tripartite network.
#' ggnetview_subgraph(
#'   obj,
#'   select_module = c("1", "2", "3"),
#'   sub_layout    = "tripartite_gephi_layout",
#'   full_args     = full_style
#' )
#' }
ggnetview_subgraph <- function(
    graph_obj,
    select_module,
    full_layout      = "gephi",
    sub_layout       = "same",
    full_args        = list(),
    sub_args         = list(),
    sub_fill         = NULL,
    sub_node_size_range = c(4, 10),
    arrow            = TRUE,
    show_stats       = TRUE,
    full_title       = "Full Network",
    sub_title        = NULL,
    widths           = c(1, 0.12, 0.62),
    seed             = 1115,
    sub_pointsize    = deprecated()
) {

  ## ---- validate ----------------------------------------------------------
  if (missing(select_module) || is.null(select_module) ||
      length(select_module) == 0L) {
    stop("`select_module` must name at least one module to magnify.")
  }
  select_module <- as.character(select_module)

  # translate deprecated (< 0.2.0) argument names
  if (lifecycle::is_present(sub_pointsize)) {
    lifecycle::deprecate_warn("0.2.0", "ggnetview_subgraph(sub_pointsize)",
                              "ggnetview_subgraph(sub_node_size_range)")
    sub_node_size_range <- sub_pointsize
  }
  full_args <- .ggnv_rename_args(full_args, fn = "ggnetview_subgraph",
                                 env = environment(), user_env = parent.frame())
  sub_args  <- .ggnv_rename_args(sub_args,  fn = "ggnetview_subgraph",
                                 env = environment(), user_env = parent.frame())

  node_tbl <- graph_obj %>%
    tidygraph::activate("nodes") %>%
    tidygraph::as_tibble()
  if (!"Modularity" %in% colnames(node_tbl)) {
    stop("`graph_obj` has no `Modularity` node column; build it with a ",
         "ggNetView `build_graph_from_*()` function first.")
  }
  available <- as.character(unique(node_tbl[["Modularity"]]))
  missing_mod <- setdiff(select_module, available)
  if (length(missing_mod) > 0L) {
    stop("Module(s) not found in `graph_obj`: ",
         paste(missing_mod, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "), ".")
  }

  ## ---- full network panel (with layout coordinates) ----------------------
  # Everything in `full_args` is forwarded to ggNetView, including
  # `add_outer = TRUE` (which ggNetView draws itself for every module).
  full_call <- utils::modifyList(
    list(
      graph_obj     = graph_obj,
      layout        = full_layout,
      seed          = seed,
      return_layout = TRUE
    ),
    full_args
  )

  full_res  <- do.call(ggNetView, full_call)
  if (!is.list(full_res) || is.null(full_res[["plot"]])) {
    stop("`full_layout = \"", full_layout, "\"` does not expose layout ",
         "coordinates. Use a standard layout such as \"gephi\", \"fr\" or ",
         "\"circle\".")
  }
  full_plot <- full_res[["plot"]]
  full_ly   <- full_res[["layout_data"]][["graph_ly_final"]]

  # Nodes belonging to the selected module(s).
  sel_ly <- full_ly[as.character(full_ly[["Modularity"]]) %in% select_module,
                    , drop = FALSE]

  # Canonical class order + palette of the FULL network. Using ggNetView's own
  # `.ggnv_class_order()` guarantees each module gets EXACTLY the colour it has
  # in the full-network panel, both for any outer boundary and the subgraph.
  # (Marking the full network is done via `full_args = list(add_outer = TRUE)`,
  # which ggNetView draws itself.)
  mod_classes  <- .ggnv_class_order(full_ly[["Modularity"]])
  full_palette <- get_palette(mod_classes)

  # Legend on the left (the magnified subgraph sits on the right).
  full_plot <- full_plot +
    ggplot2::theme(legend.position = "left") +
    ggplot2::ggtitle(full_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ## ---- extract + magnify subgraph ----------------------------------------
  sg      <- get_subgraph(graph_obj, select_module = select_module)
  sub_g   <- sg[["sub_graph_select"]]

  # Drop unused module levels so the subgraph reports EXACTLY the selected
  # modules. Without this the retained (empty) factor levels can break the
  # module count that the multipartite layouts (bipartite / tripartite /
  # quadripartite / pentapartite) require.
  sub_g <- sub_g %>%
    tidygraph::activate("nodes") %>%
    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("Modularity", "modularity2", "modularity3")),
      function(col) if (is.factor(col)) droplevels(col) else col
    ))

  if (identical(sub_layout, "same")) {
    # Inherit the exact node coordinates from the full-network layout, so the
    # magnified panel is a true zoom of that module (same relative positions),
    # reusing the edge endpoints already computed by get_location().
    gd        <- full_res[["layout_data"]][["ggplot_data"]]
    node_df   <- gd[["ggplot_node_df"]]
    edge_df   <- gd[["ggplot_edge_df"]]
    sel_names <- as.character(sel_ly[["name"]])

    nodes_sel <- node_df[as.character(node_df[["Modularity"]]) %in% select_module,
                         , drop = FALSE]
    edges_sel <- edge_df[edge_df[["from_id"]] %in% sel_names &
                         edge_df[["to_id"]]   %in% sel_names, , drop = FALSE]

    # node size by within-subgraph degree (conveys hub structure, enlarged)
    deg <- igraph::degree(tidygraph::as.igraph(sub_g))
    nodes_sel[[".mag_deg"]] <-
      as.numeric(deg[as.character(nodes_sel[["name"]])])

    # edge styling inherits from sub_args when provided
    sub_linecolor <- if (!is.null(sub_args[["edge_color"]]) &&
                         !sub_args[["edge_color"]] %in% colnames(edge_df)) sub_args[["edge_color"]] else "grey70"
    sub_linealpha <- if (!is.null(sub_args[["edge_alpha"]])) sub_args[["edge_alpha"]] else 0.5

    # `mod_classes` (computed above) keeps the palette identical to the
    # full-network panel.
    fill_scale <- if (is.null(sub_fill)) {
      scale_fill_ggnetview(mod_classes)
    } else {
      ggplot2::scale_fill_manual(values = stats::setNames(
        rep(sub_fill[1L], length(select_module)), select_module))
    }

    sub_plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data    = edges_sel,
        mapping = ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
        colour  = sub_linecolor, alpha = sub_linealpha, linewidth = 0.6
      ) +
      ggplot2::geom_point(
        data    = nodes_sel,
        mapping = ggplot2::aes(x = x, y = y, fill = Modularity, size = .mag_deg),
        shape = 21, stroke = 0.3, colour = "grey20"
      ) +
      ggplot2::scale_size(range = sub_node_size_range, guide = "none") +
      fill_scale +
      ggplot2::guides(fill = "none") +
      ggplot2::coord_fixed(clip = "off") +
      theme_ggnetview()

    sub_ly       <- nodes_sel[, c("name", "x", "y")]
    sub_xy       <- nodes_sel[, c("x", "y")]
    sub_clusters <- as.character(nodes_sel[["Modularity"]])

  } else {
    # Re-lay-out the subgraph independently with a named ggNetView layout
    # (e.g. "circle", "gephi", "fr").
    sub_call <- utils::modifyList(
      list(
        graph_obj     = sub_g,
        layout        = sub_layout,
        node_size_range = sub_node_size_range,
        seed          = seed,
        return_layout = TRUE
      ),
      sub_args
    )
    # We draw the subgraph's add_outer ourselves (below) with the full-network
    # palette, so strip it here; ggNetView's own mask colours it from
    # `modularity3` + its own scale and can mismatch the node colours in
    # module-structured layouts (e.g. circular_modules_*).
    sub_call[["module_outline"]] <- NULL
    if (!is.null(sub_fill)) {
      sub_call[["node_fill_values"]] <- stats::setNames(
        rep(sub_fill[1L], length(select_module)), select_module
      )
    } else {
      # Force the FULL network's palette so the subgraph module keeps its
      # original colour (otherwise ggNetView would recolour it from scratch,
      # giving module "14" the first palette colour instead of its own).
      sub_call[["node_fill_values"]] <- full_palette
    }
    sub_res  <- do.call(ggNetView, sub_call)
    if (is.list(sub_res) && !is.null(sub_res[["plot"]])) {
      sub_plot     <- sub_res[["plot"]]
      sub_ly       <- sub_res[["layout_data"]][["graph_ly_final"]]
      sub_xy       <- sub_ly[, c("x", "y")]
      sub_clusters <- as.character(sub_ly[["Modularity"]])
    } else {
      # Layouts such as "dendrogram"/"pie" return a bare ggplot.
      sub_plot     <- sub_res
      sub_ly       <- NULL
      sub_xy       <- NULL
      sub_clusters <- NULL
    }
  }

  if (is.null(sub_title)) {
    sub_title <- paste0("Extracted Subgraph (Module ",
                        paste(select_module, collapse = ", "), ")")
  }

  # Node / edge / component summary goes into the panel subtitle (outside the
  # plotting area) rather than a box inside the network.
  sub_subtitle <- NULL
  if (isTRUE(show_stats)) {
    ig      <- tidygraph::as.igraph(sub_g)
    sub_subtitle <- paste0(
      "Nodes: ", igraph::vcount(ig),
      "  |  Edges: ", igraph::ecount(ig),
      "  |  Components: ", igraph::components(ig)[["no"]]
    )
  }

  # add_outer on the subgraph: drawn HERE for BOTH the "same" and re-layout
  # paths, using the full-network palette, so the outer colour always matches
  # the node/full-network colours (ggNetView's own mask is bypassed above).
  # Opt in via `sub_args = list(add_outer = TRUE)`; styling reads the usual
  # ggNetView outer_* keys from `sub_args`.
  getd <- function(x, d) if (is.null(x)) d else x
  if (isTRUE(sub_args[["module_outline"]]) && !is.null(sub_xy)) {
    sub_plot <- .ggnv_add_outer_layer(
      sub_plot,
      xy_df     = sub_xy,
      clusters  = sub_clusters,
      classes   = mod_classes,
      q         = getd(sub_args[["module_outline_q"]], 0.88),
      expand    = getd(sub_args[["module_outline_expand"]], 1.02),
      bandwidth = getd(sub_args[["module_outline_bandwidth"]], 2),
      linewidth = getd(sub_args[["module_outline_width"]], 1),
      linetype  = getd(sub_args[["module_outline_linetype"]], 1),
      alpha     = getd(sub_args[["module_outline_alpha"]], 0.5),
      fb_color = "grey30", fb_linetype = "dashed",
      fb_linewidth = 0.8, fb_expand_mm = 8
    )
  }

  # Keep only the full-network legend (on the left); drop the subgraph's own.
  sub_plot <- sub_plot +
    ggplot2::labs(title = sub_title, subtitle = sub_subtitle) +
    ggplot2::guides(fill = "none", colour = "none", size = "none") +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(hjust = 0.5)
    )

  ## ---- compose panels -----------------------------------------------------
  if (isTRUE(arrow)) {
    arrow_panel <- ggplot2::ggplot() +
      ggplot2::annotate(
        "segment",
        x = 0, xend = 1, y = 0.5, yend = 0.5,
        arrow = grid::arrow(type = "closed",
                            length = grid::unit(0.16, "inches")),
        linewidth = 1.1
      ) +
      ggplot2::xlim(-0.1, 1.1) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_void()

    combined <- patchwork::wrap_plots(
      full_plot, arrow_panel, sub_plot,
      nrow = 1, widths = widths
    )
  } else {
    combined <- patchwork::wrap_plots(
      full_plot, sub_plot,
      nrow = 1, widths = widths[c(1L, 3L)]
    )
  }

  combined
}


#' Add an HDR module-outline layer to a ggNetView plot
#'
#' Internal helper shared by the full-network and subgraph panels of
#' \code{\link{ggnetview_subgraph}}. Draws the same HDR contour that
#' \code{ggNetView(module_outline = TRUE)} uses, but for an explicit set of node
#' coordinates / clusters, and coloured from a fixed class palette so colours
#' stay consistent with the rest of the figure. Falls back to a dashed ring
#' when the point cloud is too small for a KDE contour.
#'
#' @noRd
.ggnv_add_outer_layer <- function(p, xy_df, clusters, classes,
                                  q, expand, bandwidth,
                                  linewidth, linetype, alpha,
                                  fb_color, fb_linetype, fb_linewidth,
                                  fb_expand_mm) {
  maskTable <- tryCatch(
    generateMask_ggnetview(
      dims            = xy_df,
      clusters        = clusters,
      q               = q,
      expand          = expand,
      bandwidth_scale = bandwidth
    ),
    error = function(e) NULL
  )
  if (!is.null(maskTable) && nrow(maskTable) > 0L) {
    p +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color() +
      ggplot2::geom_polygon(
        data    = maskTable,
        mapping = ggplot2::aes(x = x, y = y,
                               group = interaction(cluster, polygon_id),
                               fill = cluster, colour = cluster),
        linewidth   = linewidth,
        linetype    = linetype,
        alpha       = alpha,
        show.legend = FALSE
      ) +
      scale_fill_ggnetview(classes) +
      scale_color_ggnetview(classes)
  } else {
    # Too few nodes for a KDE contour: fall back to a dashed ring.
    df <- as.data.frame(xy_df)
    df[[".mag_group"]] <- 1L
    p +
      ggforce::geom_mark_circle(
        data    = df,
        mapping = ggplot2::aes(x = x, y = y, group = .mag_group),
        fill = NA, color = fb_color, linetype = fb_linetype,
        linewidth = fb_linewidth, expand = grid::unit(fb_expand_mm, "mm")
      )
  }
}
