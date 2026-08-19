#' @importFrom lifecycle deprecated
NULL

# ---------------------------------------------------------------------------
# Deprecated-argument registry for ggNetView() (ggNetView 0.2.0)
#
# 0.2.0 renamed the plotting arguments of ggNetView() to a tidyverse-style
# `element_property` scheme (node_*, edge_*, module_label_*, module_outline_*,
# network_outline_*).  The old names are still accepted: every old argument is
# kept in the ggNetView() signature with default `deprecated()`, and the
# functions below translate old -> new with a lifecycle warning.
#
# The same registry is reused by functions that forward a `...` list to
# ggNetView() (ggnetview_modularity_heatmaps, ggnetview_subgraph), so users
# get one consistent deprecation message everywhere.
# ---------------------------------------------------------------------------

# old name -> new name(s).  `NA` = the argument was removed and is ignored.
.ggnv_deprecated_arg_map <- list(
  # ---- node ----
  `fill.by`          = "node_fill",
  `color.by`         = "node_color",
  shape              = "node_shape",
  pointsize          = "node_size_range",
  pointalpha         = "node_alpha",
  pointstroke        = "node_stroke",
  fill               = "node_fill_values",
  color              = c("node_color_values", "edge_color_values"),
  jitter             = "node_jitter",
  jitter_sd          = "node_jitter_sd",
  pointlabel         = "node_label",
  pointlabelsize     = "node_label_size",
  nodelabsize        = NA_character_,
  ring_n             = NA_character_,
  # ---- edge ----
  plot_line          = "show_edges",
  linecolor          = "edge_color",
  mapping_line       = "edge_color",
  linealpha          = "edge_alpha",
  curve              = "edge_curve",
  curvature          = "edge_curvature",
  # ---- module label ----
  label              = "module_label",
  labelsize          = "module_label_size",
  labelsegmentsize   = "module_label_segment_width",
  labelsegmentalpha  = "module_label_segment_alpha",
  label_layout       = "module_label_layout",
  label_wrap_width   = "module_label_wrap",
  label_outer_pad    = "module_label_pad",
  # ---- module outline ----
  add_outer          = "module_outline",
  q_outer            = "module_outline_q",
  expand_outer       = "module_outline_expand",
  bandwidth_scale    = "module_outline_bandwidth",
  outerwidth         = "module_outline_width",
  outerlinetype      = "module_outline_linetype",
  outeralpha         = "module_outline_alpha",
  # ---- network outline ----
  add_group_outer            = "network_outline",
  add_group_outer_expand     = "network_outline_expand",
  add_group_outer_color      = "network_outline_color",
  add_group_outer_fill       = "network_outline_fill",
  add_group_outer_fill_alpha = "network_outline_fill_alpha",
  add_group_outer_linetype   = "network_outline_linetype",
  add_group_outer_linewidth  = "network_outline_width",
  # ---- misc ----
  `layout.module`    = "layout_module",
  `group.by`         = "group_by",
  remove             = "hide_others",
  dropOthers         = "drop_others"
)

.ggnv_deprecated_since <- "0.2.0"

# ggNetView_multi_link() shares most of the registry, minus the ggNetView()-only
# arguments, plus its own inner_* / link_* / label_* names.
.ggnv_multi_link_deprecated_arg_map <- c(
  .ggnv_deprecated_arg_map[setdiff(
    names(.ggnv_deprecated_arg_map),
    c("color", "color.by", "shape", "pointalpha", "pointstroke", "pointlabel",
      "pointlabelsize", "nodelabsize", "plot_line", "curve", "curvature",
      "label", "labelsize", "labelsegmentsize", "labelsegmentalpha",
      "label_layout", "label_wrap_width", "label_outer_pad", "remove"))],
  list(
    color                      = "node_color_values",
    inner_curve                = "edge_curve",
    inner_curvature            = "edge_curvature",
    inner_curve_adaptive       = "edge_curve_adaptive",
    inner_curve_adaptive_range = "edge_curve_adaptive_range",
    inner_curve_adaptive_bins  = "edge_curve_adaptive_bins",
    link_linewidth_node        = "link_width_node",
    link_linewidth_module      = "link_width_module",
    link_linealpha_node        = "link_alpha_node",
    link_linealpha_module      = "link_alpha_module",
    label_offset               = "group_label_offset",
    label_size                 = "group_label_size"
  )
)

#' Translate deprecated ggNetView() argument names to their 0.2.0 names
#'
#' @param args Named list of arguments (any mix of old and new names).
#' @param fn   Function name used in the lifecycle message.
#' @param map Registry (old -> new) to use; defaults to the ggNetView() one.
#' @return Named list with only new names.
#' @noRd
.ggnv_rename_args <- function(args, fn = "ggNetView",
                              env = rlang::caller_env(),
                              user_env = rlang::caller_env(2),
                              map = .ggnv_deprecated_arg_map) {
  if (length(args) == 0L) return(list())
  nms <- names(args)
  if (is.null(nms) || any(nms == "")) {
    stop("All arguments forwarded to `", fn, "()` must be named.", call. = FALSE)
  }
  explicit_new <- setdiff(nms, names(map))
  out <- args[explicit_new]

  # `mapping_line` must be resolved after `linecolor` (both map to edge_color):
  # a mapping request always beats a plain colour.
  old_present <- intersect(names(map), nms)
  old_present <- c(setdiff(old_present, "mapping_line"),
                   intersect(old_present, "mapping_line"))

  for (old in old_present) {
    val <- args[[old]]
    new <- map[[old]]
    what <- sprintf("%s(%s)", fn, old)

    if (is.na(new[1L])) {
      lifecycle::deprecate_warn(
        when = .ggnv_deprecated_since, what = what,
        details = "The argument had no effect and is now ignored.",
        env = env, user_env = user_env
      )
      next
    }

    if (identical(old, "mapping_line")) {
      if (isFALSE(val) || is.null(val)) {
        # `mapping_line = FALSE` was the old default -> nothing to translate.
        lifecycle::deprecate_warn(
          when = .ggnv_deprecated_since, what = what,
          with = sprintf("%s(edge_color)", fn),
          details = c(i = "`mapping_line = TRUE` is now `edge_color = \"corr_direction\"`; a column name is passed to `edge_color` directly."),
          env = env, user_env = user_env
        )
        next
      }
      if (isTRUE(val)) val <- "corr_direction"
      lifecycle::deprecate_warn(
        when = .ggnv_deprecated_since, what = what,
        with = sprintf("%s(edge_color)", fn),
        details = c(i = "`mapping_line = TRUE` is now `edge_color = \"corr_direction\"`; a column name is passed to `edge_color` directly."),
        env = env, user_env = user_env
      )
      # mapping beats a fixed colour coming from `linecolor`
      if (!"edge_color" %in% explicit_new) out[["edge_color"]] <- val
      next
    }

    lifecycle::deprecate_warn(
      when = .ggnv_deprecated_since, what = what,
      with = sprintf("%s(%s)", fn, new[1L]),
      env = env, user_env = user_env
    )
    for (n in new) {
      if (n %in% explicit_new) next        # explicit new name wins
      if (is.null(out[[n]])) out[n] <- list(val)
    }
  }
  out
}

#' Collect the deprecated arguments that were actually supplied to a call
#'
#' Used inside ggNetView(): every deprecated argument has default
#' `deprecated()`; this returns a named list of the ones the caller supplied.
#' @noRd
.ggnv_collect_deprecated <- function(env, map = .ggnv_deprecated_arg_map) {
  old_names <- names(map)
  present <- list()
  for (nm in old_names) {
    if (!exists(nm, envir = env, inherits = FALSE)) next
    # `deprecated()` is rlang's missing-arg sentinel, so the only safe probe is
    # missing() evaluated in the calling frame.
    if (eval(call("missing", as.name(nm)), envir = env)) next
    present[nm] <- list(get(nm, envir = env, inherits = FALSE))
  }
  present
}
