# ggNetView() 0.2.0 argument scheme: new node_* / edge_* names, the
# "column name or literal value" resolver, and the lifecycle layer that keeps
# every pre-0.2.0 name working with a deprecation warning.

.arg_test_graph <- function() {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )
  # a signed correlation-like edge attribute so edge_color = "corr_direction"
  # and edge_width = "weight" have something to map
  edge_names <- colnames(tidygraph::as_tibble(tidygraph::activate(g, edges)))
  g <- g %>% tidygraph::activate(edges)
  if (!"weight" %in% edge_names) g <- tidygraph::mutate(g, weight = 1)
  g %>%
    tidygraph::mutate(
      corr_direction = rep(c("Positive", "Negative"), length.out = dplyr::n())
    ) %>%
    tidygraph::activate(nodes)
}

.layer_aes_names <- function(p) {
  unlist(lapply(p$layers, function(l) names(l$mapping)))
}

test_that("new node_* / edge_* arguments produce a ggplot", {
  g <- .arg_test_graph()
  p <- ggNetView(g, layout = "fr", seed = 1,
                 node_fill = "Modularity", node_shape = 21,
                 node_size = "Degree", node_size_range = c(2, 8),
                 node_alpha = 0.9, node_stroke = 0.4,
                 edge_color = "corr_direction", edge_width = "weight",
                 edge_alpha = 0.6, edge_linetype = 1,
                 module_label = FALSE)
  expect_s3_class(p, "ggplot")
  aes_names <- .layer_aes_names(p)
  expect_true("linewidth" %in% aes_names)   # edge_width mapped
  # edge_color mapped (ggnewscale renames the edge colour aes to colour_new)
  expect_true(any(grepl("^colou?r", aes_names)))
})

test_that("column name vs literal value is resolved per argument", {
  g <- .arg_test_graph()
  # literal colours / sizes -> constants, no extra mappings
  p <- ggNetView(g, layout = "fr", seed = 1,
                 node_fill = "steelblue", node_color = "black",
                 node_size = 3, edge_color = "grey50", edge_width = 1,
                 module_label = FALSE)
  expect_s3_class(p, "ggplot")
  aes_names <- .layer_aes_names(p)
  expect_false("fill" %in% aes_names)
  expect_false("size" %in% aes_names)
  expect_false("linewidth" %in% aes_names)

  # bad column for a numeric-only aesthetic errors clearly
  expect_error(ggNetView(g, layout = "fr", seed = 1, edge_width = "nope"),
               "edge_width")
  expect_error(ggNetView(g, layout = "fr", seed = 1, node_size = "nope"),
               "node_size")
})

test_that("solid shapes (0-20) route the fill mapping to colour", {
  g <- .arg_test_graph()
  p <- ggNetView(g, layout = "fr", seed = 1, node_shape = 16,
                 node_fill = "Modularity", module_label = FALSE)
  expect_s3_class(p, "ggplot")
  pt <- p$layers[[length(p$layers)]]
  # find the point layer
  pt <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]
  expect_true("colour" %in% names(pt$mapping))
  expect_false("fill" %in% names(pt$mapping))

  # mapping both fill and colour on a solid shape warns and drops fill
  expect_warning(
    ggNetView(g, layout = "fr", seed = 1, node_shape = 16,
              node_fill = "Modularity", node_color = "Degree",
              module_label = FALSE),
    "no fill slot"
  )
})

test_that("deprecated argument names still work and warn once each", {
  g <- .arg_test_graph()
  rlang::local_options(lifecycle_verbosity = "warning")

  p_new <- ggNetView(g, layout = "fr", seed = 1,
                     node_fill = "Modularity", node_size_range = c(2, 8),
                     edge_alpha = 0.5, edge_color = "corr_direction",
                     module_label = FALSE)

  warns <- character()
  p_old <- withCallingHandlers(
    ggNetView(g, layout = "fr", seed = 1,
              fill.by = "Modularity", pointsize = c(2, 8),
              linealpha = 0.5, mapping_line = TRUE,
              label = FALSE),
    lifecycle_warning_deprecated = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 5L)
  expect_true(any(grepl("fill.by", warns, fixed = TRUE)))
  expect_true(any(grepl("node_fill", warns, fixed = TRUE)))
  expect_true(any(grepl("edge_color", warns, fixed = TRUE)))

  # same plot either way
  expect_equal(length(p_old$layers), length(p_new$layers))
  expect_equal(p_old$labels$title, p_new$labels$title)
  expect_equal(.layer_aes_names(p_old), .layer_aes_names(p_new))
})

test_that("explicit new name wins over a deprecated alias", {
  g <- .arg_test_graph()
  rlang::local_options(lifecycle_verbosity = "warning")
  p <- suppressWarnings(
    ggNetView(g, layout = "fr", seed = 1,
              node_size_range = c(3, 9), pointsize = c(1, 2),
              module_label = FALSE)
  )
  sz <- p$scales$get_scales("size")
  expect_equal(sz$palette(c(0, 1)), c(3, 9))
})

test_that("removed arguments (ring_n, nodelabsize) warn and are ignored", {
  g <- .arg_test_graph()
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_warning(
    ggNetView(g, layout = "fr", seed = 1, ring_n = 3, module_label = FALSE),
    "ring_n"
  )
})

test_that(".ggnv_rename_args translates a forwarded ... list", {
  rlang::local_options(lifecycle_verbosity = "quiet")
  out <- ggNetView:::.ggnv_rename_args(
    list(add_outer = TRUE, fill = c(a = "red"), color = c(a = "blue"),
         linecolor = "grey", mapping_line = TRUE, node_add = 3)
  )
  expect_true(out$module_outline)
  expect_equal(out$node_fill_values, c(a = "red"))
  expect_equal(out$node_color_values, c(a = "blue"))
  expect_equal(out$edge_color_values, c(a = "blue"))
  expect_equal(out$edge_color, "corr_direction")   # mapping beats linecolor
  expect_equal(out$node_add, 3)
  expect_false("add_outer" %in% names(out))
})

test_that("ggNetView_multi / ggnetview_subgraph accept new and old names", {
  g <- .arg_test_graph()
  rlang::local_options(lifecycle_verbosity = "warning")
  mods <- as.character(unique(
    tidygraph::as_tibble(tidygraph::activate(g, nodes))$Modularity))
  mods <- setdiff(mods, "Others")
  skip_if(length(mods) == 0L)

  p_new <- ggnetview_subgraph(g, select_module = mods[1],
                              full_args = list(module_label = FALSE),
                              sub_args = list(edge_alpha = 0.4),
                              seed = 1)
  expect_s3_class(p_new, "ggplot")

  # two deprecated keys (full_args$label, sub_args$linealpha) -> two warnings
  warns <- character()
  withCallingHandlers(
    ggnetview_subgraph(g, select_module = mods[1],
                       full_args = list(label = FALSE),
                       sub_args = list(linealpha = 0.4),
                       seed = 1),
    lifecycle_warning_deprecated = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 2L)
  expect_true(any(grepl("module_label", warns, fixed = TRUE)))
  expect_true(any(grepl("edge_alpha", warns, fixed = TRUE)))
})

test_that("ggNetView_multi_link accepts new names and translates old ones", {
  skip_on_cran()
  data(otu_rare_relative, package = "ggNetView")
  data(otu_sample, package = "ggNetView")
  rlang::local_options(lifecycle_verbosity = "warning")
  # This test is about argument-name translation, not about statistics, so it
  # uses `proc = "none"` with a high |r| cutoff: that keeps every per-group
  # sub-network non-empty and the call fast. `otu_sample` splits the 18 samples
  # into three groups of six, and six samples cannot support a corrected
  # correlation network over ~20,000 pairwise tests -- with `proc = "BH"` the
  # sub-networks are legitimately empty and there is nothing to lay out.
  mat <- otu_rare_relative[seq_len(200), ]
  base <- list(mat = mat, group_info = otu_sample, transfrom.method = "none",
               r.threshold = 0.95, p.threshold = 0.05, method = "WGCNA",
               proc = "none", layout = "gephi", layout_module = "adjacent",
               top_modules = 5, seed = 1115)

  p_new <- suppressMessages(suppressWarnings(do.call(ggNetView_multi_link, c(base, list(
    node_size_range = c(1, 4), edge_color = "corr_direction",
    module_outline = "circle", link_alpha_node = 0.2, group_label_size = 3)))))
  expect_s3_class(p_new$p, "ggplot")

  base_old <- base; base_old$layout_module <- NULL
  warns <- character()
  # NB: do NOT wrap the call in suppressWarnings() here. A calling handler
  # established inside the expression runs before the ones registered outside
  # it, so an inner suppressWarnings() muffles every warning before the
  # collector below can ever see it -- `warns` would silently stay empty and
  # the assertions would be testing nothing. Swallow the unrelated ggplot2
  # warnings with a second handler instead.
  p_old <- suppressMessages(withCallingHandlers(
    do.call(ggNetView_multi_link, c(base_old, list(
      layout.module = "adjacent", pointsize = c(1, 4), mapping_line = TRUE,
      add_outer = "circle", link_linealpha_node = 0.2, label_size = 3))),
    lifecycle_warning_deprecated = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  ))
  expect_length(warns, 6L)
  expect_true(any(grepl("link_alpha_node", warns, fixed = TRUE)))
  expect_equal(length(p_old$p$layers), length(p_new$p$layers))
})

test_that("NULL node_size_range / edge_width_range fall back to defaults", {
  # callers that build the argument list programmatically (e.g.
  # ggnetview_modularity_heatmaps) may pass NULL for an unset range
  g <- .arg_test_graph()
  p <- ggNetView(g, layout = "fr", seed = 1, node_size_range = NULL,
                 edge_width_range = NULL, module_label = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$get_scales("size")$palette(c(0, 1)), c(1, 10))
})
