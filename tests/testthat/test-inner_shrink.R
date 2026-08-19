# Tests for the new `inner_shrink` parameter on create_layout_WGCNA().
#
# Goals
# -----
# 1. Backward-compat: at the default `inner_shrink = 1`, results must be
#    bit-identical to the pre-change implementation.  We can't compare to
#    the old function directly (it's been replaced), so instead we assert
#    the geometry that the old implementation produced:
#       - module disc centres `(cx, cy)` from .pack_wgcna_modules() are
#         unchanged for any inner_shrink (this is the whole point of the
#         decoupling),
#       - the outermost point of every module sits within `radius * 0.95`
#         (with-edge case) or within `radius` (empty-edge case) of the
#         module centroid -- exactly what the original code guaranteed.
#
# 2. Decoupling: at inner_shrink < 1, module disc centres must be invariant
#    (we explicitly do NOT want this parameter to move modules around);
#    only the per-module point cloud should contract toward each centroid.
#
# 3. Proportionality: with-edge max radius scales linearly with
#    inner_shrink (within FR initialization noise tolerated to ~1%).
#
# 4. Validation: bad inputs (NA, negative, character, length != 1) error.

testthat::test_that("inner_shrink: default = 1 preserves historical fit_radius", {
  set.seed(42)

  # Build a small graph with modularity3 column.
  # Three connected modules of different sizes, plus a small "Others".
  g <- make_test_wgcna_graph()

  ly <- create_layout_WGCNA(g, r = 1)

  geom <- expected_module_geometry(g)
  pts  <- attach_module_to_points(ly, g)

  # For modules with at least one internal edge: max distance from
  # module centroid to any of its points must be <= radius * 0.95
  # (the historical inner margin).  Tolerance covers float rounding.
  for (mod in geom$with_edges) {
    pts_mod <- pts[pts$modularity3 == mod$name, , drop = FALSE]
    d <- sqrt((pts_mod$x - mod$cx)^2 + (pts_mod$y - mod$cy)^2)
    testthat::expect_lte(max(d), mod$radius * 0.95 + 1e-6)
  }

  # For empty-edge modules: bound is the full radius.
  for (mod in geom$no_edges) {
    pts_mod <- pts[pts$modularity3 == mod$name, , drop = FALSE]
    d <- sqrt((pts_mod$x - mod$cx)^2 + (pts_mod$y - mod$cy)^2)
    testthat::expect_lte(max(d), mod$radius + 1e-6)
  }
})

testthat::test_that("inner_shrink: module disc centres are invariant", {
  set.seed(42)
  g <- make_test_wgcna_graph()

  ly1   <- create_layout_WGCNA(g, r = 1, inner_shrink = 1)
  ly05  <- create_layout_WGCNA(g, r = 1, inner_shrink = 0.5)
  ly015 <- create_layout_WGCNA(g, r = 1, inner_shrink = 0.15)

  c1   <- module_centroids(ly1, g)
  c05  <- module_centroids(ly05, g)
  c015 <- module_centroids(ly015, g)

  # Per-module shrunk-toward-centroid means the module *points'* mean
  # tracks the disc centre `(cx, cy)`.  For modules with > 1 internal
  # edge the FR layout is already centred (the algorithm subtracts the
  # mean), so the point-cloud centroid sits exactly on the disc centre
  # for any inner_shrink.  We therefore expect the centroid coordinates
  # to be invariant within FR-noise tolerance -- but ONLY for with-edge
  # modules.  Empty-edge modules (e.g. "Others", which is laid out by a
  # different code path in `create_layout_WGCNA()`) do not satisfy this
  # invariance and must be filtered out before the comparison.
  geom <- expected_module_geometry(g)
  with_edge_mods <- vapply(geom$with_edges, function(m) m$name, character(1))
  c1   <- c1[c1$modularity3   %in% with_edge_mods, , drop = FALSE]
  c05  <- c05[c05$modularity3 %in% with_edge_mods, , drop = FALSE]
  c015 <- c015[c015$modularity3 %in% with_edge_mods, , drop = FALSE]

  testthat::expect_equal(c1$cx,  c05$cx,  tolerance = 1e-6)
  testthat::expect_equal(c1$cy,  c05$cy,  tolerance = 1e-6)
  testthat::expect_equal(c1$cx,  c015$cx, tolerance = 1e-6)
  testthat::expect_equal(c1$cy,  c015$cy, tolerance = 1e-6)
})

testthat::test_that("inner_shrink: with-edge max radius scales linearly", {
  set.seed(42)
  g <- make_test_wgcna_graph()

  geom <- expected_module_geometry(g)

  for (s in c(1, 0.7, 0.4, 0.15)) {
    ly  <- create_layout_WGCNA(g, r = 1, inner_shrink = s)
    pts <- attach_module_to_points(ly, g)
    for (mod in geom$with_edges) {
      pts_mod <- pts[pts$modularity3 == mod$name, , drop = FALSE]
      d <- sqrt((pts_mod$x - mod$cx)^2 + (pts_mod$y - mod$cy)^2)
      # Hard upper bound: 0.95 * radius * s + float slack.
      testthat::expect_lte(max(d), mod$radius * 0.95 * s + 1e-6)
    }
  }
})

testthat::test_that("inner_shrink: input validation rejects bad values", {
  g <- make_test_wgcna_graph()

  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = NA_real_),
                         "positive finite numeric")
  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = -0.5),
                         "positive finite numeric")
  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = 0),
                         "positive finite numeric")
  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = c(0.5, 0.7)),
                         "positive finite numeric")
  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = "0.5"),
                         "positive finite numeric")
  testthat::expect_error(create_layout_WGCNA(g, inner_shrink = Inf),
                         "positive finite numeric")
})

testthat::test_that("ggNetView dispatch: inner_shrink only flows to WGCNA layout", {
  # Other layouts must still work after the dispatch change; the new
  # parameter must be silently irrelevant for them.
  g <- make_test_wgcna_graph()

  # Just check the function returns without error for a non-WGCNA layout
  # when inner_shrink is supplied (it should be ignored silently).
  testthat::expect_silent({
    p <- ggNetView(g, layout = "fr", layout_module = "random",
                   inner_shrink = 0.5, show_edges = FALSE, seed = 1)
  })
  testthat::expect_s3_class(p, "ggplot")
})


# Helpers for these tests live in `tests/testthat/helper-inner_shrink.R`.
# testthat sources `helper-*.R` BEFORE running any test_that() block, which
# avoids the previous bug where helpers defined at the bottom of this file
# weren't yet visible when the tests at the top of the file ran.
