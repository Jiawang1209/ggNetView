# Fixes on the modify_argument branch (0.2.0):
# 1. get_sample_subgraph_topology_parallel() must be callable inside
#    tryCatch()/withCallingHandlers() and must not clobber the caller's
#    future plan.
# 2. gglink_heatmap_triple(): hub anchor coordinates must not become NA when
#    ncol(Experiment) > ceiling(ncol(Environment)/2); minimal inputs (no
#    edge$weight, no node$annotation, unordered nodes) must work.

test_that("get_sample_subgraph_topology_parallel works inside tryCatch and preserves plan", {
  skip_on_cran()
  data(otu_rare_relative, package = "ggNetView")
  data(tax_tab, package = "ggNetView")
  g <- build_graph_from_mat(
    mat = otu_rare_relative[1:120, ], transfrom.method = "none",
    r.threshold = 0.7, p.threshold = 0.05, method = "WGCNA",
    proc = "bonferroni", node_annotation = tax_tab, top_modules = 4, seed = 1115)

  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  plan_before <- future::plan()

  res <- tryCatch(
    suppressWarnings(get_sample_subgraph_topology_parallel(
      graph_obj = g, mat = otu_rare_relative[1:120, 1:3],
      bootstrap = 1, seed = 1115, parallel = TRUE, n_workers = 2)),
    error = function(e) e)
  expect_false(inherits(res, "error"))
  expect_true(is.list(res))
  expect_identical(class(future::plan()), class(plan_before))

  res2 <- tryCatch(
    suppressWarnings(get_sample_subgraph_topology_parallel(
      graph_obj = g, mat = otu_rare_relative[1:120, 1:3],
      bootstrap = 1, seed = 1115, parallel = FALSE)),
    error = function(e) e)
  expect_false(inherits(res2, "error"))
  expect_identical(class(future::plan()), class(plan_before))
  expect_equal(res$stat, res2$stat)
})

test_that("cor_test2 hub anchors have no NA when Experiment is wide", {
  data(Envdf, package = "ggNetView")
  s <- ggNetView:::cor_test2(Envdf[, 1:6], Envdf[, 7:14])  # 8 > ceiling(6/2)
  expect_false(anyNA(s[[3]]$start2))
  expect_false(anyNA(s[[3]]$end2))
  # old working regime unchanged: anchors are still 1,3,5,... minus 4
  s2 <- ggNetView:::cor_test2(Envdf[, 1:6], Envdf[, 7:9])
  expect_equal(sort(unique(s2[[3]]$start2)), c(1, 3, 5) - 4)
})

test_that("gglink_heatmap_triple works with minimal inputs and validates hubs", {
  data(Envdf, package = "ggNetView")
  envS <- data.frame(Sample = rownames(Envdf), Envdf[, 1:6], check.names = FALSE)
  expS <- data.frame(Sample = rownames(Envdf), Envdf[, 7:14], check.names = FALSE)
  ev <- colnames(Envdf)[1:6]; xv <- colnames(Envdf)[7:14]
  set.seed(1)
  edge <- unique(data.frame(from = sample(xv, 20, TRUE), to = sample(ev, 20, TRUE)))
  node <- data.frame(node = sample(c(ev, xv)))   # shuffled, no weight/annotation

  p <- gglink_heatmap_triple(Environment = envS, Experiment = expS,
                             edge = edge, node = node)
  expect_s3_class(p, "ggplot")
  # renders without dropped-row warnings (NA hub coordinates would warn)
  ws <- character()
  tf <- tempfile(fileext = ".png"); grDevices::png(tf, 600, 500)
  withCallingHandlers(suppressMessages(print(p)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  grDevices::dev.off()
  expect_length(grep("Removed", ws), 0)

  expect_error(
    gglink_heatmap_triple(envS, expS, edge, data.frame(node = xv)),
    "outer circle")
  expect_error(
    gglink_heatmap_triple(envS, expS, edge, data.frame(node = c(ev, xv[1:3]))),
    "must equal the number of Experiment")
})

# 3. gglink_heatmap_triple() must build its plot with a single coordinate
#    system; a second coord_*() silently replaces the first and makes ggplot2
#    emit "Coordinate system already present".
test_that("gglink_heatmap_triple adds only one coordinate system", {
  data(Envdf, package = "ggNetView")
  envS <- data.frame(Sample = rownames(Envdf), Envdf[, 1:6], check.names = FALSE)
  expS <- data.frame(Sample = rownames(Envdf), Envdf[, 7:14], check.names = FALSE)
  ev <- colnames(Envdf)[1:6]; xv <- colnames(Envdf)[7:14]
  set.seed(1)
  edge <- unique(data.frame(from = sample(xv, 20, TRUE), to = sample(ev, 20, TRUE)))
  node <- data.frame(node = sample(c(ev, xv)))

  ms <- character()
  p <- withCallingHandlers(
    gglink_heatmap_triple(Environment = envS, Experiment = expS,
                          edge = edge, node = node),
    message = function(m) {
      ms <<- c(ms, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_length(grep("[Cc]oordinate system", ms), 0)
  # the surviving coord is the aspect-ratio-preserving one
  expect_equal(p$coordinates$ratio, 1)
})
