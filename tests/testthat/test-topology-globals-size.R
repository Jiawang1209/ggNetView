# get_network_topology() runs its null-model loop through future_lapply() under
# a sequential plan. The iteration closures capture this function's environment
# (the abundance matrix, the filtered adjacency), so future's globals accounting
# can exceed the default 500 MiB ceiling and abort a call that transfers
# nothing. The function raises the ceiling for the duration and must put the
# caller's option back (0.2.1).

small_graph <- function(seed = 1) {
  set.seed(seed)
  mat <- matrix(sample(1:10, 24, replace = TRUE), nrow = 4)
  rownames(mat) <- paste0("A", seq_len(nrow(mat)))
  colnames(mat) <- paste0("S", seq_len(ncol(mat)))
  list(mat = mat, graph = build_graph_from_mat(
    mat = mat, method = "cor", cor.method = "pearson",
    r.threshold = 0, p.threshold = 1, seed = seed))
}

test_that("get_network_topology ignores a small future.globals.maxSize", {
  d <- small_graph()
  # far below the size of the captured environment: without the override the
  # call aborts with "exceeds the maximum allowed size"
  rlang::local_options(future.globals.maxSize = 1024)
  out <- get_network_topology(
    graph_obj = d$graph, mat = d$mat, method = "cor", cor.method = "pearson",
    r.threshold = 0, p.threshold = 1, bootstrap = 2)
  expect_true(is.list(out))
  expect_true("topology" %in% names(out))
  # and the caller's option survives the call unchanged
  expect_equal(getOption("future.globals.maxSize"), 1024)
})

test_that("get_network_topology leaves future.globals.maxSize unset if it was", {
  d <- small_graph()
  rlang::local_options(future.globals.maxSize = NULL)
  expect_null(getOption("future.globals.maxSize"))
  invisible(get_network_topology(
    graph_obj = d$graph, mat = d$mat, method = "cor", cor.method = "pearson",
    r.threshold = 0, p.threshold = 1, bootstrap = 2))
  expect_null(getOption("future.globals.maxSize"))
})
