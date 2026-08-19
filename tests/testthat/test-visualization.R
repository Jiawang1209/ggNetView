test_that("ggNetView returns a ggplot object with fr layout", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(g, layout = "fr", seed = 1, module_label = FALSE)

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("ggNetView works with circle layout", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(g, layout = "circle", seed = 1, module_label = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("ggNetView works with gephi layout", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(g, layout = "gephi", seed = 1, module_label = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("ggNetView works with kk layout", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(g, layout = "kk", seed = 1, module_label = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("ggNetView respects node_fill parameter", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(
    g, layout = "fr", seed = 1,
    node_fill = "Modularity", module_label = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("ggNetView works with star layout", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  p <- ggNetView(g, layout = "star", seed = 1, module_label = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("ggNetView return_layout produces a list with layout data", {
  data(ppi_example, package = "ggNetView")
  g <- build_graph_from_df(
    df = ppi_example$ppi,
    node_annotation = ppi_example$annotation
  )

  res <- ggNetView(
    g, layout = "fr", seed = 1,
    module_label = FALSE, return_layout = TRUE
  )

  expect_true(is.list(res))
})
