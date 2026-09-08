# gglink_heatmap_triple() gained cor.method / cor.use / env_p_adjust /
# link_p_adjust / sig_breaks in 0.2.1. The correlation step used to be a fixed
# psych::corr.test() call, and that function reports raw p-values in `$p` for
# the two-matrix form whatever `adjust` says, so a correction that never
# reached the plot was easy to reintroduce. These tests pin each argument to an
# observable effect on cor_test2(), the engine behind the figure.

env_exp <- function() {
  data(Envdf, package = "ggNetView")
  list(E = Envdf[, 1:6], X = Envdf[, 7:14])
}

test_that("the default correlation arguments keep the documented labels", {
  d <- env_exp()
  base <- ggNetView:::cor_test2(d$E, d$X)
  expect_setequal(
    unique(as.character(base[[3]]$p_value)),
    c("P > 0.05", "0.01 < P <= 0.05", "0.001 <= P <= 0.01", "P < 0.001"))
  expect_setequal(unique(as.character(base[[2]]$p_value)), c("", "**", "***"))
})

test_that("cor.method and cor.use reach the correlation", {
  d <- env_exp()
  base <- ggNetView:::cor_test2(d$E, d$X)
  sp <- ggNetView:::cor_test2(d$E, d$X, cor.method = "spearman")
  expect_false(isTRUE(all.equal(base[[3]]$Value, sp[[3]]$Value)))
  expect_equal(dim(sp[[3]]), dim(base[[3]]))
  expect_s3_class(ggNetView:::cor_test2(d$E, d$X, cor.use = "everything")[[3]],
                  "data.frame")
})

test_that("env_p_adjust and link_p_adjust actually correct the p-values", {
  d <- env_exp()
  base <- ggNetView:::cor_test2(d$E, d$X)
  link_bh <- ggNetView:::cor_test2(d$E, d$X, link_p_adjust = "BH")
  env_bh <- ggNetView:::cor_test2(d$E, d$X, env_p_adjust = "BH")

  # a correction can only make a p-value larger, never smaller
  expect_true(all(link_bh[[3]]$Pvalue >= base[[3]]$Pvalue - 1e-12))
  expect_true(all(env_bh[[2]]$Pvalue >= base[[2]]$Pvalue - 1e-12))
  # and it must change something, or the argument is not reaching $p.adj
  expect_false(isTRUE(all.equal(base[[3]]$Pvalue, link_bh[[3]]$Pvalue)))
  expect_false(isTRUE(all.equal(base[[2]]$Pvalue, env_bh[[2]]$Pvalue)))
  # each correction stays on its own side of the figure
  expect_equal(link_bh[[2]]$Pvalue, base[[2]]$Pvalue)
  expect_equal(env_bh[[3]]$Pvalue, base[[3]]$Pvalue)
})

test_that("sig_breaks moves the significance cut points and is validated", {
  d <- env_exp()
  loose <- ggNetView:::cor_test2(d$E, d$X, sig_breaks = c(0.5, 0.4, 0.3))
  expect_setequal(
    unique(as.character(loose[[3]]$p_value)),
    c("P > 0.5", "0.4 < P <= 0.5", "0.3 <= P <= 0.4", "P < 0.3"))
  expect_error(ggNetView:::cor_test2(d$E, d$X, sig_breaks = c(0.01, 0.05, 0.001)),
               "strictly decreasing")
  expect_error(ggNetView:::cor_test2(d$E, d$X, sig_breaks = c(0.05, 0.01)),
               "length 3")
})

test_that("gglink_heatmap_triple forwards the correlation arguments", {
  d <- env_exp()
  envS <- data.frame(Sample = rownames(d$E), d$E, check.names = FALSE)
  expS <- data.frame(Sample = rownames(d$X), d$X, check.names = FALSE)
  ev <- colnames(d$E); xv <- colnames(d$X)
  set.seed(1)
  edge <- unique(data.frame(from = sample(xv, 20, TRUE), to = sample(ev, 20, TRUE)))
  node <- data.frame(node = sample(c(ev, xv)))

  build <- function(...) suppressMessages(ggplot2::ggplot_build(
    gglink_heatmap_triple(Environment = envS, Experiment = expS,
                          edge = edge, node = node, ...)))
  # the heatmap layer carries the environment correlations
  base <- build()
  sp <- build(cor.method = "spearman")
  expect_false(isTRUE(all.equal(base$data, sp$data)))
})
