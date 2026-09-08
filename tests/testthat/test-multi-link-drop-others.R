# `ggNetView_multi_link(drop_others = TRUE)` is a display-only switch (0.2.1):
# the "Others" bucket is dropped from the plotted graphs, but the universe that
# compare_modules_by_overlap() tests must stay the complete network. Removing
# the nodes first re-ran the hypergeometric test on a much smaller universe and
# flipped module pairs in and out of significance.

test_that("drop_others does not change the cross-group module links", {
  skip_on_cran()
  data(otu_rare_relative, package = "ggNetView")
  data(otu_sample, package = "ggNetView")

  # `proc = "none"` with a high |r| cutoff keeps every per-group sub-network
  # non-empty and the call fast; `top_modules = 5` leaves a real "Others"
  # bucket, which is what this test is about.
  base <- list(mat = otu_rare_relative[seq_len(200), ], group_info = otu_sample,
               transfrom.method = "none", r.threshold = 0.95,
               p.threshold = 0.05, method = "WGCNA", proc = "none",
               layout = "gephi", layout_module = "adjacent",
               top_modules = 5, seed = 1115)
  run <- function(drop) suppressMessages(suppressWarnings(
    do.call(ggNetView_multi_link, c(base, list(drop_others = drop)))))
  keep <- run(FALSE)
  drop <- run(TRUE)

  n_nodes <- function(res) sum(vapply(res$graph, function(g)
    nrow(tidygraph::as_tibble(g, "nodes")), numeric(1)))

  # the switch must actually remove something, or the test is vacuous
  expect_true(any(vapply(keep$graph, function(g)
    any(as.character(tidygraph::as_tibble(g, "nodes")$Modularity) == "Others"),
    logical(1))))
  expect_lt(n_nodes(drop), n_nodes(keep))

  # ... yet the module-overlap table and the module links are untouched
  expect_gt(nrow(keep$info), 0L)
  expect_equal(drop$info, keep$info)
  key <- c("link_level", "group_a", "group_b", "source", "target")
  mod_links <- function(res)
    res$link_info[res$link_info$link_level == "module", key]
  expect_equal(mod_links(drop), mod_links(keep))

  # "Others" is a display bucket, not a community: it is never an endpoint
  expect_false(any(as.character(unlist(keep$info[, c("modA", "modB")])) == "Others"))
})
