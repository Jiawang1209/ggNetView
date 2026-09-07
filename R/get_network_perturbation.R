#' Structural perturbation ("virtual attack") of a network
#'
#' Repeatedly removes nodes from a network -- at random, by targeted
#' attack on a centrality, or by knocking out a named module / node set --
#' recomputes a panel of connectivity-sensitive topology metrics after
#' each removal, and returns the resulting *perturbation curve* together
#' with a single-number robustness index (Schneider R). This is the
#' structural ("type 1") virtual-perturbation analysis: it asks how the
#' network falls apart as nodes are progressively lost.
#'
#' The random strategy is the multi-metric generalisation of the
#' node-removal robustness already computed inside
#' [get_network_topology()]; the targeted / module / manual strategies
#' answer "which nodes (or whole modules) hold the network together?".
#'
#' @section Strategies:
#' \describe{
#'   \item{`"random"`}{Remove a random subset at each fraction, repeated
#'     `bootstrap` times; mean / sd / se are reported. Seeded for
#'     reproducibility.}
#'   \item{`"targeted"`}{Rank nodes by `centrality` and remove them in
#'     order (most-central first by default). The classic intentional
#'     attack -- usually far more damaging than random failure.}
#'   \item{`"module"`}{Knock out every node whose `module_col` label is
#'     in `target`; reported as a before/after comparison.}
#'   \item{`"manual"`}{Knock out the exact node `name`s given in
#'     `target`; reported as a before/after comparison.}
#' }
#'
#' @section Metrics tracked after each removal:
#' \describe{
#'   \item{`LCC_fraction`}{Size of the largest connected component as a
#'     fraction of the *original* node count (the main attack curve).}
#'   \item{`N_components`}{Number of connected components.}
#'   \item{`Natural_connectivity`}{`log(mean(exp(eigenvalues(A))))`; a
#'     spectral robustness measure that varies smoothly and does not
#'     jump discretely the way component counts do.}
#'   \item{`Efficiency`}{Mean of `1 / shortest-path-distance`.}
#'   \item{`Mean_degree`, `Density`, `Transitivity_global`,
#'     `Modularity`}{Standard summaries recomputed on the survivor
#'     subgraph.}
#' }
#'
#' @param graph_obj A `tbl_graph` from any `build_graph_from_*()`
#'   constructor.
#' @param strategy Character. One of `"random"`, `"targeted"`,
#'   `"module"`, `"manual"`.
#' @param centrality Character. Used when `strategy = "targeted"`.
#'   One of `"degree"`, `"strength"`, `"betweenness"`, `"closeness"`,
#'   `"eigenvector"`, `"ivi"`. `"strength"` uses the `weight` edge
#'   attribute; `"ivi"` requires the `influential` package.
#' @param target Character vector. The module label(s) (when
#'   `strategy = "module"`) or node `name`s (when `strategy = "manual"`)
#'   to remove.
#' @param module_col Character (default `"Modularity"`). Node column
#'   holding module labels, used when `strategy = "module"`.
#' @param fractions Numeric vector in `(0, 1]`. Removal fractions for
#'   `"random"` / `"targeted"`. Default `seq(0.05, 1, by = 0.05)`.
#' @param decreasing Logical (default `TRUE`). For `"targeted"`, remove
#'   the most-central nodes first (`TRUE`) or least-central first
#'   (`FALSE`).
#' @param bootstrap Integer (default `100`). Number of random repetitions
#'   for `strategy = "random"`.
#' @param seed Integer (default `123`). Seed for the random strategy, so
#'   results are reproducible.
#' @param plot Logical (default `TRUE`). Attach a ready-made attack-curve
#'   ggplot of `LCC_fraction` (only for `"random"` / `"targeted"`).
#'
#' @returns A list with:
#'   \itemize{
#'     \item `curve`: long data frame (`strategy`, `fraction`, `metric`,
#'       `value`, and for random `value_sd` / `value_se`).
#'     \item `robustness_index`: data frame with the Schneider R-index
#'       (area under the `LCC_fraction` curve) per strategy.
#'     \item `plot`: ggplot of the LCC attack curve, or `NULL`.
#'   }
#'
#' @references
#' Albert R, Jeong H, Barabasi AL (2000). "Error and attack tolerance of
#' complex networks." \emph{Nature} 406:378-382.
#' Schneider CM et al. (2011). "Mitigation of malicious attacks on
#' networks." \emph{PNAS} 108(10):3838-3841.
#'
#' @seealso [get_node_influence()] for abundance-influence propagation;
#'   [press_perturbation()] for the press-perturbation approximation;
#'   [get_network_topology()] for static network-level metrics.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(ppi_example)
#' obj <- build_graph_from_df(
#'   df              = ppi_example$ppi,
#'   node_annotation = ppi_example$annotation
#' )
#' res <- get_network_perturbation(obj, strategy = "targeted",
#'                                 centrality = "degree")
#' head(res$curve)
#' res$robustness_index
#' }
get_network_perturbation <- function(
  graph_obj,
  strategy   = c("random", "targeted", "module", "manual"),
  centrality = c("degree", "strength", "betweenness", "closeness",
                 "eigenvector", "ivi"),
  target     = NULL,
  module_col = "Modularity",
  fractions  = seq(0.05, 1, by = 0.05),
  decreasing = TRUE,
  bootstrap  = 100,
  seed       = 123,
  plot       = TRUE
) {

  if (!inherits(graph_obj, "tbl_graph")) {
    stop("`graph_obj` must be a tbl_graph (e.g. produced by build_graph_from_mat()).",
         call. = FALSE)
  }
  strategy   <- match.arg(strategy)
  centrality <- match.arg(centrality)

  ig <- tidygraph::as.igraph(graph_obj)
  n0 <- igraph::vcount(ig)
  if (n0 == 0L) stop("`graph_obj` has zero vertices.", call. = FALSE)
  vnames <- igraph::V(ig)$name
  if (is.null(vnames)) vnames <- as.character(seq_len(n0))

  # ---- per-step metric panel on a survivor subgraph -----------------------
  .metrics <- function(sub, n_original) {
    nv <- igraph::vcount(sub)
    if (nv == 0L) {
      return(c(LCC_fraction = 0, N_components = 0, Natural_connectivity = 0,
               Efficiency = 0, Mean_degree = 0, Density = 0,
               Transitivity_global = 0, Modularity = 0))
    }
    comps <- igraph::components(sub)
    lcc   <- max(comps$csize) / n_original
    ncomp <- comps$no

    # natural connectivity: log(mean(exp(eigenvalues of adjacency)))
    # Use the log-sum-exp trick to avoid Inf when the largest eigenvalue is
    # large (exp() overflow) in dense subgraphs.
    A  <- as.matrix(igraph::as_adjacency_matrix(sub, sparse = FALSE))
    ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    ev_max <- max(ev)
    nat <- ev_max + log(mean(exp(ev - ev_max)))

    dd <- 1 / igraph::distances(sub)
    diag(dd) <- NA
    eff <- mean(dd, na.rm = TRUE)

    deg <- mean(igraph::degree(sub))
    den <- igraph::edge_density(sub)
    tg  <- igraph::transitivity(sub, type = "global")
    mod <- tryCatch(
      igraph::modularity(sub, igraph::membership(igraph::cluster_fast_greedy(sub))),
      error = function(e) NA_real_
    )
    c(LCC_fraction = lcc, N_components = ncomp, Natural_connectivity = nat,
      Efficiency = eff, Mean_degree = deg, Density = den,
      Transitivity_global = ifelse(is.nan(tg), 0, tg),
      Modularity = mod)
  }

  metric_names <- c("LCC_fraction", "N_components", "Natural_connectivity",
                    "Efficiency", "Mean_degree", "Density",
                    "Transitivity_global", "Modularity")

  # ---- node ranking for targeted attack -----------------------------------
  .rank_nodes <- function() {
    score <- switch(
      centrality,
      degree      = igraph::degree(ig),
      strength    = igraph::strength(ig, weights = igraph::E(ig)$weight),
      betweenness = igraph::betweenness(ig),
      closeness   = igraph::closeness(ig),
      eigenvector = igraph::eigen_centrality(ig)$vector,
      ivi = {
        if (!requireNamespace("influential", quietly = TRUE)) {
          stop("centrality = 'ivi' requires the `influential` package.",
               call. = FALSE)
        }
        gi <- get_node_ivi(graph_obj)
        iv <- gi %>% tidygraph::activate(nodes) %>% tidygraph::as_tibble()
        stats::setNames(iv$IVI, iv$name)[vnames]
      }
    )
    names(score) <- vnames
    order(score, decreasing = decreasing)
  }

  # ---- build the curve ----------------------------------------------------
  if (strategy %in% c("module", "manual")) {
    if (is.null(target) || length(target) == 0L) {
      stop("`target` must be supplied for strategy = '", strategy, "'.",
           call. = FALSE)
    }
    if (strategy == "module") {
      nd <- graph_obj %>% tidygraph::activate(nodes) %>% tidygraph::as_tibble()
      if (!module_col %in% names(nd)) {
        stop(sprintf("Column `%s` not found on the node table.", module_col),
             call. = FALSE)
      }
      rm_names <- nd$name[as.character(nd[[module_col]]) %in% as.character(target)]
    } else {
      rm_names <- intersect(target, vnames)
    }
    if (length(rm_names) == 0L) {
      stop("No nodes matched `target`; nothing to remove.", call. = FALSE)
    }
    before <- .metrics(ig, n0)
    sub    <- igraph::delete_vertices(ig, which(vnames %in% rm_names))
    after  <- .metrics(sub, n0)
    curve <- data.frame(
      strategy = strategy,
      fraction = rep(c(0, length(rm_names) / n0), each = length(metric_names)),
      metric   = rep(metric_names, 2L),
      value    = c(before[metric_names], after[metric_names]),
      value_sd = NA_real_, value_se = NA_real_,
      row.names = NULL, stringsAsFactors = FALSE
    )

  } else if (strategy == "targeted") {
    ord <- .rank_nodes()
    rows <- lapply(c(0, fractions), function(fr) {
      k <- round(n0 * fr)
      keep <- if (k <= 0L) seq_len(n0)
              else if (k >= n0) integer(0)
              else ord[(k + 1L):n0]
      sub <- igraph::induced_subgraph(ig, keep)
      m <- .metrics(sub, n0)
      data.frame(strategy = "targeted", fraction = fr, metric = metric_names,
                 value = m[metric_names], value_sd = NA_real_,
                 value_se = NA_real_, row.names = NULL,
                 stringsAsFactors = FALSE)
    })
    curve <- do.call(rbind, rows)

  } else {  # random
    .ggnv_local_seed(seed)
    rows <- lapply(c(0, fractions), function(fr) {
      k <- round(n0 * fr)
      if (k == 0L) {
        m <- .metrics(ig, n0)
        return(data.frame(strategy = "random", fraction = 0,
                          metric = metric_names, value = m[metric_names],
                          value_sd = 0, value_se = 0, row.names = NULL,
                          stringsAsFactors = FALSE))
      }
      boot <- vapply(seq_len(bootstrap), function(b) {
        rm_idx <- sample.int(n0, k)
        sub <- igraph::delete_vertices(ig, rm_idx)
        .metrics(sub, n0)[metric_names]
      }, numeric(length(metric_names)))
      mu <- rowMeans(boot)
      sdv <- apply(boot, 1L, stats::sd)
      data.frame(strategy = "random", fraction = fr, metric = metric_names,
                 value = mu, value_sd = sdv, value_se = sdv / sqrt(bootstrap),
                 row.names = NULL, stringsAsFactors = FALSE)
    })
    curve <- do.call(rbind, rows)
  }

  # ---- Schneider R-index (area under LCC_fraction curve) ------------------
  lcc <- curve[curve$metric == "LCC_fraction", c("fraction", "value")]
  lcc <- lcc[order(lcc$fraction), ]
  r_index <- mean(lcc$value)  # 1/N * sum s(Q), approximated on the grid
  robustness_index <- data.frame(strategy = strategy, R_index = r_index,
                                  stringsAsFactors = FALSE)

  # ---- optional plot ------------------------------------------------------
  p <- NULL
  if (isTRUE(plot) && strategy %in% c("random", "targeted")) {
    p <- ggnetview_perturbation_curve(curve, metric = "LCC_fraction")
  }

  list(curve = curve, robustness_index = robustness_index, plot = p)
}
