#' Build a consensus graph object from multiple adjacency matrices
#'
#' Aggregates two or more adjacency matrices -- typically produced by
#' different network inference methods (e.g. SparCC, SpiecEasi, WGCNA,
#' Hmisc, plain correlation) on the same set of features -- into a single
#' consensus network and returns a `tbl_graph` consistent with the rest
#' of the `build_graph_from_*()` family. Combining multiple methods is a
#' standard strategy for reducing method-specific bias in biological
#' network inference; see Aghayeva et al. (2024, CMiNet) and Chowdhury
#' et al. (2024, hybrid Bayesian / ML / network framework) for related
#' approaches.
#'
#' Four consensus strategies are supported through `method`. Three of
#' them (`"intersection"`, `"weighted_average"`, `"majority_vote"`)
#' operate directly on the input weights; the fourth (`"rank_fusion"`)
#' converts each method's adjacency to ranks first, which makes the
#' aggregation invariant to the very different weight scales produced
#' by different inference algorithms (correlation in [-1, 1], partial
#' correlation, TOM in [0, 1], ...). When `method = "rank_fusion"`,
#' three classical rank-aggregation algorithms are available via
#' `rank_fusion_algorithm`: Borda count (mean rank), Robust Rank
#' Aggregation (RRA, Kolde et al. 2012), and Reciprocal Rank Fusion
#' (RRF, Cormack et al. 2009).
#'
#' @section Method details:
#' \describe{
#'   \item{`"intersection"`}{An edge survives in the consensus iff it
#'     is non-zero (after `binarize`) in **every** input matrix.
#'     Consensus weight is the mean of the input weights for surviving
#'     edges. Strictest of the four strategies.}
#'   \item{`"weighted_average"`}{Each input matrix is min-max
#'     normalised to [0, 1] in absolute value (sign preserved), then
#'     combined as a weighted sum using `weights` (default uniform).
#'     Consensus weight is in [-1, 1] and inherits the sign of the
#'     methods that dominate the average.}
#'   \item{`"majority_vote"`}{Each input matrix is binarised, edges are
#'     summed, and an edge survives in the consensus iff the count is
#'     at least `min_methods`. Consensus weight is the mean of the
#'     input weights restricted to the methods that voted yes.}
#'   \item{`"rank_fusion"`}{Each input matrix's pairs are ranked by
#'     `|w|` (descending). Ranks are aggregated according to
#'     `rank_fusion_algorithm`, then min-max normalised to [-1, 1]
#'     (sign taken from the mean of the input weights). Most robust to
#'     across-method scale differences.}
#' }
#'
#' @section Rank fusion algorithms:
#' \describe{
#'   \item{`"borda"`}{Mean rank across methods. Lower mean rank = more
#'     supported edge. Simple, no assumptions, good baseline.}
#'   \item{`"rra"`}{Robust Rank Aggregation (Kolde et al. 2012,
#'     Bioinformatics). Computes a rho-score for each pair under the
#'     null hypothesis of uniformly distributed ranks; edges with very
#'     small rho consistently outrank random expectation. Uses the
#'     `RobustRankAggreg` package when available; otherwise falls back
#'     to an inline Beta-order-statistic implementation. Output score
#'     is `-log10(rho)` so larger = more confident.}
#'   \item{`"rrf"`}{Reciprocal Rank Fusion (Cormack et al. 2009).
#'     Per-edge score is `sum_m 1 / (k + rank_m)` where `k = rrf_k`
#'     (default 60, the information-retrieval convention). Robust to
#'     extreme rank values, no parameter tuning required.}
#' }
#'
#' @param adj_list Named or unnamed list of at least two square numeric
#'   adjacency matrices. All matrices must have row and column names
#'   identifying features; matrices need not share identical feature
#'   sets (see `node_handling`). Pre-thresholded inputs (zeros for
#'   non-significant pairs) are supported and recommended.
#' @param method Character. Top-level consensus strategy. One of
#'   `"rank_fusion"` (default), `"intersection"`, `"weighted_average"`,
#'   `"majority_vote"`. See \strong{Method details}.
#' @param rank_fusion_algorithm Character. Only used when `method =
#'   "rank_fusion"`. One of `"borda"` (default), `"rra"`, `"rrf"`.
#'   See \strong{Rank fusion algorithms}.
#' @param threshold Numeric or `NULL` (default). Optional final
#'   absolute-weight cutoff applied to the consensus matrix; entries
#'   with `abs(w) < threshold` are zeroed out.
#' @param binarize Character. Only used by `"intersection"` and
#'   `"majority_vote"`. One of `"none"` (default; treat any non-zero
#'   entry as an edge), `"threshold"` (apply `binarize_threshold`),
#'   `"topk"` (keep the global top-`binarize_topk` strongest edges per
#'   matrix).
#' @param binarize_threshold Numeric (default `0`). Per-method
#'   `|w| >= binarize_threshold` edge cut-off when `binarize =
#'   "threshold"`.
#' @param binarize_topk Integer or `NULL`. Per-method top-K cut-off
#'   when `binarize = "topk"`.
#' @param weights Numeric vector or `NULL` (default). Per-method
#'   weights for `"weighted_average"`. Length must match `length(
#'   adj_list)`; values are renormalised to sum to 1. `NULL` uses
#'   uniform weights.
#' @param min_methods Integer or `NULL` (default). Minimum number of
#'   methods that must support an edge for `"majority_vote"`. `NULL`
#'   defaults to a strict majority (`floor(M/2) + 1` where `M` is the
#'   number of methods).
#' @param rrf_k Numeric (default `60`). Smoothing constant for
#'   reciprocal rank fusion (`method = "rank_fusion",
#'   rank_fusion_algorithm = "rrf"`).
#' @param node_handling Character. How to align features that are not
#'   present in every input matrix. `"intersect"` (default) keeps only
#'   features common to all inputs; `"union"` keeps every feature
#'   appearing in any input and treats missing entries as zero.
#' @param module.method Character. Module detection method passed
#'   through to [build_graph_from_adj_mat()] when constructing the
#'   final `tbl_graph`. One of `"Fast_greedy"`, `"Walktrap"`,
#'   `"Edge_betweenness"`, `"Spinglass"`.
#' @param node_annotation Optional data frame attached as vertex
#'   metadata; first column must match feature names of the consensus
#'   matrix.
#' @param top_modules Integer (default `15`). Number of top-ranked
#'   modules to retain in the final `tbl_graph`; smaller modules are
#'   collapsed into `"Others"`.
#' @param seed Integer (default `1115`). Random seed.
#'
#' @returns A `tbl_graph` object whose schema matches the rest of the
#'   `build_graph_from_*()` family (vertex columns `name`,
#'   `modularity`, `modularity2`, `modularity3`, `Modularity`,
#'   `Degree`, `Strength`; edge columns `weight`, `correlation`,
#'   `corr_direction`). Plug straight into [ggNetView()] for
#'   visualisation.
#'
#' @references
#' Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). Robust rank
#' aggregation for gene list integration and meta-analysis.
#' *Bioinformatics*, 28(4), 573-580.
#'
#' Cormack, G. V., Clarke, C. L. A., & Buettcher, S. (2009). Reciprocal
#' rank fusion outperforms Condorcet and individual rank learning
#' methods. *Proceedings of the 32nd International ACM SIGIR
#' Conference*, 758-759.
#'
#' Aghayeva, R., et al. (2024). CMiNet: An R package and user-friendly
#' Shiny App for constructing consensus microbiome networks.
#'
#' Chowdhury, S., et al. (2024). A hybrid framework for disease
#' biomarker discovery in microbiome research combining Bayesian
#' networks, machine learning, and network-based methods.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 30
#' nm <- paste0("g", seq_len(n))
#'
#' # Three "methods": shared signal + per-method noise.
#' core <- matrix(0, n, n, dimnames = list(nm, nm))
#' for (i in 1:6) for (j in (i + 1):8) core[i, j] <- core[j, i] <- 0.8
#' jitter <- function() core + matrix(stats::rnorm(n * n, sd = 0.05), n, n,
#'                                    dimnames = list(nm, nm))
#'
#' adj_list <- list(method_A = jitter(), method_B = jitter(), method_C = jitter())
#'
#' # Borda rank fusion with a final |w| >= 0.3 cutoff:
#' obj <- build_graph_from_consensus(
#'   adj_list = adj_list,
#'   method = "rank_fusion",
#'   rank_fusion_algorithm = "borda",
#'   threshold = 0.3
#' )
#' obj
#' }
build_graph_from_consensus <- function(
  adj_list,
  method = c("rank_fusion", "intersection", "weighted_average", "majority_vote"),
  rank_fusion_algorithm = c("borda", "rra", "rrf"),
  threshold = NULL,
  binarize = c("none", "threshold", "topk"),
  binarize_threshold = 0,
  binarize_topk = NULL,
  weights = NULL,
  min_methods = NULL,
  rrf_k = 60,
  node_handling = c("intersect", "union"),
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
) {

  method <- match.arg(method)
  rank_fusion_algorithm <- match.arg(rank_fusion_algorithm)
  binarize <- match.arg(binarize)
  node_handling <- match.arg(node_handling)
  module.method <- match.arg(module.method)

  # ---- Input validation ----------------------------------------------------
  if (!is.list(adj_list) || length(adj_list) < 2L) {
    stop("`adj_list` must be a list of at least 2 adjacency matrices.", call. = FALSE)
  }
  for (i in seq_along(adj_list)) {
    m <- adj_list[[i]]
    if (is.data.frame(m)) m <- as.matrix(m)
    if (!is.matrix(m) || !is.numeric(m)) {
      stop(sprintf("`adj_list[[%d]]` must be a numeric matrix.", i), call. = FALSE)
    }
    if (nrow(m) != ncol(m)) {
      stop(sprintf("`adj_list[[%d]]` is not square (%d x %d).", i, nrow(m), ncol(m)),
           call. = FALSE)
    }
    if (is.null(rownames(m)) || is.null(colnames(m))) {
      stop(sprintf("`adj_list[[%d]]` must have row and column names.", i),
           call. = FALSE)
    }
    if (!identical(rownames(m), colnames(m))) {
      stop(sprintf("`adj_list[[%d]]` must have identical row and column names.", i),
           call. = FALSE)
    }
    adj_list[[i]] <- m
  }

  .ggnv_local_seed(seed)

  # ---- Step 1: align matrices to a common node set -------------------------
  mats <- .consensus_align(adj_list, node_handling)
  M <- length(mats)
  N <- nrow(mats[[1L]])

  if (N < 2L) {
    stop(sprintf("After `node_handling = '%s'` alignment, fewer than 2 features remain.",
                 node_handling), call. = FALSE)
  }

  # ---- Step 2: compute consensus matrix ------------------------------------
  cons_mat <- switch(
    method,
    intersection     = .consensus_intersection(mats, binarize, binarize_threshold, binarize_topk),
    weighted_average = .consensus_weighted_average(mats, weights),
    majority_vote    = .consensus_majority_vote(mats, binarize, binarize_threshold, binarize_topk, min_methods),
    rank_fusion      = .consensus_rank_fusion(mats, rank_fusion_algorithm, rrf_k)
  )

  # Force symmetry (numeric drift from float arithmetic can break
  # downstream igraph assumptions about undirected adjacency).
  cons_mat <- (cons_mat + t(cons_mat)) / 2
  diag(cons_mat) <- 0

  # ---- Step 3: optional final threshold ------------------------------------
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L) {
      stop("`threshold` must be a single numeric value or NULL.", call. = FALSE)
    }
    cons_mat[abs(cons_mat) < threshold] <- 0
  }

  # ---- Step 4: delegate to build_graph_from_adj_mat() ----------------------
  # Routing through the existing builder guarantees the output schema
  # (Modularity / Degree / Strength etc., weight / correlation /
  # corr_direction) matches every other ggNetView entry point, so no
  # special-casing is needed downstream.
  build_graph_from_adj_mat(
    adjacency_matrix = cons_mat,
    module.method    = module.method,
    node_annotation  = node_annotation,
    top_modules      = top_modules,
    seed             = seed
  )
}


# ===========================================================================
# Helpers (not exported)
# ===========================================================================

#' Align a list of adjacency matrices to a common node set
#' @noRd
.consensus_align <- function(adj_list, node_handling) {
  all_nodes <- lapply(adj_list, rownames)

  if (node_handling == "intersect") {
    common <- Reduce(intersect, all_nodes)
    if (length(common) < 2L) {
      stop(sprintf(
        "Intersection of feature sets has fewer than 2 features (n = %d). Try `node_handling = 'union'` or check input row/column names.",
        length(common)
      ), call. = FALSE)
    }
    lapply(adj_list, function(m) m[common, common, drop = FALSE])
  } else {
    # union: every feature appearing in any input; missing entries -> 0
    common <- Reduce(union, all_nodes)
    lapply(adj_list, function(m) {
      out <- matrix(0, nrow = length(common), ncol = length(common),
                    dimnames = list(common, common))
      shared <- intersect(common, rownames(m))
      out[shared, shared] <- m[shared, shared]
      out
    })
  }
}


#' Binarise an adjacency matrix
#' @noRd
.consensus_binarize <- function(mat, binarize, binarize_threshold, binarize_topk) {
  if (binarize == "none") {
    out <- (mat != 0) * 1
  } else if (binarize == "threshold") {
    out <- (abs(mat) >= binarize_threshold) * 1
  } else {
    # topk: keep the K strongest pairs in the upper triangle, mirror to
    # both halves so the binarised matrix stays symmetric.
    if (is.null(binarize_topk) || !is.numeric(binarize_topk) || binarize_topk < 1L) {
      stop("`binarize = 'topk'` requires `binarize_topk >= 1`.", call. = FALSE)
    }
    diag_orig <- diag(mat)
    m <- mat
    diag(m) <- 0
    upper_vals <- abs(m[upper.tri(m)])
    K <- min(as.integer(binarize_topk), length(upper_vals))
    cutoff <- sort(upper_vals, decreasing = TRUE)[K]
    upper_keep <- (abs(m) >= cutoff) & upper.tri(m)
    out <- matrix(0, nrow = nrow(m), ncol = ncol(m),
                  dimnames = dimnames(m))
    out[upper_keep] <- 1
    out <- out + t(out)
  }
  diag(out) <- 0
  out
}


#' Intersection consensus
#' @noRd
.consensus_intersection <- function(mats, binarize, binarize_threshold, binarize_topk) {
  bin_mats <- lapply(mats, .consensus_binarize, binarize, binarize_threshold, binarize_topk)
  # Element-wise AND across methods: edge kept iff every method has it.
  intersect_mask <- Reduce(`*`, bin_mats)
  # Consensus weight = mean signed weight across methods, restricted to
  # surviving edges.
  mean_signed <- Reduce(`+`, mats) / length(mats)
  intersect_mask * mean_signed
}


#' Weighted-average consensus
#' @noRd
.consensus_weighted_average <- function(mats, weights) {
  M <- length(mats)
  if (is.null(weights)) {
    weights <- rep(1, M)
  }
  if (!is.numeric(weights) || length(weights) != M) {
    stop(sprintf(
      "`weights` must be a numeric vector of length %d (one per method); got length %d.",
      M, length(weights)
    ), call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("`weights` must be non-negative.", call. = FALSE)
  }
  if (sum(weights) == 0) {
    stop("`weights` must contain at least one positive value.", call. = FALSE)
  }
  weights <- weights / sum(weights)

  # Per-method min-max normalise abs-weights to [0, 1], preserving sign.
  norm_mats <- lapply(mats, function(m) {
    a <- abs(m)
    mx <- max(a)
    if (mx > 0) sign(m) * a / mx else m
  })

  out <- matrix(0, nrow = nrow(mats[[1L]]), ncol = ncol(mats[[1L]]),
                dimnames = dimnames(mats[[1L]]))
  for (i in seq_len(M)) {
    out <- out + weights[i] * norm_mats[[i]]
  }
  out
}


#' Majority-vote consensus
#' @noRd
.consensus_majority_vote <- function(mats, binarize, binarize_threshold, binarize_topk,
                                     min_methods) {
  M <- length(mats)
  if (is.null(min_methods)) {
    min_methods <- floor(M / 2) + 1L  # strict majority
  }
  if (!is.numeric(min_methods) || min_methods < 1 || min_methods > M) {
    stop(sprintf("`min_methods` must be an integer in [1, %d]; got %s.",
                 M, format(min_methods)), call. = FALSE)
  }
  bin_mats <- lapply(mats, .consensus_binarize, binarize, binarize_threshold, binarize_topk)
  vote_count <- Reduce(`+`, bin_mats)
  keep <- (vote_count >= min_methods) * 1
  # Mean signed weight of the methods that voted yes (= mean over yes-voters).
  # Compute as: sum_m bin_m * mat_m  /  vote_count   (element-wise)
  numer <- 0
  for (i in seq_len(M)) {
    numer <- numer + bin_mats[[i]] * mats[[i]]
  }
  denom <- vote_count
  denom[denom == 0] <- 1L  # avoid 0/0; numerator is also 0 there
  keep * (numer / denom)
}


#' Rank-fusion consensus (Borda / RRA / RRF)
#' @noRd
.consensus_rank_fusion <- function(mats, algorithm, rrf_k) {
  M <- length(mats)
  N <- nrow(mats[[1L]])

  # Operate on the strict upper triangle so each undirected pair is
  # ranked exactly once.
  upper_mask <- upper.tri(mats[[1L]])
  K <- sum(upper_mask)

  # Per-method rank of |w|, descending. Ties handled by "average".
  # Pairs with |w| = 0 land at the worst rank (largest number).
  rank_mat <- matrix(NA_real_, nrow = K, ncol = M)
  for (m in seq_len(M)) {
    w <- abs(mats[[m]][upper_mask])
    # rank() ranks ascending; we negate so largest |w| -> rank 1.
    rank_mat[, m] <- rank(-w, ties.method = "average", na.last = TRUE)
  }

  if (algorithm == "borda") {
    # Mean rank across methods. Lower = better. Convert to "score where
    # higher = better" by negating.
    score <- -rowMeans(rank_mat)

  } else if (algorithm == "rrf") {
    # Reciprocal rank fusion: sum_m 1 / (k + rank_m). Higher = better.
    if (!is.numeric(rrf_k) || length(rrf_k) != 1L || rrf_k <= 0) {
      stop("`rrf_k` must be a positive numeric scalar.", call. = FALSE)
    }
    score <- rowSums(1 / (rrf_k + rank_mat))

  } else {
    # algorithm == "rra"
    # Robust Rank Aggregation (Kolde et al. 2012). Use the official
    # implementation when available; otherwise fall back to an inline
    # Beta-order-statistic computation.
    norm_rank <- rank_mat / K  # normalise to [0, 1]
    if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
      rho <- vapply(seq_len(K), function(k) {
        RobustRankAggreg::rhoScores(sort(norm_rank[k, ]))
      }, numeric(1))
    } else {
      # Inline: rho(r) = min over j of pbeta(r_(j), j, M - j + 1),
      # corrected by * M for the multiple-position test.
      rho <- vapply(seq_len(K), function(k) {
        sorted_r <- sort(norm_rank[k, ])
        candidates <- vapply(seq_len(M), function(j) {
          stats::pbeta(sorted_r[j], j, M - j + 1)
        }, numeric(1))
        min(min(candidates) * M, 1)
      }, numeric(1))
    }
    score <- -log10(rho + 1e-300)  # higher = more confident
  }

  # Min-max normalise score to [0, 1], then re-attach the sign of the
  # mean correlation across methods. Pairs that everyone gave |w| = 0
  # can produce score = NaN/-Inf; clamp those to the smallest finite
  # score, or to 0 if no scores are finite (otherwise `min(numeric(0))`
  # is `Inf` and every score gets stamped with `Inf`).
  .finite_min <- if (any(is.finite(score))) min(score[is.finite(score)]) else 0
  score[!is.finite(score)] <- .finite_min
  score_range <- range(score, finite = TRUE)
  if (diff(score_range) > 0) {
    score_norm <- (score - score_range[1L]) / diff(score_range)
  } else {
    score_norm <- rep(0, length(score))
  }

  # Sign from the mean signed weight of each pair.
  mean_signed <- Reduce(`+`, mats) / M
  signed_pair <- sign(mean_signed[upper_mask])
  # Pairs with no signal across methods get sign 0 -> consensus weight 0.
  signed_score <- signed_pair * score_norm

  # Reshape back to symmetric N x N matrix.
  out <- matrix(0, nrow = N, ncol = N,
                dimnames = dimnames(mats[[1L]]))
  out[upper_mask] <- signed_score
  out + t(out)
}
