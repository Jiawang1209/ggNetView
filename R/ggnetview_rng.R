# ---------------------------------------------------------------------------
# Internal helpers for RNG hygiene and multiple-testing correction.
#
# Both address contracts that ggNetView's exported functions are expected to
# honour but previously did not:
#
#   * `.ggnv_local_seed()` seeds the RNG for the duration of ONE call and
#     restores the caller's stream on exit, so `seed = ` no longer changes the
#     state of the user's session (CRAN policy, and a prerequisite for the
#     package's reproducibility claim).
#   * `.ggnv_adjust_p_matrix()` corrects a symmetric correlation p-value matrix
#     over its n(n-1)/2 unique off-diagonal tests instead of over all n^2 cells.
# ---------------------------------------------------------------------------

#' Seed the RNG locally, restoring the caller's stream on exit
#'
#' Drop-in replacement for a bare `set.seed(seed)` at the top of a function.
#' The seed takes effect immediately and for the rest of the calling function,
#' exactly as `set.seed()` would, but `.Random.seed` is restored when that
#' function returns -- including when it exits via an error. Numerical results
#' are therefore unchanged; only the leak into the caller's session is removed.
#'
#' @param seed Single value passed to [set.seed()]. `NULL` is a no-op.
#' @param envir Frame whose exit should restore the stream; defaults to the
#'   caller, which is what every call site wants.
#' @return Invisibly `TRUE` when a seed was set, `FALSE` when `seed` was `NULL`.
#' @noRd
.ggnv_local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(FALSE))

  # Capture the generator kind as well as the state. `.Random.seed[1]` encodes
  # the kind, so a caller that switches RNGkind() (e.g. to L'Ecuyer-CMRG for
  # parallel streams) would otherwise have that switch silently restored along
  # with the state. Restoring the kind first and the state second puts both
  # back exactly as they were: RNGkind() itself reinitialises .Random.seed, so
  # the explicit assignment has to come after it.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

  restore <- bquote({
    RNGkind(kind = .(old_kind[1L]),
            normal.kind = .(old_kind[2L]),
            sample.kind = .(old_kind[3L]))
    if (.(had_seed)) {
      assign(".Random.seed", .(old_seed), envir = globalenv())
    } else {
      # The session had not drawn a random number yet: leave it that way.
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    }
  })

  # Register the restore handler on the CALLING function's frame, so it fires
  # when that function returns rather than when this helper returns.
  do.call(base::on.exit, list(restore, add = TRUE, after = TRUE), envir = envir)

  set.seed(as.integer(seed)[1L])
  invisible(TRUE)
}

#' Adjust a symmetric p-value matrix over its unique off-diagonal tests
#'
#' A correlation p-value matrix holds each pairwise test twice (once per
#' triangle) plus an uninformative diagonal (p = 0 for `WGCNA::corAndPvalue()`
#' and `psych::corr.test()`, `NA` for `Hmisc::rcorr()`). Passing all n^2 cells
#' to [stats::p.adjust()] inflates the test count and, because the n diagonal
#' zeros occupy the lowest ranks, makes Benjamini-Hochberg anti-conservative
#' for the most significant edges. Correct over the upper triangle only, then
#' mirror.
#'
#' @param p_mat Symmetric numeric matrix of p-values.
#' @param proc_method Method passed to [stats::p.adjust()].
#' @return Matrix of the same shape and dimnames; the diagonal is set to 1 so
#'   that self-loops can never pass a `p <= threshold` test.
#' @noRd
.ggnv_adjust_p_matrix <- function(p_mat, proc_method) {
  p_mat <- as.matrix(p_mat)
  n <- nrow(p_mat)
  if (n != ncol(p_mat)) {
    stop("`p_mat` must be square.", call. = FALSE)
  }
  if (n < 2L) {
    out <- p_mat
    if (n == 1L) out[1L, 1L] <- 1
    return(out)
  }

  ut  <- upper.tri(p_mat, diag = FALSE)
  adj <- p_mat
  adj[ut] <- stats::p.adjust(p_mat[ut], method = proc_method)
  adj[lower.tri(adj, diag = FALSE)] <- t(adj)[lower.tri(adj, diag = FALSE)]
  diag(adj) <- 1
  dimnames(adj) <- dimnames(p_mat)
  adj
}
