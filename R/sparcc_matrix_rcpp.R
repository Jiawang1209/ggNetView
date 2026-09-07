#' SparCC correlation matrix (Rcpp implementation)
#'
#' Fast C++ implementation of SparCC for inferring correlations from
#' compositional count data. Input: samples x taxa matrix.
#' Output: taxa x taxa correlation matrix.
#'
#' @useDynLib ggNetView, .registration = TRUE
#'
#' @param data Numeric matrix. Rows = samples, columns = taxa (ASVs/OTUs).
#' @param iter Integer. Number of outer iterations (Dirichlet resampling).
#'   Default 20.
#' @param inner_iter Integer. Max inner iterations for pair exclusion.
#'   Default 10.
#' @param th Numeric. Threshold for excluding highly correlated pairs.
#'   Default 0.1.
#' @param nthreads Integer. Number of OpenMP threads for the parallel C++ SparCC
#'   inner loop (0 = use default). Only effective when OpenMP is available.
#'
#' @return Numeric matrix of taxa x taxa correlations (median across iterations).
#' @export
#'
#' @examples
#' \dontrun{
#' # data: samples x taxa count matrix
#' cor_mat <- sparcc_matrix_rcpp(asv_mat, iter = 20, inner_iter = 10, th = 0.1)
#' # With 4 threads (when OpenMP available)
#' cor_mat <- sparcc_matrix_rcpp(asv_mat, iter = 20, nthreads = 4)
#' }
sparcc_matrix_rcpp <- function(data, iter = 20, inner_iter = 10, th = 0.1, nthreads = 0L) {
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("data must be a numeric matrix", call. = FALSE)
  }
  if (any(data < 0, na.rm = TRUE)) {
    stop("data contains negative values", call. = FALSE)
  }
  iter <- as.integer(iter)[1L]
  inner_iter <- as.integer(inner_iter)[1L]
  if (iter < 1L || inner_iter < 1L) {
    stop("iter and inner_iter must be positive integers", call. = FALSE)
  }
  th <- as.numeric(th)[1L]
  if (is.na(th) || th < 0) {
    stop("th must be a non-negative numeric", call. = FALSE)
  }
  nthreads <- as.integer(nthreads)[1L]
  if (is.na(nthreads) || nthreads < 0) {
    stop("nthreads must be a non-negative integer", call. = FALSE)
  }

  out <- sparcc_matrix_cpp(data, iter, inner_iter, th, nthreads)
  storage.mode(out) <- "double"
  dimnames(out) <- list(colnames(data), colnames(data))
  out
}


#' SparCC p-value matrix (Rcpp implementation)
#'
#' Fast C++-accelerated SparCC p-values via bootstrap and permutation.
#' Uses \code{sparcc_matrix_rcpp} internally for speed. Input: samples x taxa matrix.
#' Output: taxa x taxa p-value matrix.
#'
#' @param data Numeric matrix. Rows = samples, columns = taxa (ASVs/OTUs).
#' @param iter Integer. Number of outer iterations (Dirichlet resampling).
#'   Default 20.
#' @param inner_iter Integer. Max inner iterations for pair exclusion.
#'   Default 10.
#' @param th Numeric. Threshold for excluding highly correlated pairs.
#'   Default 0.1.
#' @param R Integer. Number of bootstrap/permutation replicates. Default 20.
#' @param ncpus Integer. Number of CPUs for parallel boot. Default 1.
#' @param seed Integer or NULL (default). If supplied, makes the bootstrap /
#'   permutation reproducible. When \code{ncpus > 1} the L'Ecuyer-CMRG generator
#'   is used so parallel workers get independent, deterministic RNG streams; the
#'   previous RNG kind is restored on exit.
#'
#' @return Numeric matrix of taxa x taxa p-values. Diagonal = 0.
#'   NaN indicates observed correlation outside bootstrap CI.
#' @export
#'
#' @examples
#' \dontrun{
#' p_mat <- sparcc_pvalue_rcpp(asv_mat, R = 20, ncpus = 4, seed = 1115)
#' }
sparcc_pvalue_rcpp <- function(data, iter = 20, inner_iter = 10, th = 0.1, R = 20, ncpus = 1, seed = NULL) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("package 'boot' is required for SparCC p-values", call. = FALSE)
  }
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("data must be a numeric matrix", call. = FALSE)
  }
  if (any(data < 0, na.rm = TRUE)) {
    stop("data contains negative values", call. = FALSE)
  }
  iter <- as.integer(iter)[1L]
  inner_iter <- as.integer(inner_iter)[1L]
  R <- as.integer(R)[1L]
  ncpus <- as.integer(ncpus)[1L]
  if (iter < 1L || inner_iter < 1L || R < 1L || ncpus < 1L) {
    stop("iter, inner_iter, R and ncpus must be positive integers", call. = FALSE)
  }
  th <- as.numeric(th)[1L]
  if (is.na(th) || th < 0) {
    stop("th must be a non-negative numeric", call. = FALSE)
  }

  # Reproducible parallel resampling: L'Ecuyer-CMRG gives each worker an
  # independent, deterministic stream. Restore the caller's RNG kind on exit.
  if (!is.null(seed)) {
    # Capture the caller's RNG kind *and* state up front; the handler restores
    # both when this function returns, so neither the L'Ecuyer switch below nor
    # the seed leaks into the user's session.
    .ggnv_local_seed(seed)
    if (ncpus > 1L) {
      RNGkind("L'Ecuyer-CMRG")
    }
    # Reseed under the (possibly switched) kind so the stream is the documented,
    # deterministic one -- identical to the previous behaviour.
    set.seed(as.integer(seed)[1L])
  }

  triu_local <- function(x) x[upper.tri(x)]
  triu2diag_local <- function(x, diagval = 0) {
    e <- length(x)
    n <- 0.5 * (sqrt(8 * e + 1) + 1)
    mat <- matrix(0, n, n)
    mat[upper.tri(mat)] <- x
    mat <- mat + t(mat)
    diag(mat) <- diagval
    mat
  }

  statisticboot <- function(dat, indices) {
    triu_local(sparcc_matrix_rcpp(dat[indices, , drop = FALSE], iter = iter, inner_iter = inner_iter, th = th))
  }
  statisticperm <- function(dat, indices) {
    permuted <- apply(dat[indices, , drop = FALSE], 2, sample)
    triu_local(sparcc_matrix_rcpp(permuted, iter = iter, inner_iter = inner_iter, th = th))
  }

  # multicore (fork) not available on Windows; use snow for cross-platform parallel
  parallel_type <- if (.Platform$OS.type == "windows") "snow" else "multicore"
  res <- boot::boot(data, statisticboot, R = R, parallel = parallel_type, ncpus = ncpus)
  null_av <- boot::boot(data, statisticperm, sim = "permutation", R = R, parallel = parallel_type, ncpus = ncpus)

  nparams <- ncol(res$t)
  tmeans <- colMeans(null_av$t)
  niters <- nrow(res$t)
  ind95 <- max(1, round(0.025 * niters)):round(0.975 * niters)
  boot_ord <- apply(res$t, 2, sort)
  boot_ord95 <- boot_ord[ind95, , drop = FALSE]
  outofrange <- vapply(
    seq_along(res$t0),
    function(i) {
      aitvar <- res$t0[i]
      rg <- range(boot_ord95[, i])
      rg[1] > aitvar || rg[2] < aitvar
    },
    logical(1)
  )
  bs_above <- vapply(
    seq_len(nparams),
    function(i) length(which(res$t[, i] > tmeans[i])),
    numeric(1)
  )
  is_above <- bs_above > res$R / 2
  pvals <- ifelse(is_above, 2 * (1 - bs_above / res$R), 2 * bs_above / res$R)
  pvals[pvals > 1] <- 1
  pvals[outofrange] <- NaN

  out <- triu2diag_local(pvals, diagval = 0)
  dimnames(out) <- list(colnames(data), colnames(data))
  out
}
