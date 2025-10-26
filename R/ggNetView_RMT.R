#' One-shot, non-interactive RMT threshold scan with auto selection
#'
#' @param mat  Numeric, symmetric, square matrix (e.g., correlation or weighted adjacency).
#' @param nr_thresholds Integer. Number of thresholds to scan (default 51).
#' @param interval Optional numeric length-2 vector [min, max] for scanning range over |mat| upper triangle (default auto).
#' @param unfold.method "gaussian" or "spline" (default "gaussian").
#' @param bandwidth Bandwidth for gaussian unfolding (passed to stats::density), default "nrd0".
#' @param nr.fit.points Integer, number of support points for spline unfolding (default 51).
#' @param discard.outliers Logical. Remove eigenvalue outliers via IQR before unfolding (default TRUE).
#' @param discard.zeros Logical. Drop all-zero rows/cols after thresholding (default TRUE).
#' @param min.mat.dim Integer. Early stop if effective dim < this (default 40).
#' @param max.ev.spacing Numeric. Cutoff for spacing tail when computing NNSD metrics (default 3).
#' @param save_plots Logical. Save PNGs (non-interactive) into out_dir (default FALSE).
#' @param out_dir Character. Output directory for plots (default "RMT_plots").
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @returns A list with: chosen_threshold, chosen_reason, tested_thresholds, scores (data.frame),
#'         unfolded (last-step unfolding), meta (matrix info & params), plots (file paths if saved).
#' @export
#'
#' @examples NULL
#'
#'
#' # m <- cor(scale(matrix(rnorm(40000), 200, 200)))
#' # res <- ggNetView_RMT(m, save_plots=TRUE)
#' # res$chosen_threshold
ggNetView_RMT <- function(
    mat,
    nr_thresholds = 51,
    interval = NULL,
    unfold.method = c("gaussian", "spline"),
    bandwidth = "nrd0",
    nr.fit.points = 51,
    discard.outliers = TRUE,
    discard.zeros = TRUE,
    min.mat.dim = 40,
    max.ev.spacing = 3,
    save_plots = FALSE,
    out_dir = "RMT_plots",
    verbose = TRUE
) {
  ## ---------- validation ----------
  stopifnot(is.matrix(mat))
  if (!is.numeric(mat)) stop("'mat' must be numeric.")
  if (nrow(mat) != ncol(mat)) stop("'mat' must be square.")
  if (!isSymmetric(mat)) stop("'mat' must be symmetric (real-valued).")
  unfold.method <- match.arg(unfold.method)

  N <- nrow(mat)
  if (verbose) message(sprintf("[Info] Matrix dimension: %d x %d", N, N))
  if (N < 100) warning("Matrix is relatively small (<100): RMT statistics may be unstable.")

  nz0  <- sum(mat != 0)
  sprs <- signif(sum(mat == 0) / (N * N), 4)
  if (verbose) message(sprintf("[Info] #non-zeros: %d | Sparseness: %.4f", nz0, sprs))

  ## ---------- threshold range ----------
  ut_abs <- abs(mat[upper.tri(mat, diag = FALSE)])
  min_cell <- min(ut_abs)
  max_cell <- max(ut_abs)
  if (is.null(interval)) {
    thresholds <- seq(min_cell, max_cell, length.out = nr_thresholds)
  } else {
    if (!is.numeric(interval) || length(interval) != 2)
      stop("'interval' must be a numeric vector of length 2.")
    if ((min(interval) < min_cell && min(interval) != 0) || max(interval) > max_cell) {
      stop(sprintf("'interval' must lie within the absolute upper-triangle range: [%.4g, %.4g]",
                   min_cell, max_cell))
    }
    thresholds <- seq(min(interval), max(interval), length.out = nr_thresholds)
  }

  ## ---------- plotting helpers (ggplot2, non-interactive PNG) ----------
  plot_files <- list()
  if (save_plots) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    gg_png <- function(p, fname, width=10, height=8, dpi=150) {
      ggplot2::ggsave(filename = file.path(out_dir, fname), plot = p,
                      width = width, height = height, dpi = dpi, units = "in")
    }
  }

  ## ---------- internal helpers ----------
  rank_qr <- function(A, tol = NULL) qr(A, tol = tol)$rank

  remove_outliers_IQR <- function(x, factor = 1.5) {
    q25 <- stats::quantile(x, 0.25, names = FALSE)
    q75 <- stats::quantile(x, 0.75, names = FALSE)
    iqr <- q75 - q25
    x[x >= (q25 - factor * iqr) & x <= (q75 + factor * iqr)]
  }

  unfold_gauss <- function(ev, bw = "nrd0") {
    ev <- sort(ev) / max(abs(ev))
    dens <- stats::density(ev, bw = bw, adjust = 1, n = 512, kernel = "gaussian")
    mids <- ev[-length(ev)] + 0.5 * diff(ev)
    scale_y <- stats::approx(dens$x, dens$y, xout = mids, rule = 2)$y
    spacing <- diff(ev)
    spacing <- spacing * scale_y
    spacing <- spacing / mean(spacing)
    uf_ev <- cumsum(spacing)
    uf_ev <- uf_ev / max(uf_ev)
    list(eigenvalues = ev, unfolded.ev = uf_ev, ev.spacing = spacing)
  }

  unfold_spline <- function(ev, nfit = 51) {
    ev <- sort(ev) / max(abs(ev))
    nfit <- min(nfit, floor(0.75 * length(ev)))
    cdf  <- stats::ecdf(ev)
    support <- seq(min(ev), max(ev), length.out = nfit)
    sf  <- stats::splinefun(support, cdf(support), method = "hyman")
    uf_ev <- sf(ev)
    spacing <- diff(uf_ev)
    spacing <- spacing / mean(spacing)
    list(eigenvalues = ev, unfolded.ev = uf_ev, ev.spacing = spacing)
  }

  ev_unfold <- function(A, method, bw, nfit, drop_outliers = TRUE) {
    ev <- eigen(A, only.values = TRUE)$values
    ev <- as.numeric(ev)
    if (drop_outliers) {
      ev0 <- ev
      ev  <- unique(remove_outliers_IQR(ev))
      removed <- length(ev0) - length(ev)
    } else {
      removed <- 0
    }
    u <- if (method == "gaussian") unfold_gauss(ev, bw) else unfold_spline(ev, nfit)
    u$nr.outliers.removed <- removed
    u
  }

  wigner_pdf <- function(x) (pi/2) * x * exp(-pi * x^2 / 4)
  exp_pdf    <- function(x) exp(-x)

  kld_discrete <- function(obs_pdf, exp_pdf) {
    obs_pdf <- pmax(obs_pdf, 0)
    exp_pdf <- pmax(exp_pdf, 0)
    s1 <- sum(obs_pdf); s2 <- sum(exp_pdf)
    if (s1 == 0 || s2 == 0) return(Inf)
    p <- obs_pdf / s1
    q <- exp_pdf / s2
    idx <- which(p > 0 & q > 0)
    sum(p[idx] * log(p[idx] / q[idx]))
  }

  nnsd_metrics <- function(sp, max_cut = 3, nbins = 51) {
    sp <- sp[is.finite(sp) & sp > 0]
    if (!is.null(max_cut)) sp <- sp[sp <= max_cut]
    if (length(sp) < 10) return(list(
      N = length(sp), kld_exp = NA, kld_wig = NA, ks_p = NA, sse_exp = NA, frac_small = NA
    ))
    brks <- seq(min(sp), max(sp), length.out = nbins)
    h <- hist(sp, breaks = brks, plot = FALSE)
    mids <- h$mids
    obs  <- h$density
    pdf_exp <- exp_pdf(mids)
    pdf_wig <- wigner_pdf(mids)
    kld_exp <- kld_discrete(obs, pdf_exp)
    kld_wig <- kld_discrete(obs, pdf_wig)
    ks_p <- tryCatch({
      suppressWarnings(stats::ks.test(unique(sp), "pexp", 1)$p.value)
    }, error = function(e) NA_real_)
    sse_exp <- {
      dens <- stats::density(sp, bw = "nrd0", n = 512)
      x <- seq(min(sp), max(sp), length.out = 1000)
      y_obs <- stats::approx(dens$x, dens$y, xout = x, rule = 2)$y
      A <- exp(-min(sp)) - exp(-max(sp))
      Nseg <- 20
      xs <- numeric(Nseg + 1)
      xs[1] <- min(sp)
      for (i in 1:Nseg) xs[i+1] <- -log(exp(-xs[i]) - A/Nseg)
      trap_int <- function(xx, yy) {
        idx <- 2:length(xx)
        as.double(t(xx[idx] - xx[idx-1]) %*% (yy[idx] + yy[idx-1])) / 2
      }
      area_obs <- numeric(Nseg)
      for (i in 1:Nseg) {
        xsec <- x[x >= xs[i] & x <= xs[i+1]]
        xsec <- unique(c(xs[i], xsec, xs[i+1]))
        ysec <- stats::approx(x, y_obs, xout = xsec, rule = 2)$y
        area_obs[i] <- trap_int(xsec, ysec)
      }
      sum((area_obs - A / Nseg)^2)
    }
    eps <- max_cut / 1000
    frac_small <- mean(sp < eps)
    list(N = length(sp), kld_exp = kld_exp, kld_wig = kld_wig,
         ks_p = ks_p, sse_exp = sse_exp, frac_small = frac_small)
  }

  effective_mat <- function(A, tol = 0) {
    dg <- diag(A); diag(A) <- 0
    keep_row <- rowSums(abs(A) > tol) > 0
    keep_col <- colSums(abs(A) > tol) > 0
    idx <- which(keep_row & keep_col)
    if (length(idx) == 0) {
      B <- matrix(0, 0, 0)
    } else {
      B <- A[idx, idx, drop = FALSE]
      diag(B) <- dg[idx]
    }
    B
  }

  denoise <- function(A, thr, keep_diag = TRUE) {
    if (keep_diag) dg <- diag(A)
    A[abs(A) < abs(thr)] <- 0
    if (keep_diag) diag(A) <- dg
    A
  }

  ## ---------- main scan loop (non-interactive) ----------
  res_tbl <- vector("list", length(thresholds))
  last_unfold <- NULL

  for (i in seq_along(thresholds)) {
    thr <- thresholds[i]
    if (verbose) message(sprintf("  [Scan] %3d/%d | threshold = %.4g", i, length(thresholds), thr))

    A <- denoise(mat, thr, keep_diag = TRUE)
    effA <- if (discard.zeros) effective_mat(A) else A

    effN <- nrow(effA)
    if (verbose) message(sprintf("         Effective dimension: %d", effN))
    if (effN == 0 || effN < min.mat.dim) {
      if (verbose) message("         Too small. Stop scanning.")
      res_tbl <- res_tbl[seq_len(i)]
      thresholds <- thresholds[seq_len(i)]
      break
    }

    u <- ev_unfold(effA, method = unfold.method, bw = bandwidth,
                   nfit = nr.fit.points, drop_outliers = discard.outliers)
    sp <- u$ev.spacing
    met <- nnsd_metrics(sp, max_cut = max.ev.spacing, nbins = 51)

    nr_zero <- sum(A == 0)
    ev_vals <- eigen(effA, only.values = TRUE)$values
    ev_r <- round(ev_vals, 8)
    tbl <- table(ev_r)
    uniq_ev <- length(tbl)
    max_mult <- max(tbl)

    res_tbl[[i]] <- data.frame(
      threshold = thr,
      eff_dim = effN,
      nr_zeros = nr_zero,
      nr_spacings = met$N,
      kld_exp = met$kld_exp,
      kld_wig = met$kld_wig,
      ks_p = met$ks_p,
      sse_exp = met$sse_exp,
      frac_small = met$frac_small,
      nr_unique_ev = uniq_ev,
      max_ev_mult = max_mult,
      stringsAsFactors = FALSE
    )

    if (save_plots) {
      # Per-threshold NNSD with theoretical curves (ggplot2)
      df_sp <- data.frame(spacing = sp[sp > 0 & sp <= max.ev.spacing])
      rng <- range(df_sp$spacing)
      p <- ggplot2::ggplot(df_sp, ggplot2::aes(x = spacing)) +
        ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)),
                                bins = 51, fill = "darkolivegreen3", color = "grey30") +
        ggplot2::stat_function(fun = exp_pdf, linewidth = 1, linetype = "solid") +
        ggplot2::stat_function(fun = wigner_pdf, linewidth = 1, linetype = "dashed") +
        ggplot2::labs(title = sprintf("NNSD (threshold = %.4g)", thr),
                      x = "eigenvalue spacing", y = "density",
                      subtitle = "solid: Exponential | dashed: Wigner surmise") +
        ggplot2::coord_cartesian(xlim = rng)
      gg_png(p, sprintf("NNSD_%03d.png", i))
    }

    last_unfold <- u
  }

  score_df <- do.call(rbind, res_tbl)
  rownames(score_df) <- NULL
  if (nrow(score_df) == 0) stop("No valid threshold result (matrix may be too sparse or 'min.mat.dim' too large).")

  ## ---------- auto selection (no interaction) ----------
  scl <- function(x) if (all(is.na(x))) rep(0, length(x)) else as.numeric(scale(x))
  w1 <- 0.35; w2 <- 0.25; w3 <- 0.20; w4 <- 0.10; w5 <- 0.10
  score <- w1 * scl(score_df$ks_p) +
    w2 * scl(-score_df$kld_exp) +
    w3 * scl(-score_df$sse_exp) +
    w4 * scl(-score_df$frac_small) +
    w5 * scl(score_df$eff_dim)
  score_df$auto_score <- score

  cand1_idx <- which(score_df$kld_exp == min(score_df$kld_exp, na.rm = TRUE))
  cand1_idx <- cand1_idx[which.max(score_df$eff_dim[cand1_idx])]

  best_idx <- which.max(score_df$auto_score)
  chosen_idx <- if (abs(score_df$threshold[best_idx] - score_df$threshold[cand1_idx]) < 1e-8) cand1_idx else best_idx
  chosen_threshold <- score_df$threshold[chosen_idx]
  chosen_reason <- sprintf(
    "Highest composite score (weights: ks_p %.2f, KLD_exp %.2f, SSE %.2f, small-spacing ratio %.2f, effective dimension %.2f).",
    w1, w2, w3, w4, w5
  )

  ## ---------- summary plots (ggplot2) ----------
  if (save_plots) {
    df <- score_df

    p1 <- ggplot2::ggplot(df, ggplot2::aes(threshold, ks_p)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = chosen_threshold, linetype = 2) +
      ggplot2::labs(title = "KS p-value vs threshold", x = "threshold", y = "p-value")
    gg_png(p1, "KS_p_vs_threshold.png", width=12, height=8, dpi=160)

    p2 <- ggplot2::ggplot(df, ggplot2::aes(threshold, sse_exp)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = chosen_threshold, linetype = 2) +
      ggplot2::labs(title = "SSE (NNSD vs Exponential) vs threshold", x = "threshold", y = "SSE")
    gg_png(p2, "SSE_vs_threshold.png", width=12, height=8, dpi=160)

    df_long <- data.frame(
      threshold = rep(df$threshold, 2),
      value = c(df$kld_wig, df$kld_exp),
      curve = rep(c("KLD(Wigner)", "KLD(Exponential)"), each = nrow(df))
    )
    p3 <- ggplot2::ggplot(df_long, ggplot2::aes(threshold, value, color = curve)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = chosen_threshold, linetype = 2) +
      ggplot2::labs(title = "KLD to Wigner / Exponential vs threshold", x = "threshold", y = "KLD") +
      ggplot2::scale_color_manual(values = c("steelblue4", "firebrick"))
    gg_png(p3, "KLD_vs_threshold.png", width=12, height=8, dpi=160)

    p4 <- ggplot2::ggplot(df, ggplot2::aes(threshold, auto_score)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = chosen_threshold, linetype = 2) +
      ggplot2::labs(title = "Auto score vs threshold (higher is better)", x = "threshold", y = "score")
    gg_png(p4, "AutoScore_vs_threshold.png", width=12, height=8, dpi=160)

    plot_files <- list(
      nnsd_each = list.files(out_dir, pattern = "^NNSD_\\d+\\.png$", full.names = TRUE),
      ks = file.path(out_dir, "KS_p_vs_threshold.png"),
      sse = file.path(out_dir, "SSE_vs_threshold.png"),
      kld = file.path(out_dir, "KLD_vs_threshold.png"),
      score = file.path(out_dir, "AutoScore_vs_threshold.png")
    )
  }

  ## ---------- return ----------
  list(
    chosen_threshold = chosen_threshold,
    chosen_reason    = chosen_reason,
    tested_thresholds = thresholds,
    scores = score_df,
    unfolded = last_unfold,
    meta = list(
      N = N,
      sparseness = sprs,
      nz0 = nz0,
      params = list(
        nr_thresholds = nr_thresholds,
        interval = if (is.null(interval)) c(min_cell, max_cell) else interval,
        unfold.method = unfold.method,
        bandwidth = bandwidth,
        nr.fit.points = nr.fit.points,
        discard.outliers = discard.outliers,
        discard.zeros = discard.zeros,
        min.mat.dim = min.mat.dim,
        max.ev.spacing = max.ev.spacing
      )
    ),
    plots = plot_files
  )
}
