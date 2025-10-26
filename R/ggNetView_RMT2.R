#' ggNetView_RMT2: 一次性、非交互的随机矩阵阈值扫描与自动选阈值
#'
#' @param mat 数值型、实对称矩阵（如相关矩阵、邻接矩阵）
#' @param nr_thresholds 整数。阈值扫描的格点数，默认 51
#' @param interval 长度为 2 的数值向量。阈值搜索区间（基于 |mat| 的上三角），默认自动[min, max]
#' @param unfold.method "gaussian" 或 "spline"，默认 "gaussian"
#' @param bandwidth 仅对 gaussian 展开有效。传递给 density() 的 bw，默认 "nrd0"
#' @param nr.fit.points 仅对 spline 展开有效。样条拟合采样点数，默认 51
#' @param discard.outliers 是否在展开前剔除特征值异常点（IQR：1.5），默认 TRUE
#' @param discard.zeros 在每个阈值下是否丢弃全零的行列（保持对角），默认 TRUE
#' @param min.mat.dim 有效矩阵维度（非零行列）低于此值提前停止扫描，默认 40
#' @param max.ev.spacing 计算 NNSD 时的最大间距截断（去除极大 spacing 的长尾），默认 3
#' @param save_plots 是否保存诊断图，默认 FALSE（非交互；保存为 png）
#' @param out_dir 图像输出目录（当 save_plots=TRUE 时必须有写权限），默认 "RMT_plots"
#' @param verbose 是否打印简要进度信息，默认 TRUE
#'
#'
#' @examples NULL
#' # m <- cor(scale(matrix(rnorm(40000), 200, 200)))
#' # res <- ggNetView_RMT(m, save_plots=TRUE)
#' # res$chosen_threshold
#'
ggNetView_RMT2 <- function(
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
  ## ------------------------- 基本校验 -------------------------
  stopifnot(is.matrix(mat))
  if (!is.numeric(mat)) stop("'mat' 必须是数值矩阵。")
  if (nrow(mat) != ncol(mat)) stop("'mat' 必须是方阵。")
  if (!isSymmetric(mat)) stop("'mat' 必须是实对称矩阵（symmetric=TRUE）。")
  unfold.method <- match.arg(unfold.method)

  N <- nrow(mat)
  if (verbose) message(sprintf("[Info] 输入矩阵维度: %d x %d", N, N))
  if (N < 100) warning("矩阵较小（<100），RMT 统计可能不稳。")

  ## 矩阵稀疏度与非零计数（参考）
  nz0  <- sum(mat != 0)
  sprs <- signif(sum(mat == 0) / (N * N), 4)
  if (verbose) message(sprintf("[Info] 非零元素数: %d | Sparseness: %.4f", nz0, sprs))

  ## 阈值搜索区间（基于上三角绝对值，不含对角）
  ut_abs <- abs(mat[upper.tri(mat, diag = FALSE)])
  min_cell <- min(ut_abs)
  max_cell <- max(ut_abs)
  if (is.null(interval)) {
    thresholds <- seq(min_cell, max_cell, length.out = nr_thresholds)
  } else {
    if (!is.numeric(interval) || length(interval) != 2)
      stop("'interval' 必须是长度为 2 的数值向量。")
    if ((min(interval) < min_cell && min(interval) != 0) || max(interval) > max_cell) {
      msg <- sprintf("interval 必须落在 |mat| 的取值范围内：[%.4g, %.4g]",
                     min_cell, max_cell)
      stop(msg)
    }
    thresholds <- seq(min(interval), max(interval), length.out = nr_thresholds)
  }

  ## 图像保存目录
  plot_files <- list()
  if (save_plots) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    .png <- function(fname, width=1000, height=1000, res=120) {
      png(filename = file.path(out_dir, fname), width = width, height = height, res = res)
    }
  }

  ## ------------------------- 内部小工具 -------------------------
  .rank_qr <- function(A, tol = NULL) {
    # 基于 QR 分解的数值秩（避免引入 Matrix 依赖）
    q <- qr(A, tol = tol)
    q$rank
  }

  .remove_outliers_IQR <- function(x, factor = 1.5) {
    q25 <- stats::quantile(x, 0.25, names = FALSE)
    q75 <- stats::quantile(x, 0.75, names = FALSE)
    iqr <- q75 - q25
    x[x >= (q25 - factor * iqr) & x <= (q75 + factor * iqr)]
  }

  .unfold_gauss <- function(ev, bw = "nrd0") {
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

  .unfold_spline <- function(ev, nfit = 51) {
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

  .ev_unfold <- function(A, method, bw, nfit, drop_outliers = TRUE) {
    ev <- eigen(A, only.values = TRUE)$values
    ev <- as.numeric(ev)
    if (drop_outliers) {
      ev0 <- ev
      ev  <- unique(.remove_outliers_IQR(ev))
      # 记录删掉多少个，仅供参考
      removed <- length(ev0) - length(ev)
    } else {
      removed <- 0
    }
    if (method == "gaussian") {
      u <- .unfold_gauss(ev, bw)
    } else {
      u <- .unfold_spline(ev, nfit)
    }
    u$nr.outliers.removed <- removed
    u
  }

  .wigner_pdf <- function(x) (pi/2) * x * exp(-pi * x^2 / 4)
  .exp_pdf    <- function(x) exp(-x)

  .kld <- function(obs_x, obs_pdf, exp_pdf) {
    # 需要概率密度非负且归一化
    obs_pdf <- pmax(obs_pdf, 0)
    exp_pdf <- pmax(exp_pdf, 0)
    s1 <- sum(obs_pdf); s2 <- sum(exp_pdf)
    if (s1 == 0 || s2 == 0) return(Inf)
    p <- obs_pdf / s1
    q <- exp_pdf / s2
    idx <- which(p > 0 & q > 0)
    sum(p[idx] * log(p[idx] / q[idx]))
  }

  .nnsd_metrics <- function(sp, max_cut = 3, nbins = 51) {
    sp <- sp[is.finite(sp) & sp > 0]
    if (!is.null(max_cut)) sp <- sp[sp <= max_cut]
    if (length(sp) < 10) return(list(
      N = length(sp), kld_exp = NA, kld_wig = NA, ks_p = NA, sse_exp = NA,
      frac_small = NA
    ))
    # 直方图估计
    brks <- seq(min(sp), max(sp), length.out = nbins)
    h <- hist(sp, breaks = brks, plot = FALSE)
    mids <- h$mids
    obs  <- h$density
    # 目标分布
    pdf_exp <- .exp_pdf(mids)
    pdf_wig <- .wigner_pdf(mids)
    # KLD
    kld_exp <- .kld(mids, obs, pdf_exp)
    kld_wig <- .kld(mids, obs, pdf_wig)
    # KS 与指数
    ks_p <- tryCatch({
      suppressWarnings(stats::ks.test(unique(sp), "pexp", 1)$p.value)
    }, error = function(e) NA_real_)
    # SSE to exponential（等面积分段）
    sse_exp <- {
      # 用 N 等面积区间比较 observed vs exp 的面积（粗略版）
      dens <- stats::density(sp, bw = "nrd0", n = 512)
      x <- seq(min(sp), max(sp), length.out = 1000)
      y_obs <- stats::approx(dens$x, dens$y, xout = x, rule = 2)$y
      y_exp <- .exp_pdf(x)
      A <- exp(-min(sp)) - exp(-max(sp))
      # 构造等面积 N=20 段
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
    # 小间距比例（监测“0堆积”）
    eps <- max_cut / 1000
    frac_small <- mean(sp < eps)
    list(
      N = length(sp),
      kld_exp = kld_exp, kld_wig = kld_wig,
      ks_p = ks_p, sse_exp = sse_exp,
      frac_small = frac_small
    )
  }

  .effective_mat <- function(A, tol = 0) {
    # 丢弃全零行列（保持对角）
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

  .denoise <- function(A, thr, keep_diag = TRUE) {
    if (keep_diag) dg <- diag(A)
    A[abs(A) < abs(thr)] <- 0
    if (keep_diag) diag(A) <- dg
    A
  }

  ## --------------------- 阈值扫描主循环（无交互） ---------------------
  res_tbl <- vector("list", length(thresholds))
  last_unfold <- NULL

  for (i in seq_along(thresholds)) {
    thr <- thresholds[i]
    if (verbose) message(sprintf("  [Scan] %3d/%d | threshold = %.4g", i, length(thresholds), thr))

    # 阈值去噪（保持对角）
    A <- .denoise(mat, thr, keep_diag = TRUE)

    # 丢零行列得到有效矩阵（可选）
    effA <- if (discard.zeros) .effective_mat(A) else A

    effN <- nrow(effA)
    if (verbose) message(sprintf("          Effective dimension: %d", effN))
    if (effN == 0 || effN < min.mat.dim) {
      if (verbose) message("          Too small. Stop scanning.")
      # 提前截断：把余下都设为 NA
      res_tbl <- res_tbl[seq_len(i)]
      thresholds <- thresholds[seq_len(i)]
      break
    }

    # 展开（unfold）
    u <- .ev_unfold(effA, method = unfold.method, bw = bandwidth,
                    nfit = nr.fit.points, drop_outliers = discard.outliers)
    sp <- u$ev.spacing
    # NNSD 指标
    met <- .nnsd_metrics(sp, max_cut = max.ev.spacing, nbins = 51)

    # 基本计数
    nr_zero   <- sum(A == 0)
    uniq_ev   <- NA_integer_
    max_mult  <- NA_integer_
    # 粗略估计（不保存特征向量，按四舍五入聚类）
    ev_vals <- eigen(effA, only.values = TRUE)$values
    ev_r <- round(ev_vals, 8)
    tbl <- table(ev_r)
    uniq_ev <- length(tbl)
    max_mult <- max(tbl)

    # 记录
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

    # 若需要保存图
    if (save_plots) {
      # NNSD 直方 + 理论曲线
      .png(sprintf("NNSD_%03d.png", i))
      h <- hist(sp[sp>0 & sp <= max.ev.spacing], breaks = 51, col = "darkolivegreen2",
                main = sprintf("NNSD (thr=%.4g)", thr), xlab = "spacing", freq = FALSE)
      curve(.exp_pdf, from = 0, to = max(h$breaks), add = TRUE, col = "red", lwd = 2)
      curve(.wigner_pdf, from = 0, to = max(h$breaks), add = TRUE, col = "blue", lwd = 2)
      legend("topright", c("Exp", "Wigner"), col = c("red","blue"), lty = 1, lwd = 2, bty = "n")
      dev.off()

      # KS pvalue, SSE, KLD 这些会在结束后画整合图
    }

    last_unfold <- u
  }

  score_df <- do.call(rbind, res_tbl)
  rownames(score_df) <- NULL

  if (nrow(score_df) == 0) stop("没有有效的阈值结果（可能 min.mat.dim 过大或矩阵过稀）。")

  ## --------------------- 自动选阈值（无交互） ---------------------
  # 目标：希望间距分布接近 Exp（Poisson），远离 Wigner；同时矩阵不能太小，spacing 数量充分，
  #      小间距堆积(frac_small)不要过大；SSE 越小越好。
  # 做一个简单综合评分（全部转为“越大越好”再加权）：
  #  score = w1*scale(ks_p) + w2*scale(-kld_exp) + w3*scale(-sse_exp) + w4*scale(-frac_small) + w5*scale(eff_dim)
  #  若某列有 NA，则不参与该部分评分。
  scl <- function(x) if (all(is.na(x))) rep(0, length(x)) else as.numeric(scale(x))
  w1 <- 0.35; w2 <- 0.25; w3 <- 0.20; w4 <- 0.10; w5 <- 0.10
  score <- w1 * scl(score_df$ks_p) +
    w2 * scl(-score_df$kld_exp) +
    w3 * scl(-score_df$sse_exp) +
    w4 * scl(-score_df$frac_small) +
    w5 * scl(score_df$eff_dim)

  score_df$auto_score <- score

  # 备选：交叉点启发（KLD_wig > KLD_exp 或 -logLik 到 Exp < 到 Wigner）
  # 这里用 KLD：挑选 KLD_exp 最小、且 eff_dim >= 中位数 的候选
  cand1_idx <- which(score_df$kld_exp == min(score_df$kld_exp, na.rm = TRUE))
  cand1_idx <- cand1_idx[which.max(score_df$eff_dim[cand1_idx])]

  # 终选：优先综合分最高；如和 cand1 距离很近（|thr差|<1e-8），用 cand1
  best_idx <- which.max(score_df$auto_score)
  chosen_idx <- if (abs(score_df$threshold[best_idx] - score_df$threshold[cand1_idx]) < 1e-8) {
    cand1_idx
  } else {
    best_idx
  }
  chosen_threshold <- score_df$threshold[chosen_idx]
  chosen_reason <- sprintf(
    "综合分最高（权重: ks_p %.2f, KLD_exp %.2f, SSE %.2f, 小间距比例 %.2f, 有效维度 %.2f）",
    w1, w2, w3, w4, w5
  )

  ## --------------------- 汇总图（非交互） ---------------------
  if (save_plots) {
    # p-value vs threshold
    .png("KS_p_vs_threshold.png", 1200, 900, 140)
    plot(score_df$threshold, score_df$ks_p, type = "b", col = "red",
         main = "KS p-value vs threshold", xlab = "threshold", ylab = "p-value")
    abline(v = chosen_threshold, col = "darkgreen", lty = 2)
    dev.off()

    # SSE vs threshold
    .png("SSE_vs_threshold.png", 1200, 900, 140)
    plot(score_df$threshold, score_df$sse_exp, type = "b", col = "red",
         main = "SSE (NNSD vs Exp) vs threshold", xlab = "threshold", ylab = "SSE")
    abline(v = chosen_threshold, col = "darkgreen", lty = 2)
    dev.off()

    # KLD vs threshold
    .png("KLD_vs_threshold.png", 1200, 900, 140)
    matplot(score_df$threshold, cbind(score_df$kld_wig, score_df$kld_exp),
            type = "b", pch = c(1,2), col = c("blue","red"),
            main = "KLD to Wigner / Exp vs threshold", xlab = "threshold", ylab = "KLD")
    legend("topright", c("KLD(Wigner)", "KLD(Exp)"), col = c("blue","red"), pch = c(1,2), bty = "n")
    abline(v = chosen_threshold, col = "darkgreen", lty = 2)
    dev.off()

    # auto score vs threshold
    .png("AutoScore_vs_threshold.png", 1200, 900, 140)
    plot(score_df$threshold, score_df$auto_score, type = "b", col = "purple",
         main = "Auto score vs threshold", xlab = "threshold", ylab = "score (higher is better)")
    abline(v = chosen_threshold, col = "darkgreen", lty = 2)
    dev.off()

    plot_files <- list(
      nnsd_each = list.files(out_dir, pattern = "^NNSD_\\d+\\.png$", full.names = TRUE),
      ks = file.path(out_dir, "KS_p_vs_threshold.png"),
      sse = file.path(out_dir, "SSE_vs_threshold.png"),
      kld = file.path(out_dir, "KLD_vs_threshold.png"),
      score = file.path(out_dir, "AutoScore_vs_threshold.png")
    )
  }

  ## --------------------- 返回 ---------------------
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
        interval = interval %||% c(min_cell, max_cell),
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

`%||%` <- function(a, b) if (is.null(a)) b else a
