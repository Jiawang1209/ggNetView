create_layout_heart_centered <- function(
    graph_obj,
    r = 0.15,                    # 层间尺度步长（整体大小递增）
    node_add = 8L,          # 每层基础分段数（默认更紧凑：8, 16, 24, ...）
    orientation = c("up","down","left","right"),
    scale = T,
    anchor_dist = 10,
    angle = 0,                   # 在 orientation 基础上的微调（弧度）
    y_squash = 1.0               # 纵向压缩（<1更“薄”，>1更“胖”）
){
  # 心形同心布局（无重复，质心居中；首点=整体中心(0,0)）
  # 第1层：1个点(0,0)；第k层(k>=2)：base_per_ring*(k-1)个点

  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 1) 节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  if (n <= 0) return(data.frame(x=numeric(), y=numeric(), ring=integer(), idx=integer()))

  # 2) 层配额：1, base*(1), base*2, ... 截断到 n
  ring_sizes <- c(1L)
  total <- 1L; k <- 2L
  while (total < n) {
    add  <- node_add * (k - 1L)
    take <- min(add, n - total)
    ring_sizes <- c(ring_sizes, take)
    total <- total + take
    k <- k + 1L
  }

  # 3) 经典心形参数方程（未缩放、未压缩）
  heart_raw <- function(t){
    x <- 16 * (sin(t))^3
    y <- 13 * cos(t) - 5 * cos(2*t) - 2 * cos(3*t) - cos(4*t)
    cbind(x, y)
  }

  # 4) 给定尺度 L，生成“质心在原点”的心形轮廓（用于采样）
  #    - 先致密采样 -> 算质心 -> 平移到原点
  #    - y 方向可用 y_squash 调薄/调胖
  heart_centered_dense <- function(L, squash = 1.0, dense = 1024L){
    tt <- seq(0, 2*pi, length.out = dense + 1L)[- (dense + 1L)]  # 不含 2π 端点
    xy <- heart_raw(tt)
    xy[,2] <- xy[,2] * squash
    xy <- xy * L
    cx <- mean(xy[,1]); cy <- mean(xy[,2])
    xy <- cbind(xy[,1] - cx, xy[,2] - cy)
    xy
  }

  # 5) 按弧长“等距”的中点采样（避免端点/重复）
  #    返回该层 m 个点（m <= M_full）
  sample_on_heart <- function(L, m, M_full, squash = 1.0){
    if (m <= 0) return(matrix(numeric(0), ncol = 2))
    if (L == 0)  return(matrix(c(0,0), ncol = 2))

    poly <- heart_centered_dense(L, squash = squash, dense = 2048L)
    # 计算累积弧长
    seg  <- sqrt(rowSums((poly[c(2:nrow(poly),1), ] - poly)^2))
    s    <- c(0, cumsum(seg))
    P    <- s[length(s)]
    # 将周长分为 M_full 段，取每段中点位置
    s_mid_all <- ((0:(M_full - 1L)) + 0.5) / M_full * P
    s_target  <- s_mid_all[seq_len(m)]

    # 线性插值：根据 s_target 找到所在边并插值坐标
    idx <- findInterval(s_target, s, all.inside = TRUE)
    # s[idx] 到 s[idx+1] 之间
    s0  <- s[idx]; s1 <- s[idx + 1L]
    p0  <- poly[idx, , drop = FALSE]
    p1  <- poly[(idx %% nrow(poly)) + 1L, , drop = FALSE]
    w   <- (s_target - s0) / pmax(s1 - s0, .Machine$double.eps)
    xy  <- p0 + (p1 - p0) * w
    xy
  }

  out <- list()
  # 第1层：中心点
  out[[1]] <- data.frame(x = 0, y = 0, ring = 1L, idx = 1L)

  # 6) 逐层：等弧长中点采样，保证不重叠
  if (length(ring_sizes) >= 2) {
    for (ring in 2:length(ring_sizes)) {
      m <- ring_sizes[ring]
      L <- (ring - 1) * r
      if (m <= 0) next
      M_full <- node_add * (ring - 1L)

      xy <- sample_on_heart(L, m, M_full, squash = y_squash)

      # 统一旋转
      if (theta_shift != 0) {
        Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                       sin(theta_shift),  cos(theta_shift)), nrow = 2)
        xy <- t(Rm %*% t(xy))
      }

      out[[ring]] <- data.frame(x = xy[,1], y = xy[,2],
                                ring = ring, idx = seq_len(m))
    }
  }

  layout_df <- dplyr::bind_rows(out) %>%
    dplyr::select(x, y)

  # （可选）去重校验
  # eps <- 1e-12
  # key <- paste(round(layout_df$x/eps), round(layout_df$y/eps))
  # stopifnot(!any(duplicated(key)))

  layout_df
}
