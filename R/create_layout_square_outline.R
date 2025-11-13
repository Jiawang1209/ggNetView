create_layout_square_outline <- function(
    graph_obj,
    r = 6,
    node_add = NULL,
    orientation = c("up","down","left","right"),
    anchor_dist = 10,
    angle = 0  # 弧度
){
  # 方向与旋转角
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 节点
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  if (n == 0) return(tibble::tibble(name = character(), x = numeric(), y = numeric()))

  # 把 r(顶点距)换算为半边长 a（确保几何是“正方向”的正方形）
  a <- r / sqrt(2)

  # 四个角（未旋转前）：按 TL, TR, BR, BL 排序
  TL <- c(-a,  a); TR <- c( a,  a)
  BR <- c( a, -a); BL <- c(-a, -a)

  # 边内取点（不含端点）
  edge_points <- function(p1, p2, k){
    if (k <= 0) return(matrix(numeric(0), ncol = 2))
    t <- seq_len(k) / (k + 1)  # 避开角
    cbind(p1[1] + (p2[1]-p1[1]) * t,
          p1[2] + (p2[2]-p1[2]) * t)
  }

  # 小样本优先占角（越少越把“方”的感觉保住）
  if (n <= 4) {
    pts <- switch(n,
                  `1` = rbind(TR),
                  `2` = rbind(TL, TR),
                  `3` = rbind(TL, TR, BR),
                  `4` = rbind(TL, TR, BR, BL)
    )
    coords <- as.data.frame(pts); names(coords) <- c("x","y")
  } else {
    # 先占四角，再分配剩余 m 个点到四条等长边
    m <- n - 4
    base <- m %/% 4
    rem  <- m %% 4
    # 尽量对称：上/下优先，再右/左（上、右、下、左）
    counts <- rep(base, 4)
    if (rem >= 1) counts[1] <- counts[1] + 1  # 上
    if (rem >= 2) counts[3] <- counts[3] + 1  # 下
    if (rem >= 3) counts[2] <- counts[2] + 1  # 右

    X <- rbind(
      TL,
      edge_points(TL, TR, counts[1]), TR,
      edge_points(TR, BR, counts[2]), BR,
      edge_points(BR, BL, counts[3]), BL,
      edge_points(BL, TL, counts[4])
    )
    coords <- as.data.frame(X); names(coords) <- c("x","y")
  }

  # 旋转
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(coords[, c("x","y")])
    coords[, c("x","y")] <- t(Rm %*% t(xy))
  }

  tibble::tibble(
    name = if ("name" %in% names(node_df)) node_df$name else as.character(seq_len(nrow(coords))),
    x = coords$x,
    y = coords$y
  )
}
