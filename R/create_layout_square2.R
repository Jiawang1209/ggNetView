create_layout_square2 <- function(
    graph_obj,
    node_add = 8,
    r = 0.1,
    scale = T,
    orientation = c("up","down","left","right"),
    angle = 0 # 在 orientation 基础上的微调（弧度）
){
  # 旋转角度
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 获取节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)

  # 每圈点数分配
  ring_counts <- (function(n, node_add){
    counts <- 1
    total  <- 1
    i <- 2
    while (total < n) {
      add <- node_add * (i - 1)
      if (total + add <= n) {
        counts <- c(counts, add)
        total  <- total + add
      } else {
        counts <- c(counts, n - total)
        total  <- n
      }
      i <- i + 1
    }
    counts
  })(n, node_add)

  layout_df_info <- data.frame(
    number_circle = seq_along(ring_counts),
    number = ring_counts
  )

  # 把 u ∈ [0,1) 映射到正方形边
  square_param_to_xy <- function(u, a){
    u <- u %% 1
    s <- 0.25
    x <- numeric(length(u))
    y <- numeric(length(u))

    # 顶边 (-a,a)→(a,a)
    idx0 <- u < s
    t0 <- u[idx0] / s
    x[idx0] <- -a + 2*a*t0
    y[idx0] <-  a

    # 右边 (a,a)→(a,-a)
    idx1 <- (u >= s) & (u < 2*s)
    t1 <- (u[idx1] - s) / s
    x[idx1] <-  a
    y[idx1] <-  a - 2*a*t1

    # 底边 (a,-a)→(-a,-a)
    idx2 <- (u >= 2*s) & (u < 3*s)
    t2 <- (u[idx2] - 2*s) / s
    x[idx2] <-  a - 2*a*t2
    y[idx2] <- -a

    # 左边 (-a,-a)→(-a,a)
    idx3 <- u >= 3*s
    t3 <- (u[idx3] - 3*s) / s
    x[idx3] <- -a
    y[idx3] <- -a + 2*a*t3

    data.frame(x = x, y = y)
  }

  # 生成布局
  ly <- data.frame(x = 0, y = 0)  # 中心点

  for (index in 2:nrow(layout_df_info)) {
    m <- layout_df_info$number[index]
    half_len <- (index - 1) * r
    u <- (0:(m-1)) / m   # 没有任何偏移，整齐对齐
    coords <- square_param_to_xy(u, a = half_len)
    ly <- dplyr::bind_rows(ly, coords)
  }

  # 开始旋转
  # 统一旋转（绕原点）
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(ly[, c("x","y")])
    ly[, c("x","y")] <- t(Rm %*% t(xy))
  }

  return(ly)
}
