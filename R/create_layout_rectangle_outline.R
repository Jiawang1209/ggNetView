# 只生成“最外圈”的等腰直角三角形轮廓布局
# r: 半边长；直角顶点在 (-r, -r)，另两顶点在 ( r, -r) 与 (-r,  r)
create_layout_rectangle_outline <- function(
    graph_obj,
    r = 6,
    node_add = NULL,
    orientation = c("up","down","left","right"),
    angle = 0  # 弧度；在 orientation 基础上的微调
){
  # 方向与旋转
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  if (n == 0) {
    return(tibble::tibble(name = character(), x = numeric(), y = numeric()))
  }

  a <- r  # 半边长

  # 三个顶点（未旋转前）：直角在左下
  V0 <- c(-a, -a)  # 左下（直角）
  V1 <- c( a, -a)  # 右下
  V2 <- c(-a,  a)  # 左上
  # 三条边：E01（底边，长=2a）、E02（左边，长=2a）、E12（斜边，长=2a*sqrt(2)）

  # 在边上取“内点”（不含端点）
  edge_points <- function(p1, p2, k){
    if (k <= 0) return(matrix(numeric(0), ncol = 2))
    t <- seq_len(k) / (k + 1)  # 避开端点
    cbind(p1[1] + (p2[1]-p1[1]) * t,
          p1[2] + (p2[2]-p1[2]) * t)
  }

  if (n == 1) {
    coords <- matrix(V0, ncol = 2, byrow = TRUE)
  } else if (n == 2) {
    coords <- rbind(V0, V1)
  } else {
    # 先放三个顶点
    m <- n - 3
    # 按边长比例分配剩余 m 个点：E01:E02:E12 = 1 : 1 : sqrt(2)
    w <- c(1, 1, sqrt(2))
    p <- w / sum(w)
    k_base <- floor(m * p)
    rem <- m - sum(k_base)
    # 按小数部分从大到小分配余数，保持比例更精确
    frac <- m * p - k_base
    if (rem > 0) {
      order_idx <- order(frac, decreasing = TRUE)
      k_base[order_idx[seq_len(rem)]] <- k_base[order_idx[seq_len(rem)]] + 1
    }
    k01 <- k_base[1]  # 底边内点数
    k02 <- k_base[2]  # 左边内点数
    k12 <- k_base[3]  # 斜边内点数

    X <- rbind(
      V0,
      edge_points(V0, V1, k01), V1,
      edge_points(V1, V2, k12), V2,
      edge_points(V2, V0, k02)
    )
    coords <- X
  }

  coords <- as.data.frame(coords)
  names(coords) <- c("x","y")

  # 整体旋转（绕原点）
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(coords[, c("x","y")])
    coords[, c("x","y")] <- t(Rm %*% t(xy))
  }

  # 返回，带 name 便于 join
  ly <- data.frame(
    x = coords$x,
    y = coords$y
  )
}
