create_layout_quadripartite_gephi_layout <- function(
    graph_obj,
    r = 1,
    node_add = 7,
    scale = T,
    orientation = c("up","down","left","right"),
    angle = 0
){
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # ---- 正方形的四个锚点 ----
  radius <- r * 6
  anchors <- list(
    c(-radius, -radius),  # 左下
    c(-radius,  radius),  # 左上
    c( radius,  radius),  # 右上
    c( radius, -radius)   # 右下
  )

  # 平移到质心在原点
  A <- do.call(rbind, anchors)
  A_centered <- sweep(A, 2, colMeans(A), "-")
  anchors <- split(A_centered, row(A_centered))

  # ---- 获取节点和模块 ----
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_levels <- node_df$Modularity %>% droplevels() %>% levels() %>% as.character()
  module_list <- node_df %>% dplyr::group_split(Modularity)
  n_vec <- purrr::map_int(module_list, nrow)

  if (length(n_vec) < 4) {
    stop("Quadripartite layout 需要至少 4 个模块（来自列 Modularity）。")
  }
  if (length(n_vec) > 4) {
    message("检测到超过 4 个模块，仅使用前 4 个模块进行正方形布局。")
    module_list <- module_list[1:4]
    n_vec <- n_vec[1:4]
    mod_levels <- mod_levels[1:4]
  }

  # ---- 同心圆节点分层 ----
  circle_layout <- function(n, node_add){
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
  }

  n_vec_node <- purrr::map(n_vec, ~{
    data.frame(
      number_circle = seq_along(circle_layout(.x, node_add)),
      number_node   = circle_layout(.x, node_add)
    )
  })

  # ---- 同心圆布局函数 ----
  concentric_from_anchor <- function(cx, cy, info_df, r_step){
    ly <- data.frame(x = cx, y = cy)
    offset <- 0
    prev_n <- info_df$number_node
    if (nrow(info_df) >= 2) {
      for (index in 2:nrow(info_df)) {
        if (index == 2) {
          l <- 2*pi * (0:(prev_n[index]-1)) / prev_n[index]
        } else {
          offset <- (offset + pi/prev_n[index]) %% (2*pi)
          l <- offset + 2*pi * (0:(prev_n[index]-1)) / prev_n[index]
        }
        x <- cx + sin(l) * (index - 1) * r_step
        y <- cy + cos(l) * (index - 1) * r_step
        ly <- dplyr::bind_rows(ly, data.frame(x = x, y = y))
      }
    }
    ly
  }

  # ---- 按模块生成布局 ----
  ly_list <- list()
  for (i in 1:4) {
    cx <- anchors[[i]][1]; cy <- anchors[[i]][2]
    ly_i <- concentric_from_anchor(cx, cy, n_vec_node[[i]], r_step = r)
    ly_i$group <- mod_levels[i]
    ly_list[[i]] <- ly_i
  }

  ly <- dplyr::bind_rows(ly_list)

  # ---- 统一旋转 ----
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(ly[, c("x","y")])
    ly[, c("x","y")] <- t(Rm %*% t(xy))
  }

  rownames(ly) <- NULL
  ly
}
