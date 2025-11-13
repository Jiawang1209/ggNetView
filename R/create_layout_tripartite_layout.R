create_layout_tripartite_layout <- function(
    graph_obj,
    scale = T,
    r = 1,
    anchor_dist = 10,
    node_add = 7,
    orientation = c("up","down","left","right"),
    angle = 0
){
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 基础半径与等腰直角三角形三个锚点（先按“up”放置，再统一旋转）
  radius <- r * 6

  anchors <- list(
    c(-anchor_dist, 0),  # 左
    c( anchor_dist, 0),  # 右
    c( 0,  anchor_dist)  # 上（此处是直角顶点，三角形为等腰直角）
  )

  # 取节点与模块
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_levels <- node_df$Modularity %>% droplevels() %>% levels() %>% as.character()
  module_list <- node_df %>% dplyr::group_split(Modularity)

  n_vec <- purrr::map_int(module_list, nrow)

  if (length(n_vec) < 3) {
    stop("Tripartite layout requires at least 2 modules in `Modularity`")
  }
  if (length(n_vec) > 3) {
    message("`tripartite layout more than 2 modules detected.")
  }

  # 计算每个模块同心圆分桶
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

  # 生成“以某锚点为中心”的同心圆布局（带交错 offset）
  concentric_from_anchor <- function(cx, cy, info_df, r_step){
    ly <- data.frame(x = cx, y = cy)  # 第一圈放 1 个
    offset <- 0
    prev_n <- info_df$number_node
    if (nrow(info_df) >= 2) {
      for (index in 2:nrow(info_df)) {
        if (index == 2) {
          l <- 2 * pi * (0:(prev_n[index]-1)) / prev_n[index]
        } else {
          offset <- (offset + pi/prev_n[index]) %% (2*pi)
          l <- offset + 2 * pi * (0:(prev_n[index]-1)) / prev_n[index]
        }
        x <- cx + sin(l) * (index - 1) * r_step
        y <- cy + cos(l) * (index - 1) * r_step
        ly <- dplyr::bind_rows(ly, data.frame(x = x, y = y))
      }
    }
    ly
  }

  # 为三个模块依次生成
  ly_list <- list()
  for (i in 1:3) {
    cx <- anchors[[i]][1]; cy <- anchors[[i]][2]
    ly_i <- concentric_from_anchor(cx, cy, n_vec_node[[i]], r_step = r)
    ly_i$group <- mod_levels[i]  # 可用于后续按模块对齐/合并
    ly_list[[i]] <- ly_i
  }

  ly <- dplyr::bind_rows(ly_list)

  # ---- 统一旋转（绕原点）----
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(ly[, c("x","y")])
    ly[, c("x","y")] <- t(Rm %*% t(xy))
  }

  rownames(ly) <- NULL
  ly
}
