create_layout_rings <- function(graph_obj,
                                ring_n = NULL,
                                r = 0.25,
                                orientation = c("up","down","left","right"),
                                angle = 0){
  graph_obj = graph_obj %>%
    tidygraph::arrange(group)

  # 旋转角度
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # 获取节点
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  # 对每一个分组，进行统计个数和排序
  node_df_stat <- node_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n = n())

  # 获取边
  graph_obj %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble()

  if (is.null(ring_n)) {
    # 自动给n区分出来
    ring_n = dim(node_df_stat)[1]
  }

  circle_list <- list()

  for (i in 1:ring_n) {
    # 设置外圈的点
    n_points <- node_df_stat$n[i]
    radius <- r
    center_x <- 0
    center_y <- 0

    # 计算每一个点的角度
    angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points+1)]

    # 计算坐标
    x <- center_x + radius * cos(angles)
    y <- center_y + radius * sin(angles)

    # 创建数据框

    circle_df <- data.frame(
      id = 1:n_points,
      x = x,
      y = y
    )

    r = r + 0.5
    circle_list[[i]] <- circle_df

  }

  circle_df_layout <- do.call(rbind, circle_list)

  # 自定一个布局
  layout_manual <- ggraph::create_layout(graph_obj, layout = "circle")

  layout_manual_2 <- layout_manual %>%
    dplyr::mutate(x = circle_df_layout$x,
                  y = circle_df_layout$y)

  ly <- ggraph::create_layout(graph_obj,
                      layout = "manual",
                      x = layout_manual_2$x,
                      y = layout_manual_2$y)

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
