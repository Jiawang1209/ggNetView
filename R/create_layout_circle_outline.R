create_layout_circle_outline <- function(
    graph_obj,
    node_add = NULL,
    r = 6,
    orientation = c("up","down","left","right"),
    angle = 0
){

  # 旋转角度
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # set radius
  radius = r

  # 获取节点
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  # 获取边
  graph_obj %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble()

  # as igraph object
  igraph_obj <- tidygraph::as.igraph(graph_obj)


  n_points <- node_df %>% dplyr::pull(name) %>% length()
  #
  # 计算每一个点的角度
  angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points+1)]
  center_x <- 0
  center_y <- 0
  #
  # # 计算坐标
  x <- center_x + radius * cos(angles)
  y <- center_y + radius * sin(angles)

  ly <- data.frame(
    x = x,
    y = y
  )


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
