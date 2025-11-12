create_layout_bipartite_layout <- function(
    graph_obj,
    r = 6,
    node_add = NULL,
    orientation = c("up","down","left","right"),
    angle = 0
  ){

  # 既然，要画二分网络，那么你肯定要有二分网络

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

  Modularity_name <- node_df$Modularity  %>% droplevels() %>% levels() %>% as.character()


  # 那么肯定是左边一个 右边一个
  module_list <- purrr::map(graph_obj %>%
                              tidygraph::activate(nodes) %>%
                              tidygraph::as_tibble() %>%
                              dplyr::group_split(Modularity),
                            ~.x)

  module_list

  module_number <- purrr::map(module_list, ~dim(.x)[1])

  names(module_number) <- Modularity_name

  ly <- purrr::map2(module_number, seq_along(module_number), function(n_points, i) {
    # 中心点
    if (i == 1) {
      center_x <- 0 - radius*1.5
    }else{
      center_x <- 0 + radius*1.5
    }

    center_y <- 0

    angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points + 1)]
    data.frame(
      x = center_x + radius * cos(angles),
      y = center_y + radius * sin(angles)
    )
  }) %>%
    do.call(rbind, .)

  # ly %>%
  #   dplyr::mutate(number = 1:n())  %>%
  #   ggplot() +
  #   geom_point(aes(x = x, y = y, color = Group), size = 5) +
  #   geom_text(aes(x = x, y = y, color = Group, label = number), color ="#000000", size = 5)

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
