create_layout_fr1 <- function(
    graph_obj,
    node_add = NULL,
    r=NULL,
    scale = T,
    orientation = c("up","down","left","right"),
    angle = 0
){

  # 旋转角度
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

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

  # layout
  ly <- igraph::layout_with_fr(igraph_obj, niter = 1000) %>%
    as.data.frame() %>%
    purrr::set_names(c("x", "y"))


  # ly <- as.data.frame(ly) %>% dplyr::select(1,2)


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
