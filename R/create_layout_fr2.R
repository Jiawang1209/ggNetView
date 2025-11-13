create_layout_fr2 <- function(
    graph_obj,
    node_add = NULL,
    r=NULL,
    scale = T,
    anchor_dist = 10,
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
  e <- igraph::as_edgelist(igraph_obj, names=FALSE)

  ly <- qgraph::qgraph.layout.fruchtermanreingold(e,
                                                 vcount = igraph::vcount(igraph_obj)) %>%
    as.data.frame() %>%
    purrr::set_names(c("x", "y"))

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
