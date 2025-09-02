#' Create gephi layout
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#'   The network object to be visualized.
#' @param node_add Integer (default = 7).
#'   Number of nodes to add in each layer of the layout.
#' @param r Numeric (default = 1).
#'   Radius increment for concentric or layered layouts.
#' @param orientation Character string.
#'   custom orientation; one of
#'   "up","down","left","right"
#' @param angle Integer  (default = 0).
#'  change  orientation angle
#'
#' @returns layout.
#' @description internal function
#' @keywords internal
#'
#' @examples NULL
create_layout_gephi <- function(
    graph_obj,
    node_add = 7,
    r = 0.1,
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

  # 节点个数
  n <- dim(node_df)[1]

  # 现在开始正式写布局
  circle_layout <- function(n, node_add){
    # 初始化
    counts <- 1
    total <- 1
    i <- 2
    while (total < n) {
      add <- node_add * (i-1)
      if (total + add <= n) {
        counts <- c(counts, add)
        total <- total + add
      }else{
        # 最后一圈
        counts <- c(counts, n-total)
        total <- n
      }
      i <- i + 1
    }
    return(counts)
  }

  # 获取到layout的具体信息 每一层，有多少个点
  layout_df_info <- data.frame(
    number_circle = seq_along(circle_layout(n = n, node_add = node_add)),
    number = circle_layout(n = n, node_add = node_add)
  )
  layout_df_info

  # 然后我们开始写每一圈的真实布局
  # 第一圈没得说
  ly <- data.frame(x = 0, y = 0)
  offset <- 0
  prev_n <- layout_df_info$number

  # 第二圈 到 最后一圈
  for (index in 2:(dim(layout_df_info)[1])) {
    if (index == 2) {
      # index = 2
      # l <- seq(0, 2*pi, length.out = prev_n[index])
      l <- 2* pi * (0:(prev_n[index]-1)) / prev_n[index]
    }else{
      # 第三圈开始 错开半个身位置
      offset <- pi/prev_n[index] %% (2*pi) + offset
      # l <- offset + seq(0, 2*pi, length.out = prev_n[index])
      l <- offset + (2* pi * (0:(prev_n[index]-1)) / prev_n[index])
    }

    x <- sin(l) * (index-1)*r
    y <- cos(l) * (index-1)*r
    ly_tmp <- data.frame(x = x,
                         y = y)

    ly <- dplyr::bind_rows(ly,ly_tmp)

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
