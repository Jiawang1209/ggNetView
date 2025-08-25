#' Create petal layout
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#'   The network object to be visualized.
#' @param node_add Integer (default = 7).
#'   Number of nodes to add in each layer of the layout.
#' @param r Numeric (default = 1).
#'   Radius increment for concentric or layered layouts.
#' @param petals Numeric (default = 6).
#' @param amp Numeric (default = 0.35).
#'
#' @returns Layout
#' @description internal function
#' @keywords internal
#'
#' @examples NULL
create_layout_petal <- function(graph_obj, node_add = 7, r = 0.1,
                                petals = 6, amp = 0.35){
  # 取节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)

  # 每圈点数（延续你之前的配额逻辑）
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

  # 生成布局
  ly <- data.frame(x = 0, y = 0)  # 第1圈中心点

  # 从第二圈开始：按“花瓣曲线”等角度取点（无偏移，保证对称）
  for (index in 2:nrow(layout_df_info)) {
    m <- layout_df_info$number[index]
    R_base <- (index - 1) * r            # 这一圈的基准半径

    theta <- 2 * pi * (0:(m - 1)) / m    # 等角度采样
    radius <- R_base * (1 + amp * cos(petals * theta))
    # 为避免amp过大产生半径过小甚至负数，建议 amp ∈ [0, 0.9]

    x <- radius * cos(theta)
    y <- radius * sin(theta)

    ly <- dplyr::bind_rows(ly, data.frame(x = x, y = y))
  }

  return(ly)
}
