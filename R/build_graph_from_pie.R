#' Build a pie graph object from a data frame
#'
#' @param df Data frame.
#'   Edge list with columns \code{from}, \code{to}, and optionally \code{weight}.
#'   If \code{weight} is absent, an unweighted graph is constructed.
#' @param node_annotation Data Frame
#' The annotation file of nodes in network
#' @param directed  Logical (default: \code{FALSE}).
#'   Whether edges between nodes are directed.
#' @param seed Integer, optional.
#'   Random seed for reproducibility; if \code{NULL}, no seed is set.
#'
#' @returns An graph object representing the correlation network.
#'   Node/edge attributes include correlation statistics and (optionally) module labels.
#' @export
#'
#' @examples NULL
build_graph_from_pie <- function(df,
                                 node_annotation = NULL,
                                 directed = F,
                                 seed = 1115){
  set.seed(seed)

  # 构建igraph对象
  g <- igraph::graph_from_data_frame(
    d = df,
    vertices = node_annotation,
    directed = directed
  )

  # 删除自相关
  g <- igraph::simplify(g)

  # 删除孤立节点
  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))

  # 构建ggraph对象
  graph_obj <- tidygraph::as_tbl_graph(g)

  return(graph_obj)

}
