#' Build a pie graph object from a data frame
#'
#' @param df Data frame.
#' Edge list with columns \code{from}, \code{to}, and optionally \code{weight}.
#' If \code{weight} is absent, an unweighted graph is constructed.
#' @param node_annotation Data Frame
#' The annotation file of nodes in network
#' @param directed  Logical (default: \code{FALSE}).
#'   Whether edges between nodes are directed.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns An graph object representing the correlation network.
#'   Node/edge attributes include correlation statistics and (optionally) module labels.
#' @export
#'
#' @examples
#' data(ppi_example)
#' obj <- build_graph_from_pie(
#'   df              = ppi_example$ppi,
#'   node_annotation = ppi_example$annotation
#' )
#' obj
build_graph_from_pie <- function(df,
                                 node_annotation = NULL,
                                 directed = FALSE,
                                 seed = 1115){
  .ggnv_local_seed(seed)


  g <- igraph::graph_from_data_frame(
    d = df,
    vertices = node_annotation,
    directed = directed
  )


  g <- igraph::simplify(g)


  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))


  graph_obj <- tidygraph::as_tbl_graph(g)

  return(graph_obj)

}
