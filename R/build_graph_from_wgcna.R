#' Build a graph object from WGCNA result
#'
#' @param wgcna_tom WGCNA TOM matrix.
#'
#' @param module Module data frame
#' @param directed Logical (default: \code{FALSE}).
#'   Whether edges between nodes are directed.
#' @param seed Integer, optional.
#'   Random seed for reproducibility; if \code{NULL}, no seed is set.
#'
#' @returns An graph object representing the correlation network.
#' Node/edge attributes from WGCNA TOM matrix (optionally) module labels.
#'
#' @export
#'
#' @examples NULL
build_graph_from_wgcna <- function(wgcna_tom,
                                   module = NULL,
                                   directed = F,
                                   seed = 1115){

  set.seed(seed)

  # 构建igraph对象
  g <- igraph::graph_from_data_frame(
    d = wgcna_tom,
    vertices = module,
    directed = directed
  )

  graph_obj <- tidygraph::as_tbl_graph(wgcna_tom) %>%
    tidygraph::left_join(module, by = c("name" = "ID")) %>%
    tidygraph::mutate(modularity2 = factor(Module),
                      modularity2 = factor(modularity2),
                      modularity3 = as.character(modularity2),
                      degree = tidygraph::centrality_degree(mode = "out"),
                      strength = tidygraph::centrality_degree(weights = weight)) %>%
    tidygraph::arrange(modularity2, desc(degree))

  return(graph_obj)
}
