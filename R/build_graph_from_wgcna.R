#' Build a graph object from WGCNA result
#'
#' @param wgcna_tom WGCNA TOM matrix.
#'
#' @param module Module data frame
#' @param directed Logical (default: \code{FALSE}).
#' Whether edges between nodes are directed.
#' @param seed nteger (default = 1115).
#' Random seed for reproducibility.
#' @param node_annotation Data Frame
#' The annotation file of nodes in network
#'
#' @returns An graph object representing the correlation network.
#' Node/edge attributes from WGCNA TOM matrix (optionally) module labels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # `wgcna_tom` is a long-format TOM edge list (from, to, weight)
#' # produced by `trans_TOM_in_WGCNA()`, and `module` is the WGCNA
#' # module-assignment data frame with columns (ID, Module).
#' obj <- build_graph_from_wgcna(
#'   wgcna_tom = wgcna_tom,
#'   module    = module
#' )
#' }
build_graph_from_wgcna <- function(wgcna_tom,
                                   module = NULL,
                                   node_annotation = NULL,
                                   directed = FALSE,
                                   seed = 1115){

  .ggnv_local_seed(seed)

  # Build the igraph with `module` as the authoritative vertex list.
  #
  # Subtle but important: downstream sparsification (e.g. a `threshold`
  # or `top_k` pass inside `trans_TOM_in_WGCNA()`) can remove every edge
  # incident to some gene.  If we inferred the vertex set from the edge
  # list alone (as the previous implementation effectively did via
  # `as_tbl_graph(wgcna_tom)`), those edgeless genes would silently
  # disappear from the plot even though `module` explicitly lists them.
  #
  # Passing only the ID column of `module` as `vertices` pins the full
  # vertex set without also smuggling in `Module`/other columns as
  # vertex attributes -- that lets the standard `left_join(module, ...)`
  # below attach `Module` cleanly, without triggering name collisions
  # (`Module.x` / `Module.y`) that would break `factor(Module)`.
  g <- igraph::graph_from_data_frame(
    d        = wgcna_tom,
    vertices = if (is.null(module)) NULL else module[, "ID", drop = FALSE],
    directed = directed
  )

  if (is.null(node_annotation)) {
    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::left_join(module, by = c("name" = "ID")) %>%
      tidygraph::mutate(modularity = factor(Module),
                        modularity2 = factor(modularity),
                        modularity3 = as.character(modularity2),
                        Modularity = modularity2,
                        Degree = tidygraph::centrality_degree(mode = "out"),
                        Strength = tidygraph::centrality_degree(weights = weight)
      ) %>%
      tidygraph::arrange(Modularity, desc(Degree))
  }else{

    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::left_join(module, by = c("name" = "ID")) %>%
      tidygraph::mutate(modularity = factor(Module),
                        modularity2 = factor(modularity),
                        modularity3 = as.character(modularity2),
                        Modularity = modularity2,
                        Degree = tidygraph::centrality_degree(mode = "out"),
                        Strength = tidygraph::centrality_degree(weights = weight)
      ) %>%
      tidygraph::arrange(Modularity, desc(Degree)) %>%
      tidygraph::left_join(node_annotation %>%
                             purrr::set_names(c("name", colnames(node_annotation)[-1])),
                           by = "name")
  }

  return(graph_obj)
}
