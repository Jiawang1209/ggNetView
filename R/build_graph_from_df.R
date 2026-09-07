#' Build a graph object from a data frame
#'
#' @param df Data frame.
#'   Edge list with columns \code{from}, \code{to}, and optionally \code{weight}.
#'   If \code{weight} is absent, an unweighted graph is constructed.
#' @param node_annotation Data Frame
#' The annotation file of nodes in network
#' @param directed Logical (default: \code{FALSE}).
#'   Whether edges between nodes are directed.
#' @param top_modules Integer.
#'   Number of top-ranked modules to select.
#' @param seed Integer (default = 1115).
#'   Random seed for reproducibility.
#' @param module.method Character
#'   Module analysis methods contains "Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"
#' @returns An graph object representing the correlation network.
#'   Node/edge attributes include correlation statistics and (optionally) module labels.
#' @export
#'
#' @examples
#' data(ppi_example)
#' obj <- build_graph_from_df(
#'   df              = ppi_example$ppi,
#'   node_annotation = ppi_example$annotation,
#'   directed        = FALSE,
#'   module.method   = "Fast_greedy",
#'   top_modules     = 5
#' )
#' obj
#'
build_graph_from_df <- function(df,
                                node_annotation = NULL,
                                directed = FALSE,
                                module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                top_modules = 15,
                                seed = 1115){

  module.method <- match.arg(module.method)

  # `cluster_fast_greedy()` and `cluster_spinglass()` only support undirected
  # graphs in igraph. Fail fast with a clear message instead of letting igraph
  # surface a cryptic error deep inside the switch() below.
  if (isTRUE(directed) && module.method %in% c("Fast_greedy", "Spinglass")) {
    stop(sprintf(
      "module.method = '%s' only supports undirected graphs; set `directed = FALSE` or pick another method (e.g. 'Walktrap', 'Edge_betweenness').",
      module.method
    ))
  }

  .ggnv_local_seed(seed)

  # ---- Normalize edge data frame column names ------------------------------
  # igraph::graph_from_data_frame() takes the 1st/2nd columns as from/to by
  # position (column names don't matter), but turns column 3+ into edge
  # attributes using their literal names. Downstream code in this function
  # (igraph::E(g)$weight, centrality_degree(weights = weight), ...) looks up
  # the attribute named exactly "weight", and igraph attribute names are
  # case-sensitive. So we tolerate any casing the user supplies.
  df <- as.data.frame(df)

  if (ncol(df) < 2) {
    stop("`df` must have at least 2 columns (from, to).")
  }

  w_idx <- which(tolower(names(df)) == "weight")
  if (length(w_idx) >= 1) {
    if (length(w_idx) > 1) {
      warning(sprintf(
        "Multiple weight-like columns found (%s); using the first one and renaming it to 'weight'.",
        paste(names(df)[w_idx], collapse = ", ")
      ))
    }
    names(df)[w_idx[1]] <- "weight"
  } else if (ncol(df) >= 3) {
    message(sprintf(
      "No column named 'weight' (any case) found; treating the 3rd column ('%s') as weight.",
      names(df)[3]
    ))
    names(df)[3] <- "weight"
  }
  # --------------------------------------------------------------------------

  g <- igraph::graph_from_data_frame(
    d = df,
    vertices = node_annotation,
    directed = directed
  )


  g <- igraph::simplify(g)


  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))


  # ---- Handle weighted vs. unweighted graphs ------------------------------
  # When `df` carries no weight column, igraph::E(g)$weight is NULL and the
  # subsequent abs()/centrality_degree(weights = weight) calls would fail.
  # Fall back to an unweighted graph (every edge weight = 1) in that case.
  if (is.null(igraph::E(g)$weight)) {
    igraph::E(g)$correlation <- 1
    igraph::E(g)$weight      <- 1
  } else {
    igraph::E(g)$correlation <- igraph::E(g)$weight
    igraph::E(g)$weight      <- abs(igraph::E(g)$weight)
  }
  # --------------------------------------------------------------------------


  membership_vec <- switch(
    module.method,
    Fast_greedy = igraph::membership(igraph::cluster_fast_greedy(g)),
    Walktrap = igraph::membership(igraph::cluster_walktrap(g)),
    Edge_betweenness = igraph::membership(igraph::cluster_edge_betweenness(g)),
    Spinglass = igraph::membership(igraph::cluster_spinglass(g))
  )

  igraph::V(g)$modularity  <- membership_vec
  igraph::V(g)$modularity2 <- as.character(membership_vec)

  table(igraph::V(g)$modularity2) %>% sort(., decreasing = TRUE)

  # max model length
  max_model <- length(table(igraph::V(g)$modularity2) %>% sort(., decreasing = TRUE))

  if (max_model < top_modules) {

    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(max_model)] %>% names()

  }else if (max_model >= top_modules) {

    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(top_modules)] %>% names()
  }

  igraph::V(g)$modularity2 <- ifelse(igraph::V(g)$modularity2 %in% modularity_top_15, igraph::V(g)$modularity2, "Others")


  graph_obj <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::mutate(modularity = factor(modularity),
                      modularity2 = factor(modularity2),
                      modularity3 = as.character(modularity2),
                      Modularity = modularity2,
                      Degree = tidygraph::centrality_degree(mode = "out"),
                      Segree = tidygraph::centrality_degree(mode = "out"),
                      Strength = tidygraph::centrality_degree(weights = weight)
    ) %>%
    tidygraph::arrange(modularity2, desc(Degree))

  return(graph_obj)
}
