#' Build a correlation-based network from a adjacency matrix
#'
#' @param adjacency_matrix Numeric matrix.
#' A numeric matrix with adjacency matrix.
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#' @param node_annotation Data frame.
#' Optional node annotation table, containing metadata such as taxonomy or functional categories.
#' @param top_modules Integer.
#' Number of top-ranked modules to retain for downstream visualization or analysis.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns A graph object representing the correlation-based microbial network.
#' Node/edge attributes include correlation statistics and (optionally) module labels.
#' @export
#'
#' @examples
#' data(adjacency_matrix_example)
#' set.seed(1)
#' idx <- sample(ncol(adjacency_matrix_example), 80)
#' adj <- adjacency_matrix_example[idx, idx]
#' obj <- build_graph_from_adj_mat(
#'   adjacency_matrix = adj,
#'   module.method    = "Fast_greedy",
#'   top_modules      = 5
#' )
#' obj
build_graph_from_adj_mat <- function(adjacency_matrix,
                                     module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                     node_annotation = NULL,
                                     top_modules = 15,
                                     seed = 1115){

  .ggnv_local_seed(seed)

  # argument check
  if (is.data.frame(adjacency_matrix)){
    adjacency_matrix <- as.matrix(adjacency_matrix)
  }

  if (!is.matrix(adjacency_matrix) || !is.numeric(adjacency_matrix)) {
    stop("`adjacency_matrix` must be numeric matrix.", call. = FALSE)
  }

  if (any(!is.finite(adjacency_matrix))) {
    stop("`adjacency_matrix` dont contain any NA/NaN/Inf. please check mat.", call. = FALSE)
  }

  if (is.null(colnames(adjacency_matrix))) {
    stop("`adjacency_matrix` must contains colnames.", call. = FALSE)
  }

  if (anyDuplicated(colnames(adjacency_matrix))) {
    dup <- unique(colnames(adjacency_matrix)[duplicated(colnames(adjacency_matrix))])
    stop(sprintf("`adjacency_matrix` must contain only colname. The duplicated colname: %s", paste(dup, collapse = ", ")), call. = FALSE)
  }

  if (length(top_modules) != 1L || !is.numeric(top_modules) || top_modules < 1) {
    stop("`top_modules` must be a single numeric value >= 1.", call. = FALSE)
  }
  top_modules <- as.integer(top_modules)
  if (length(seed) != 1L || !is.numeric(seed)) {
    stop("`seed` must be a single numeric.", call. = FALSE)
  }
  seed <- as.integer(seed)


  if (!is.null(node_annotation)) {
    if (!is.data.frame(node_annotation)) {
      stop("`annotation` must be a data.frame / tibble.", call. = FALSE)
    }
    if (ncol(node_annotation) < 2) {
      stop("`annotation` requires at least two columns (the first is the name, the rest are the annotation columns to be merged).", call. = FALSE)
    }
  }

  # argument check
  module.method <- match.arg(module.method)


  # create igraph object
  g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE, mode = 'undirected')

  # remove self correlation
  g <- igraph::simplify(g)

  # delect single node
  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))

  ## set weight
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  igraph::E(g)$corr_direction <- ifelse(igraph::E(g)$correlation > 0, "Positive", "Negative")

  # membership
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
    # no others

  }else if (max_model >= top_modules) {

    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(top_modules)] %>% names()
  }

  igraph::V(g)$modularity2 <- ifelse(igraph::V(g)$modularity2 %in% modularity_top_15, igraph::V(g)$modularity2, "Others")

  modularity_top_final <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% names()
  modularity_top_final <- c(setdiff(modularity_top_final, "Others"), "Others")

  if (is.null(node_annotation)) {
    # create ggraph_obj
    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::mutate(modularity = factor(modularity),
                        modularity2 = factor(modularity2, levels = modularity_top_final, ordered = TRUE),
                        modularity3 = as.character(modularity2),
                        Modularity = modularity2,
                        Degree = tidygraph::centrality_degree(mode = "out"),
                        Strength = tidygraph::centrality_degree(weights = weight)
      ) %>%
      tidygraph::arrange(Modularity, desc(Degree))
  }else{
    # create ggraph_obj
    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::mutate(modularity = factor(modularity),
                        # modularity2 = factor(modularity2),
                        modularity2 = factor(modularity2, levels = modularity_top_final, ordered = TRUE),
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
