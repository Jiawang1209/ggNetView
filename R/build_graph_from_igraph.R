#' Build a graph object from a graph object
#'
#' @param igraph a igraph object
#' @param use_existing_modules Logical (default = TRUE).
#' If \code{TRUE}, reuse module information already stored on igraph vertices when available.
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#' @param module_attr Character or NULL (default = NULL).
#' Optional vertex attribute name containing precomputed module labels. When \code{NULL},
#' the function will auto-detect one of \code{"Modularity"}, \code{"modularity2"},
#' \code{"modularity3"}, or \code{"modularity"}.
#' @param node_annotation Data Frame
#' The annotation file of nodes in network
#' @param top_modules Integer.
#' Number of top-ranked modules to select.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns An graph object
#' @export
#'
#' @examples
#' data(ppi_example)
#' ig <- igraph::graph_from_data_frame(
#'   d        = ppi_example$ppi,
#'   vertices = ppi_example$annotation,
#'   directed = FALSE
#' )
#' obj <- build_graph_from_igraph(igraph = ig, module.method = "Fast_greedy")
#' levels(get_graph_nodes(obj)$Modularity)
build_graph_from_igraph <- function(igraph,
                                    use_existing_modules = TRUE,
                                    module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                    module_attr = NULL,
                                    node_annotation = NULL,
                                    top_modules = 15,
                                    seed = 1115
                                    ){
  g <- igraph
  # check argument
  module.method <- match.arg(module.method)

  if (!inherits(g, "igraph")) {
    stop("`igraph` must be an igraph object.", call. = FALSE)
  }

  if (igraph::gorder(g) == 0) {
    stop("`igraph` must contain at least one node.", call. = FALSE)
  }

  vertex_attrs <- igraph::vertex_attr_names(g)
  module_candidates <- c("Modularity", "modularity2", "modularity3", "modularity")

  # If the user explicitly named `module_attr` but it is not on the graph, fail
  # loudly instead of silently re-running community detection — the latter
  # ignores the user's intent and was the source of a hidden bug.
  if (!is.null(module_attr)) {
    if (length(module_attr) != 1L || !is.character(module_attr) || is.na(module_attr)) {
      stop("`module_attr` must be a single non-NA character string naming a vertex attribute.",
           call. = FALSE)
    }
    if (!module_attr %in% vertex_attrs) {
      stop(
        "`module_attr = \"", module_attr,
        "\"` is not a vertex attribute on the supplied graph. ",
        "Available vertex attributes: ",
        if (length(vertex_attrs) == 0L) "<none>" else paste(vertex_attrs, collapse = ", "), ".",
        call. = FALSE
      )
    }
    detected_module_attr <- module_attr
  } else {
    detected_module_attr <- module_candidates[module_candidates %in% vertex_attrs][1]
  }

  use_existing_modules <- isTRUE(use_existing_modules) &&
    !is.na(detected_module_attr) &&
    length(detected_module_attr) == 1L &&
    detected_module_attr %in% vertex_attrs

  if (!"weight" %in% igraph::edge_attr_names(g)) {
    igraph::E(g)$weight <- 1
  }

  if (igraph::any_multiple(g) || any(igraph::which_loop(g))) {
    g <- igraph::simplify(
      g,
      remove.multiple = TRUE,
      remove.loops    = TRUE,
      edge.attr.comb  = list(weight = "sum", "first")
    )
  }

  if (use_existing_modules) {
    module_values <- igraph::vertex_attr(g, detected_module_attr)
    if (all(is.na(module_values))) {
      stop(sprintf("Vertex attribute `%s` exists but contains only NA module values.", detected_module_attr), call. = FALSE)
    }
    module_values <- as.character(module_values)
    module_values[is.na(module_values) | trimws(module_values) == ""] <- "Others"

    igraph::V(g)$modularity <- factor(module_values)
    igraph::V(g)$modularity2 <- module_values
  } else {
    # membership
    # honour the documented `seed`: community detection (e.g. Spinglass) is
    # stochastic, so seed before it to keep module assignment reproducible.
    .ggnv_local_seed(seed)
    membership_vec <- switch(
      module.method,
      Fast_greedy = igraph::membership(igraph::cluster_fast_greedy(g)),
      Walktrap = igraph::membership(igraph::cluster_walktrap(g)),
      Edge_betweenness = igraph::membership(igraph::cluster_edge_betweenness(g)),
      Spinglass = igraph::membership(igraph::cluster_spinglass(g))
    )

    igraph::V(g)$modularity  <- membership_vec
    igraph::V(g)$modularity2 <- as.character(membership_vec)
  }

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
