#' Build a graph object from separate node and edge tables
#'
#' Provides an explicit two-table interface where the user supplies a node
#' table and an edge table separately. Compared with [build_graph_from_df()],
#' which takes the edge list as the primary input and treats the node
#' annotation as optional, here the node table is **required and
#' authoritative**: nodes that have no incident edges are kept in the
#' resulting graph (no orphan-node removal). This is the right entry point
#' when the user already has a curated vertex set and wants every node to
#' appear in the visualization, even if it ends up isolated after edge
#' filtering.
#'
#' @param node Data frame.
#'   Node table whose **first column is the node identifier** that matches
#'   the `from` / `to` values in `edge`. Remaining columns are attached as
#'   vertex attributes (taxonomy, gene metadata, ...).
#' @param edge Data frame.
#'   Edge list with at least two columns; the **first two columns are
#'   interpreted as `from` and `to`** (column names are not enforced).
#'   An optional `weight` column (any case) is used as edge weight; if no
#'   `weight`-named column is found and `ncol(edge) >= 3`, the third
#'   column is used.
#' @param directed Logical (default: `FALSE`). Whether edges between nodes
#'   are directed.
#' @param module.method Character. Network community detection method.
#'   Options: `"Fast_greedy"`, `"Walktrap"`, `"Edge_betweenness"`,
#'   `"Spinglass"`. `"Fast_greedy"` and `"Spinglass"` require an undirected
#'   graph.
#' @param top_modules Integer. Number of top-ranked modules (by node count)
#'   to retain; smaller modules are collapsed into `"Others"`.
#' @param seed Integer (default: `1115`). Random seed for reproducibility.
#'
#' @returns A `tbl_graph` object with vertex columns `name`, `modularity`,
#'   `modularity2`, `modularity3`, `Modularity`, `Degree`, `Segree`,
#'   `Strength` plus any extra columns supplied via `node`, and edge
#'   columns `weight` (positive), `correlation` (signed) plus any extra
#'   columns supplied via `edge`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' node <- data.frame(
#'   name  = c("A", "B", "C", "D", "E"),
#'   group = c("g1", "g1", "g2", "g2", "g3")
#' )
#' edge <- data.frame(
#'   from   = c("A", "B", "C"),
#'   to     = c("B", "C", "D"),
#'   weight = c(0.8, 0.6, 0.7)
#' )
#' # Note: node "E" is isolated. Unlike build_graph_from_df(), the result
#' # of build_graph_from_node_edge() retains "E".
#' obj <- build_graph_from_node_edge(node = node, edge = edge)
#' obj
#' }
build_graph_from_node_edge <- function(node,
                                       edge,
                                       directed = FALSE,
                                       module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                       top_modules = 15,
                                       seed = 1115) {

  module.method <- match.arg(module.method)

  # `cluster_fast_greedy()` and `cluster_spinglass()` only support undirected
  # graphs in igraph. Fail fast with a clear message instead of letting
  # igraph surface a cryptic error from inside the switch() below.
  if (isTRUE(directed) && module.method %in% c("Fast_greedy", "Spinglass")) {
    stop(sprintf(
      "module.method = '%s' only supports undirected graphs; set `directed = FALSE` or pick another method (e.g. 'Walktrap', 'Edge_betweenness').",
      module.method
    ), call. = FALSE)
  }

  .ggnv_local_seed(seed)

  # ---- Argument validation -------------------------------------------------
  if (missing(node) || is.null(node)) {
    stop("`node` is required: pass a data.frame whose first column is the node identifier.", call. = FALSE)
  }
  if (missing(edge) || is.null(edge)) {
    stop("`edge` is required: pass a data.frame with at least `from` and `to` columns.", call. = FALSE)
  }

  node <- as.data.frame(node)
  edge <- as.data.frame(edge)

  if (ncol(node) < 1L) {
    stop("`node` must have at least one column (the node identifier).", call. = FALSE)
  }
  if (ncol(edge) < 2L) {
    stop("`edge` must have at least 2 columns (from, to).", call. = FALSE)
  }
  if (anyDuplicated(node[[1L]])) {
    dup <- unique(node[[1L]][duplicated(node[[1L]])])
    stop(sprintf(
      "`node` first column must contain unique IDs. Duplicated: %s",
      paste(utils::head(dup, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  # ---- Normalize edge weight column name -----------------------------------
  # Mirrors build_graph_from_df(): tolerate any casing of "weight"; fall back
  # to the 3rd column if no weight-like column is present.
  w_idx <- which(tolower(names(edge)) == "weight")
  if (length(w_idx) >= 1L) {
    if (length(w_idx) > 1L) {
      warning(sprintf(
        "Multiple weight-like columns found in `edge` (%s); using the first one and renaming to 'weight'.",
        paste(names(edge)[w_idx], collapse = ", ")
      ))
    }
    names(edge)[w_idx[1L]] <- "weight"
  } else if (ncol(edge) >= 3L) {
    message(sprintf(
      "No column named 'weight' (any case) in `edge`; treating column 3 ('%s') as weight.",
      names(edge)[3L]
    ))
    names(edge)[3L] <- "weight"
  }

  # ---- Cross-check: every endpoint must exist in node[[1]] -----------------
  node_ids <- as.character(node[[1L]])
  edge_ids <- unique(c(as.character(edge[[1L]]), as.character(edge[[2L]])))
  missing_ids <- setdiff(edge_ids, node_ids)
  if (length(missing_ids) > 0L) {
    stop(sprintf(
      "%d node ID(s) appear in `edge` but not in `node[[1]]`. First few: %s",
      length(missing_ids),
      paste(utils::head(missing_ids, 5L), collapse = ", ")
    ), call. = FALSE)
  }

  # ---- Build igraph anchored on `node` (preserve isolated nodes) -----------
  g <- igraph::graph_from_data_frame(
    d        = edge,
    vertices = node,
    directed = directed
  )

  # IMPORTANT: deliberately skip
  #   `igraph::delete_vertices(g, which(igraph::degree(g) == 0))`
  # The whole point of this two-table interface is that the node table is
  # authoritative and isolated nodes must be retained.
  g <- igraph::simplify(g)

  # ---- Edge attributes ------------------------------------------------------
  if (is.null(igraph::E(g)$weight)) {
    igraph::E(g)$correlation <- 1
    igraph::E(g)$weight      <- 1
  } else {
    igraph::E(g)$correlation <- igraph::E(g)$weight
    igraph::E(g)$weight      <- abs(igraph::E(g)$weight)
  }

  # ---- Module detection -----------------------------------------------------
  # Edge-cluster algorithms can crash on a graph with zero edges. Guard with
  # a trivial "each node = its own module" fallback so that the rest of the
  # pipeline (factor levels, top_modules truncation, tbl_graph conversion)
  # still runs cleanly.
  if (igraph::ecount(g) == 0L) {
    membership_vec <- seq_len(igraph::vcount(g))
    names(membership_vec) <- igraph::V(g)$name
  } else {
    membership_vec <- switch(
      module.method,
      Fast_greedy      = igraph::membership(igraph::cluster_fast_greedy(g)),
      Walktrap         = igraph::membership(igraph::cluster_walktrap(g)),
      Edge_betweenness = igraph::membership(igraph::cluster_edge_betweenness(g)),
      Spinglass        = igraph::membership(igraph::cluster_spinglass(g))
    )
  }

  igraph::V(g)$modularity  <- membership_vec
  igraph::V(g)$modularity2 <- as.character(membership_vec)

  # ---- top_modules truncation ----------------------------------------------
  max_model <- length(table(igraph::V(g)$modularity2) %>% sort(., decreasing = TRUE))
  if (max_model < top_modules) {
    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(max_model)] %>% names()
  } else {
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(top_modules)] %>% names()
  }
  igraph::V(g)$modularity2 <- ifelse(
    igraph::V(g)$modularity2 %in% modularity_top_15,
    igraph::V(g)$modularity2,
    "Others"
  )

  # ---- Materialize tbl_graph -----------------------------------------------
  graph_obj <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::mutate(
      modularity  = factor(modularity),
      modularity2 = factor(modularity2),
      modularity3 = as.character(modularity2),
      Modularity  = modularity2,
      Degree      = tidygraph::centrality_degree(mode = "out"),
      Segree      = tidygraph::centrality_degree(mode = "out"),
      Strength    = tidygraph::centrality_degree(weights = weight)
    ) %>%
    tidygraph::arrange(modularity2, dplyr::desc(Degree))

  return(graph_obj)
}
