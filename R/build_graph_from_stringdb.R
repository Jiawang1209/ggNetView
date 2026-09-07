#' Build a graph object from a STRING database export
#'
#' Loads a STRING-formatted protein-protein interaction file (the TSV
#' export with columns `node1, node2, ..., combined_score`) or a
#' user-supplied data frame following the same schema, and converts it
#' into a `tbl_graph` consistent with the rest of the
#' `build_graph_from_*()` family. All STRING evidence channels in the
#' input -- `homology`, `coexpression`, `experimentally_determined_interaction`,
#' `automated_textmining`, etc. -- are preserved as edge attributes so
#' downstream visualizations can colour or filter edges by evidence type.
#'
#' STRING's TSV exports use a header line that starts with `#node1`. The
#' leading `#` is stripped automatically; all other column names are kept
#' verbatim.
#'
#' @param stringdb Either a file path (`.tsv` / `.tsv.gz`) or a
#'   data.frame already loaded in memory. STRING-style headers (a leading
#'   `#` on the first column name) are handled automatically.
#' @param node1_col,node2_col,score_col Character. Column names to use as
#'   the edge `from`, `to`, and weight columns. Defaults match STRING's
#'   TSV layout (`"node1"`, `"node2"`, `"combined_score"`). Override if
#'   you want to weight by an individual evidence channel -- e.g.
#'   `score_col = "experimentally_determined_interaction"`.
#' @param score_threshold Numeric or `NULL` (default).
#'   If non-NULL, edges with `score_col < score_threshold` are dropped.
#'   STRING confidence cut-offs are typically `0.4` (medium), `0.7`
#'   (high), `0.9` (highest). Pass `NULL` to keep every edge in the input.
#' @param node_annotation Optional data frame; first column must match the
#'   protein/gene names appearing in `node1_col` / `node2_col`. If
#'   provided, these annotations are attached as vertex attributes via
#'   [build_graph_from_df()].
#' @param directed Logical (default: `FALSE`). STRING interactions are
#'   undirected by convention.
#' @param module.method Character. Network community detection method.
#'   Options: `"Fast_greedy"`, `"Walktrap"`, `"Edge_betweenness"`,
#'   `"Spinglass"`.
#' @param top_modules Integer. Number of top-ranked modules to retain;
#'   smaller modules are collapsed into `"Others"`.
#' @param seed Integer (default: `1115`). Random seed for reproducibility.
#'
#' @returns A `tbl_graph` object compatible with [ggNetView()]. Edge
#'   attributes include `weight` (the chosen score), `correlation` (signed
#'   copy of the score, equal to `weight` for STRING since scores are
#'   non-negative), and every other column present in the input
#'   (e.g. `combined_score`, `coexpression`, `homology`,
#'   `experimentally_determined_interaction`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From a STRING TSV file:
#' obj <- build_graph_from_stringdb(
#'   stringdb        = "string_interactions.tsv",
#'   score_threshold = 0.7
#' )
#'
#' # From an in-memory data.frame, no filtering:
#' df <- data.frame(
#'   node1          = c("AANAT", "AANAT", "ABCA1"),
#'   node2          = c("CRY1",  "TPH1",  "SIRT1"),
#'   combined_score = c(0.608,   0.675,   0.520)
#' )
#' obj <- build_graph_from_stringdb(stringdb = df)
#' }
build_graph_from_stringdb <- function(stringdb,
                                      node1_col = "node1",
                                      node2_col = "node2",
                                      score_col = "combined_score",
                                      score_threshold = NULL,
                                      node_annotation = NULL,
                                      directed = FALSE,
                                      module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                      top_modules = 15,
                                      seed = 1115) {

  module.method <- match.arg(module.method)

  # ---- Read input ----------------------------------------------------------
  if (is.character(stringdb) && length(stringdb) == 1L) {
    if (!file.exists(stringdb)) {
      stop(sprintf("File not found: %s", stringdb), call. = FALSE)
    }
    # `comment = ""` keeps readr from treating the leading `#` as a
    # comment marker and silently dropping the header. We strip it
    # ourselves below to keep the operation explicit.
    df <- readr::read_tsv(stringdb, show_col_types = FALSE, comment = "")
    df <- as.data.frame(df)
  } else if (is.data.frame(stringdb)) {
    df <- as.data.frame(stringdb)
  } else {
    stop("`stringdb` must be a file path or a data.frame.", call. = FALSE)
  }

  if (nrow(df) == 0L) {
    stop("`stringdb` contains zero rows.", call. = FALSE)
  }

  # ---- Strip leading `#` from first column name (STRING export style) ------
  # STRING TSV exports use `#node1` to mark the header line. We strip the
  # `#` (and any whitespace immediately after it) without touching any
  # other column name -- this preserves the user's original schema.
  if (ncol(df) >= 1L) {
    names(df)[1L] <- sub("^#\\s*", "", names(df)[1L])
  }

  # ---- Validate required columns -------------------------------------------
  required <- c(node1_col, node2_col, score_col)
  missing_required <- setdiff(required, names(df))
  if (length(missing_required) > 0L) {
    stop(sprintf(
      "Column(s) not found in `stringdb`: %s.\n  Available columns: %s",
      paste(missing_required, collapse = ", "),
      paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }

  # ---- Coerce scores to numeric and warn on suspicious scale ---------------
  # `as.numeric()` on a non-numeric column would normally emit a warning;
  # we suppress here because we re-validate below and surface a cleaner
  # error if the result is unusable.
  score_vals <- suppressWarnings(as.numeric(df[[score_col]]))
  if (all(is.na(score_vals))) {
    stop(sprintf(
      "`%s` cannot be coerced to numeric. Check the column or pick a different `score_col`.",
      score_col
    ), call. = FALSE)
  }

  # STRING TSV exports are normally in [0, 1]; the STRING REST API returns
  # scores in [0, 1000]. Detect the latter and surface a soft message -- we
  # do NOT auto-divide because we don't want to silently mutate the user's
  # data.
  score_max <- max(score_vals, na.rm = TRUE)
  if (score_max > 1) {
    message(sprintf(
      "Detected `%s` values > 1 (max = %.3f). STRING TSV exports are usually in [0, 1]. If your file came from the STRING REST API, you may want to divide by 1000 before passing it in.",
      score_col, score_max
    ))
  }

  # ---- Optional score filter -----------------------------------------------
  if (!is.null(score_threshold)) {
    if (!is.numeric(score_threshold) || length(score_threshold) != 1L) {
      stop("`score_threshold` must be a single numeric value or NULL.", call. = FALSE)
    }
    keep <- !is.na(score_vals) & score_vals >= score_threshold
    n_dropped <- sum(!keep)
    df <- df[keep, , drop = FALSE]
    score_vals <- score_vals[keep]
    if (n_dropped > 0L) {
      message(sprintf(
        "Dropped %d edges with `%s` < %s.",
        n_dropped, score_col, format(score_threshold)
      ))
    }
    if (nrow(df) == 0L) {
      stop("No edges left after filtering -- try a lower `score_threshold`.", call. = FALSE)
    }
  }

  # ---- Drop rows with non-numeric scores (post-filter) ---------------------
  # If `score_threshold` is NULL we haven't culled NA scores yet. Drop them
  # here so downstream weighted centrality calculations don't propagate NA.
  na_score <- is.na(score_vals)
  if (any(na_score)) {
    message(sprintf(
      "Dropped %d edges with non-numeric `%s` value(s).",
      sum(na_score), score_col
    ))
    df <- df[!na_score, , drop = FALSE]
    score_vals <- score_vals[!na_score]
  }
  if (nrow(df) == 0L) {
    stop("No usable edges remain after coercion / filtering.", call. = FALSE)
  }

  # ---- Build the edge data frame --------------------------------------------
  # Keep all original columns (with the user's original names) so STRING
  # evidence channels survive as edge attributes. Two transformations:
  #
  #   1. Move `node1_col` / `node2_col` to positions 1 / 2 -- `igraph` picks
  #      `from` / `to` positionally.
  #   2. Add a `weight` column equal to a numeric copy of `score_col` so
  #      downstream weighted centrality and module detection can find it
  #      by literal name. The original `score_col` stays in place under
  #      its own name (e.g. `combined_score`).
  others <- setdiff(names(df), c(node1_col, node2_col))
  df <- df[, c(node1_col, node2_col, others), drop = FALSE]
  df$weight <- score_vals

  # `cluster_fast_greedy()` and `cluster_spinglass()` only accept undirected
  # graphs. Fail fast with a clear message instead of letting igraph surface
  # a cryptic error inside the switch() below.
  if (isTRUE(directed) && module.method %in% c("Fast_greedy", "Spinglass")) {
    stop(sprintf(
      "module.method = '%s' only supports undirected graphs; set `directed = FALSE` or pick another method (e.g. 'Walktrap', 'Edge_betweenness').",
      module.method
    ), call. = FALSE)
  }

  .ggnv_local_seed(seed)

  # ---- Build the igraph object ---------------------------------------------
  # We deliberately do NOT call build_graph_from_df() here. That function
  # uses `igraph::simplify(g)` with the default `edge.attr.comb`
  # (`list(weight = "sum", "ignore")`), which silently drops every edge
  # attribute that is not named `weight`. For STRING data that is exactly
  # the wrong default -- it would erase the evidence channels (homology /
  # coexpression / experimentally_determined_interaction / ...) that this
  # function is specifically designed to preserve.
  #
  # We therefore inline the same pipeline build_graph_from_df() runs, but
  # with an explicit `edge.attr.comb = list(weight = "first", "first")` so
  # every edge attribute survives simplification. (`"first"` rather than
  # `"sum"` because STRING evidence values are not additive across
  # duplicate rows -- the typical STRING TSV has no duplicate (node1,
  # node2) pairs anyway, so this only matters in pathological inputs.)
  g <- igraph::graph_from_data_frame(
    d        = df,
    vertices = node_annotation,
    directed = directed
  )

  g <- igraph::simplify(
    g,
    remove.multiple = TRUE,
    remove.loops    = TRUE,
    edge.attr.comb  = list(weight = "first", "first")
  )

  g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

  # ---- Edge weight handling -------------------------------------------------
  if (is.null(igraph::E(g)$weight)) {
    igraph::E(g)$correlation <- 1
    igraph::E(g)$weight      <- 1
  } else {
    igraph::E(g)$correlation <- igraph::E(g)$weight
    igraph::E(g)$weight      <- abs(igraph::E(g)$weight)
  }

  # ---- Module detection -----------------------------------------------------
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
