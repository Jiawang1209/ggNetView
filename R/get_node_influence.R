#' Propagate a virtual perturbation from source node(s) across the network
#'
#' Injects a virtual perturbation at one or more `source` nodes and lets
#' it spread along the (optionally signed) weighted edges, returning how
#' strongly every other node is affected. This is the abundance-influence
#' ("type 2") virtual-perturbation analysis: it treats edge weights as
#' interaction strengths and asks "if I nudge species A, how far and how
#' strongly does the ripple reach?".
#'
#' Propagation uses a Katz / random-walk-with-restart diffusion,
#' \eqn{influence = (I - \alpha W)^{-1} s}, where `W` is the
#' column-normalised (signed) adjacency, `s` places `delta` on the source
#' node(s), and `alpha` is the decay. This always converges (`alpha` is
#' capped below the inverse spectral radius) and is well-defined even on
#' disconnected graphs.
#'
#' @section Important interpretation note:
#' Correlation / co-occurrence networks encode **association, not
#' causation**, and give no edge direction. The score returned here is a
#' *structural influence estimate* -- a weighted measure of how reachable
#' each node is from the source -- and should be read as a qualitative
#' ranking, **not** as a quantitative ecological-dynamics prediction. For
#' a perturbation read with a (still approximate) mechanistic flavour, see
#' [press_perturbation()].
#'
#' @param graph_obj A `tbl_graph` from any `build_graph_from_*()`
#'   constructor.
#' @param source Character vector of node `name`s to perturb.
#' @param delta Numeric (default `1`). Magnitude of the injected
#'   perturbation placed on each source node.
#' @param alpha Numeric in `(0, 1)` (default `0.5`). Diffusion decay; the
#'   function automatically caps it just below `1 / spectral-radius(W)` so
#'   the series converges.
#' @param signed Logical (default `TRUE`). Use the signed `correlation`
#'   edge attribute (so anticorrelated neighbours receive negative
#'   influence). `FALSE` uses `|weight|` only.
#' @param drop_source Logical (default `TRUE`). Zero out the source
#'   node(s)' own influence in the returned column so the ranking reflects
#'   downstream spread only.
#' @param overwrite Logical (default `TRUE`). If an `Influence` column
#'   already exists on the input graph, controls whether to overwrite it
#'   (silent overwrite when `TRUE`; warning + return unchanged when `FALSE`).
#'
#' @returns The input `tbl_graph` with one new node column, `Influence`
#'   (signed when `signed = TRUE`). Larger magnitude = more strongly
#'   affected. Map it straight onto a figure with
#'   `ggNetView(..., node_fill = "Influence")`.
#'
#' @seealso [get_network_perturbation()] for structural attacks;
#'   [press_perturbation()] for the press-perturbation approximation.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(ppi_example)
#' obj <- build_graph_from_df(
#'   df              = ppi_example$ppi,
#'   node_annotation = ppi_example$annotation
#' )
#' src <- get_graph_nodes(obj)$name[1]
#' obj2 <- get_node_influence(obj, source = src)
#' obj2 %>%
#'   tidygraph::activate(nodes) %>%
#'   tidygraph::as_tibble() %>%
#'   dplyr::arrange(dplyr::desc(abs(Influence))) %>%
#'   utils::head(5)
#' }
get_node_influence <- function(
  graph_obj,
  source,
  delta       = 1,
  alpha       = 0.5,
  signed      = TRUE,
  drop_source = TRUE,
  overwrite   = TRUE
) {

  if (!inherits(graph_obj, "tbl_graph")) {
    stop("`graph_obj` must be a tbl_graph (e.g. produced by build_graph_from_mat()).",
         call. = FALSE)
  }
  if (missing(source) || length(source) == 0L) {
    stop("`source` must name at least one node.", call. = FALSE)
  }
  if (!(alpha > 0 && alpha < 1)) {
    stop("`alpha` must lie strictly in (0, 1).", call. = FALSE)
  }

  ig <- tidygraph::as.igraph(graph_obj)
  n  <- igraph::vcount(ig)
  if (n == 0L) stop("`graph_obj` has zero vertices.", call. = FALSE)
  vnames <- igraph::V(ig)$name
  if (is.null(vnames)) vnames <- as.character(seq_len(n))

  bad <- setdiff(source, vnames)
  if (length(bad) > 0L) {
    stop(sprintf("Unknown source node(s): %s", paste(bad, collapse = ", ")),
         call. = FALSE)
  }

  # signed (or absolute) weighted adjacency
  attr_use <- if (isTRUE(signed) &&
                  "correlation" %in% igraph::edge_attr_names(ig)) "correlation" else "weight"
  W <- as.matrix(igraph::as_adjacency_matrix(ig, attr = attr_use, sparse = FALSE))
  if (!isTRUE(signed)) W <- abs(W)
  W[!is.finite(W)] <- 0
  dimnames(W) <- list(vnames, vnames)

  # column-normalise by absolute column sums (keeps spectral radius bounded)
  cs <- colSums(abs(W)); cs[cs == 0] <- 1
  Wn <- sweep(W, 2L, cs, "/")

  # cap alpha just under 1 / spectral radius for guaranteed convergence
  sr <- max(Mod(eigen(Wn, only.values = TRUE)$values))
  a  <- if (sr > 0) min(alpha, 0.95 / sr) else alpha

  s <- numeric(n); names(s) <- vnames
  s[source] <- delta

  infl <- solve(diag(n) - a * Wn, s)
  names(infl) <- vnames
  if (isTRUE(drop_source)) infl[source] <- 0

  # Honour `overwrite = FALSE`: don't silently clobber an existing `Influence`
  # column (mirrors the contract of get_node_centrality()).
  existing_cols <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble() %>%
    names()
  if ("Influence" %in% existing_cols && !isTRUE(overwrite)) {
    warning("Column `Influence` already exists on the graph; returning it unchanged ",
            "(set `overwrite = TRUE` to replace).", call. = FALSE)
    return(graph_obj)
  }

  graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::mutate(Influence = as.numeric(infl[vnames]))
}
