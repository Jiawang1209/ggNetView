#' Compute the Integrated Value of Influence (IVI) per node
#'
#' Computes the Integrated Value of Influence (IVI) of Salavaty et al.
#' (2020) and attaches the per-node score as a new vertex column on the
#' input `tbl_graph`. IVI integrates six centrality measures across
#' three levels -- local (degree, local H-index), semi-local
#' (ClusterRank), and global (betweenness, closeness, collective
#' influence) -- and is more robust than any single centrality for
#' identifying influential nodes in complex biological networks
#' (microbiome co-occurrence, gene co-expression, PPI, ...).
#'
#' This function is a thin wrapper over `influential::ivi()`. The
#' `influential` package is the reference implementation maintained by
#' the IVI paper's first author and is required at runtime; install
#' with `install.packages("influential")` if not already available.
#' Routing through the canonical implementation guarantees the IVI
#' values match the published algorithm exactly.
#'
#' @section Forwarded arguments and version compatibility:
#' Different `influential` releases have exposed slightly different
#' argument lists for `ivi()`. Recent releases (the one this package
#' is developed against) accept a `scale` string with values `"range"`
#' / `"z-scale"` / `"none"`; older releases used a `scaled` logical
#' instead. To stay forward- and backward-compatible, this wrapper
#' inspects `formals(influential::ivi)` at call time and translates
#' the user-facing `scale` argument to whichever signature the
#' installed version expects. Arguments the installed version does not
#' recognise are silently dropped from the call.
#'
#' @param graph_obj A `tbl_graph` produced by any
#'   `build_graph_from_*()` constructor.
#' @param weights Numeric vector or `NULL` (default).
#'   Edge weights to pass through to `influential::ivi()`. `NULL` uses
#'   the unweighted IVI definition. Pass the string `"weight"` to use
#'   the graph's existing `weight` edge attribute as IVI weights
#'   (interpreted as distances by `influential::ivi`).
#' @param mode Character (default `"all"`). Edge mode passed to
#'   `influential::ivi()`. One of `"all"`, `"in"`, `"out"`. `"all"` is
#'   the right choice for undirected graphs.
#' @param directed Logical (default `FALSE`). Whether the graph should
#'   be treated as directed. Passed to `influential::ivi()`.
#' @param d Integer (default `3`). Distance horizon used by the
#'   collective-influence component. Larger `d` looks farther from each
#'   node when scoring it but increases run-time roughly linearly. The
#'   IVI paper recommends `d = 3` as a good trade-off.
#' @param scale Character (default `"range"`). How to scale the IVI
#'   values returned by `influential::ivi()`:
#'   \describe{
#'     \item{`"range"`}{Normalise to `[1, 100]`. Use this when
#'       exploring a single network; lets you see the full spread of
#'       node influences.}
#'     \item{`"z-scale"`}{z-score standardisation. Use this when
#'       comparing IVI across multiple networks, or when you want a
#'       defensible numeric threshold (`z > 1.645` is a common cutoff
#'       for "significantly influential").}
#'     \item{`"none"`}{Return the raw IVI scores with no scaling.}
#'   }
#'   On older versions of `influential` that exposed the boolean
#'   `scaled` argument instead of `scale`, this wrapper automatically
#'   translates: `"none"` -> `scaled = FALSE`, anything else ->
#'   `scaled = TRUE`.
#' @param ncores Integer (default `1`). Number of parallel cores for
#'   `influential::ivi()`'s internal cluster-rank step. The default
#'   (`1`) forces serial execution -- this is deliberately conservative
#'   so the wrapper works inside `R CMD check` (which caps tests at 2
#'   cores via `_R_CHECK_LIMIT_CORES_`) and inside vignette builds
#'   (where leftover parallel clusters are a common source of "9
#'   simultaneous processes spawned" / "invalid connection" errors).
#'   For interactive production work on a real network, raising this
#'   to e.g. `4` or `parallel::detectCores() - 1` gives a meaningful
#'   speed-up.
#' @param overwrite Logical (default `TRUE`). If an `IVI` column already
#'   exists on the input graph, controls whether to overwrite it (silent
#'   overwrite when `TRUE`; warning + return unchanged when `FALSE`).
#'
#' @returns A `tbl_graph` whose node table is augmented with a single
#'   new column, `IVI`. Other vertex / edge columns are preserved
#'   verbatim. Larger `IVI` = more influential.
#'
#' @references
#' Salavaty, A., Ramialison, M., & Currie, P. D. (2020). Integrated
#' Value of Influence: An Integrative Method for the Identification of
#' the Most Influential Nodes within Networks. *Patterns*, 1(5),
#' 100052.
#'
#' @seealso [get_node_centrality()] for the underlying per-node
#'   centrality panel; [ggnetview_zipi()] for the complementary Zi-Pi
#'   role classification (Guimera & Amaral 2005).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # `obj` is a tbl_graph from any build_graph_from_*() constructor.
#' obj_aug <- get_node_ivi(obj)
#'
#' # Top 10 most influential nodes:
#' obj_aug %>%
#'   tidygraph::activate(nodes) %>%
#'   tidygraph::as_tibble() %>%
#'   dplyr::arrange(dplyr::desc(IVI)) %>%
#'   utils::head(10)
#'
#' # Map IVI to point fill in ggNetView(): the fill aesthetic uses a
#' # discrete scale, so bin the continuous IVI into ordered quartile
#' # factor levels first, then pass the bin column as `node_fill`.
#' obj_plot <- obj_aug %>%
#'   tidygraph::activate(nodes) %>%
#'   tidygraph::mutate(IVI_bin = cut(
#'     IVI,
#'     breaks = stats::quantile(IVI, probs = seq(0, 1, 0.25), na.rm = TRUE),
#'     labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
#'     include.lowest = TRUE,
#'     ordered_result = TRUE
#'   ))
#' ggNetView(obj_plot, layout = "fr", node_fill = "IVI_bin")
#'
#' # Use z-score scaling when comparing across networks or applying a
#' # threshold (e.g. "significantly influential" = z > 1.645):
#' obj_z <- get_node_ivi(obj, scale = "z-scale")
#' }
get_node_ivi <- function(
  graph_obj,
  weights = NULL,
  mode = c("all", "in", "out"),
  directed = FALSE,
  d = 3,
  scale = c("range", "z-scale", "none"),
  ncores = 1L,
  overwrite = TRUE
) {

  if (!inherits(graph_obj, "tbl_graph")) {
    stop("`graph_obj` must be a tbl_graph (e.g. produced by build_graph_from_mat()).",
         call. = FALSE)
  }

  if (!requireNamespace("influential", quietly = TRUE)) {
    stop(
      "`get_node_ivi()` requires the `influential` package. ",
      "Install with: install.packages(\"influential\").",
      call. = FALSE
    )
  }

  mode  <- match.arg(mode)
  scale <- match.arg(scale)

  ig <- tidygraph::as.igraph(graph_obj)
  if (igraph::vcount(ig) == 0L) {
    stop("`graph_obj` has zero vertices.", call. = FALSE)
  }

  # Resolve `weights = "weight"` shorthand into the actual edge-weight
  # vector. Anything else (NULL or numeric vector) is passed through
  # to influential::ivi() unchanged.
  resolved_weights <- weights
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% igraph::edge_attr_names(ig)) {
      stop(sprintf(
        "Edge attribute `%s` not found on the graph; available: %s",
        weights,
        paste(igraph::edge_attr_names(ig), collapse = ", ")
      ), call. = FALSE)
    }
    resolved_weights <- igraph::edge_attr(ig, weights)
  }

  # ---- Defensive call to influential::ivi() --------------------------------
  # The signature of influential::ivi() has shifted across releases:
  # newer versions take `scale` (a character: "range" / "z-scale" /
  # "none") while older versions took `scaled` (a logical). Inspect
  # the installed version's formals and dispatch accordingly so the
  # same user-facing API works across the whole package history.
  ivi_formals <- names(formals(influential::ivi))

  desired_args <- list(
    graph    = ig,
    weights  = resolved_weights,
    directed = directed,
    mode     = mode,
    d        = d
  )

  if ("scale" %in% ivi_formals) {
    # Modern API: pass the user's choice straight through.
    desired_args$scale <- scale
  } else if ("scaled" %in% ivi_formals) {
    # Legacy API: translate. "none" disables scaling, every other
    # option turns scaling on (the legacy ivi() didn't distinguish
    # between range-normalised and z-scaled output).
    desired_args$scaled <- (scale != "none")
    if (scale == "z-scale") {
      message(
        "Installed `influential` version exposes a boolean `scaled` ",
        "argument and does not support `scale = 'z-scale'`; falling ",
        "back to `scaled = TRUE` (range-normalised output)."
      )
    }
  }

  if ("verbose" %in% ivi_formals) {
    desired_args$verbose <- FALSE
  }
  # `influential::ivi()` defaults `ncores = "default"` which spawns
  # `parallel::detectCores() - 1` workers. That blows past the 2-core
  # cap R CMD check imposes via `_R_CHECK_LIMIT_CORES_`, and also
  # causes intermittent "invalid connection" failures during vignette
  # builds. Forcing the user-controlled `ncores` argument (default 1)
  # eliminates both classes of failure.
  if ("ncores" %in% ivi_formals) {
    desired_args$ncores <- ncores
  }

  call_args <- desired_args[intersect(names(desired_args), ivi_formals)]
  ivi_vals  <- do.call(influential::ivi, call_args)

  # ---- Align IVI values to graph_obj's node order --------------------------
  vnames <- igraph::V(ig)$name
  if (!is.null(vnames) && !is.null(names(ivi_vals)) &&
      all(vnames %in% names(ivi_vals))) {
    ivi_aligned <- as.numeric(ivi_vals[vnames])
  } else {
    ivi_aligned <- as.numeric(ivi_vals)
    if (length(ivi_aligned) != igraph::vcount(ig)) {
      stop(sprintf(
        "Unexpected length from influential::ivi(): got %d, expected %d.",
        length(ivi_aligned), igraph::vcount(ig)
      ), call. = FALSE)
    }
  }

  # Honour `overwrite = FALSE`: don't silently clobber an existing `IVI` column
  # (mirrors the contract of get_node_centrality()).
  existing_cols <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble() %>%
    names()
  if ("IVI" %in% existing_cols && !isTRUE(overwrite)) {
    warning("Column `IVI` already exists on the graph; returning it unchanged ",
            "(set `overwrite = TRUE` to replace).", call. = FALSE)
    return(graph_obj)
  }

  graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::mutate(IVI = ivi_aligned)
}
