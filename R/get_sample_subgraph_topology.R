#' Create sample-level subgraph topology from a global network
#'
#' For each sample (each column in \code{mat}), OTUs with non-zero abundance
#' are treated as present (\code{OTU != 0}). The function extracts a subgraph
#' from \code{graph_obj} using \code{igraph::subgraph()}, then computes topology
#' for each sample-specific subgraph via \code{get_network_topology()}.
#'
#' @param graph_obj An graph object from \code{build_graph_from_mat} or
#'   \code{build_graph_from_df}.
#' @param mat Numeric matrix used to build \code{graph_obj}. Rows are OTUs and
#'   columns are samples.
#' @param transfrom.method Character. Passed to \code{get_network_topology()}.
#' @param r.threshold Numeric. Passed to \code{get_network_topology()}.
#' @param p.threshold Numeric. Passed to \code{get_network_topology()}.
#' @param method Character. Passed to \code{get_network_topology()}.
#' @param cor.method Character. Passed to \code{get_network_topology()}.
#' @param proc Character. Passed to \code{get_network_topology()}.
#' @param SpiecEasi.method Character. Passed to \code{get_network_topology()};
#'   one of \code{"mb"} or \code{"glasso"} for SpiecEasi inverse-covariance.
#' @param sparcc_R Integer. Passed to \code{get_network_topology()} for SparCC p-values. Default 20.
#' @param bootstrap Numeric (default = 100). Passed to
#'   \code{get_network_topology()}.
#' @param seed Integer (default = 1115). Random seed for reproducibility,
#'   consistent with \code{get_sample_subgraph_topology_parallel()}.
#'
#' @returns A list with:
#' \itemize{
#'   \item \code{subgraph_list}: sample-wise \code{tbl_graph} subgraphs
#'   \item \code{topology}: merged sample-wise topology table
#'   \item \code{Robustness}: merged sample-wise robustness table
#'   \item \code{sample_stat}: sample-wise node/edge counts and status
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # `mat` is the OTU/feature-by-sample matrix used to build `graph_obj`.
#' obj <- build_graph_from_mat(mat = mat, method = "cor")
#' res <- get_sample_subgraph_topology(
#'   graph_obj = obj,
#'   mat       = mat,
#'   bootstrap = 10
#' )
#' head(res$topology)
#' }
get_sample_subgraph_topology <- function(graph_obj,
                                         mat = NULL,
                                         transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy", "rrarefy_relative"),
                                         r.threshold = 0.7,
                                         p.threshold = 0.05,
                                         method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
                                         cor.method = c("pearson", "kendall", "spearman"),
                                         proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                         SpiecEasi.method = c("mb", "glasso"),
                                         sparcc_R = 20,
                                         bootstrap = 100,
                                         seed = 1115) {
  .ggnv_local_seed(seed)
  transfrom.method <- match.arg(transfrom.method)
  method <- match.arg(method)
  cor.method <- match.arg(cor.method)
  proc <- match.arg(proc)

  if (!inherits(graph_obj, "tbl_graph")) {
    stop("`graph_obj` must be a `tbl_graph` object.", call. = FALSE)
  }
  if (is.null(mat)) {
    stop("`mat` is required for sample-level subgraph topology.", call. = FALSE)
  }
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("`mat` must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(mat))) {
    stop("`mat` must contain rownames (OTU IDs).", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    stop("`mat` must contain colnames (sample IDs).", call. = FALSE)
  }

  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()
  if (!"name" %in% colnames(node_df)) {
    stop("`graph_obj` node data must contain `name`.", call. = FALSE)
  }

  ig <- tidygraph::as.igraph(graph_obj)
  node_name <- as.character(node_df$name)

  sample_ids <- colnames(mat)
  subgraph_list <- list()
  topology_list <- list()
  robustness_list <- list()
  sample_stat <- list()

  for (sid in sample_ids) {
    present_otu <- rownames(mat)[mat[, sid] != 0]
    present_otu <- intersect(present_otu, node_name)
    vids <- which(node_name %in% present_otu)

    if (length(vids) == 0) {
      sample_stat[[sid]] <- data.frame(
        Sample = sid,
        Node = 0,
        Edge = 0,
        Status = "No present OTUs in graph"
      )
      next
    }

    ig_sub <- igraph::subgraph(ig, vids)
    graph_sub <- tidygraph::as_tbl_graph(ig_sub)
    subgraph_list[[sid]] <- graph_sub

    mat_sub <- mat[present_otu, , drop = FALSE]

    topo_res <- tryCatch(
      get_network_topology(
        graph_obj = graph_sub,
        mat = mat_sub,
        transfrom.method = transfrom.method,
        r.threshold = r.threshold,
        p.threshold = p.threshold,
        method = method,
        cor.method = cor.method,
        proc = proc,
        SpiecEasi.method = SpiecEasi.method,
        sparcc_R = sparcc_R,
        bootstrap = bootstrap,
        seed = seed + match(sid, sample_ids) - 1L
      ),
      error = function(e) e
    )

    if (inherits(topo_res, "error")) {
      sample_stat[[sid]] <- data.frame(
        Sample = sid,
        Node = igraph::vcount(ig_sub),
        Edge = igraph::ecount(ig_sub),
        Status = paste0("Failed: ", topo_res$message)
      )
      next
    }

    topology_list[[sid]] <- topo_res$topology %>%
      dplyr::mutate(Sample = sid, .before = 1)

    robustness_list[[sid]] <- topo_res$Robustness %>%
      dplyr::mutate(Sample = sid, .before = 1)

    sample_stat[[sid]] <- data.frame(
      Sample = sid,
      Node = igraph::vcount(ig_sub),
      Edge = igraph::ecount(ig_sub),
      Status = "OK"
    )
  }

  topology_df <- if (length(topology_list) > 0) {
    do.call(rbind, topology_list)
  } else {
    data.frame()
  }

  robustness_df <- if (length(robustness_list) > 0) {
    do.call(rbind, robustness_list)
  } else {
    data.frame()
  }

  sample_stat_df <- if (length(sample_stat) > 0) {
    do.call(rbind, sample_stat)
  } else {
    data.frame()
  }

  return(list(
    subgraph_list = subgraph_list,
    topology = topology_df,
    Robustness = robustness_df,
    sample_stat = sample_stat_df
  ))
}
