#' Create sample-level subgraph topology in parallel
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
#' @param parallel Logical (default = FALSE). Whether to enable parallel
#'   computation across samples.
#' @param n_workers Integer (default = NULL). Number of workers when
#'   \code{parallel = TRUE}. If \code{NULL}, use
#'   \code{future::availableCores() - 1}.
#' @param seed Integer (default = 1115). Random seed for reproducibility.
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
#' res <- get_sample_subgraph_topology_parallel(
#'   graph_obj = obj,
#'   mat       = mat,
#'   bootstrap = 10,
#'   parallel  = TRUE,
#'   n_workers = 2
#' )
#' head(res$topology)
#' }
get_sample_subgraph_topology_parallel <- function(graph_obj,
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
                                                  parallel = FALSE,
                                                  n_workers = NULL,
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

  if (isTRUE(parallel)) {
    if (is.null(n_workers)) {
      n_workers <- max(1, future::availableCores() - 1)
    }
    old_size <- getOption("future.globals.maxSize")
    old_plan <- future::plan()
    on.exit({
      options(future.globals.maxSize = old_size)
      future::plan(old_plan)
    }, add = TRUE)
    options(future.globals.maxSize = 6 * 1024^3)
    future::plan(future::multisession, workers = n_workers)
    # NOTE: progress display is provided by progressr::with_progress() around
    # the future_lapply() call below.  Do NOT call
    # `progressr::handlers(global = TRUE)` here: global handler registration
    # errors when the function is invoked inside tryCatch()/withCallingHandlers()
    # ("should not be called with handlers on the stack") and permanently
    # mutates the user's global progressr state.
  }
  # The sequential branch uses a plain lapply(), so no future::plan() setup is
  # needed there (and setting plan(sequential) would clobber the caller's plan
  # without restoring it).

  .compute_one_sample <- function(sid) {
    present_otu <- rownames(mat)[mat[, sid] != 0]
    present_otu <- intersect(present_otu, node_name)
    vids <- which(node_name %in% present_otu)

    if (length(vids) == 0) {
      return(list(
        sample = sid,
        subgraph = NULL,
        topology = NULL,
        robustness = NULL,
        stat = data.frame(
          Sample = sid,
          Node = 0,
          Edge = 0,
          Status = "No present OTUs in graph"
        )
      ))
    }

    ig_sub <- igraph::subgraph(ig, vids)
    graph_sub <- tidygraph::as_tbl_graph(ig_sub)
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
      return(list(
        sample = sid,
        subgraph = graph_sub,
        topology = NULL,
        robustness = NULL,
        stat = data.frame(
          Sample = sid,
          Node = igraph::vcount(ig_sub),
          Edge = igraph::ecount(ig_sub),
          Status = paste0("Failed: ", topo_res$message)
        )
      ))
    }

    list(
      sample = sid,
      subgraph = graph_sub,
      topology = topo_res$topology %>% dplyr::mutate(Sample = sid, .before = 1),
      robustness = topo_res$Robustness %>% dplyr::mutate(Sample = sid, .before = 1),
      stat = data.frame(
        Sample = sid,
        Node = igraph::vcount(ig_sub),
        Edge = igraph::ecount(ig_sub),
        Status = "OK"
      )
    )
  }

  if (isTRUE(parallel)) {
    progressr::with_progress({
      p <- progressr::progressor(steps = length(sample_ids))
      sample_res <- future.apply::future_lapply(sample_ids, function(sid) {
        out <- .compute_one_sample(sid)
        p()
        out
      }, future.seed = TRUE)
    })
  } else {
    sample_res <- lapply(sample_ids, .compute_one_sample)
  }

  names(sample_res) <- sample_ids

  subgraph_list <- sample_res %>%
    purrr::map("subgraph") %>%
    purrr::compact()

  topology_list <- sample_res %>%
    purrr::map("topology") %>%
    purrr::compact()

  robustness_list <- sample_res %>%
    purrr::map("robustness") %>%
    purrr::compact()

  stat_list <- sample_res %>%
    purrr::map("stat")

  topology_df <- if (length(topology_list) > 0) {
    dplyr::bind_rows(topology_list)
  } else {
    data.frame()
  }

  robustness_df <- if (length(robustness_list) > 0) {
    dplyr::bind_rows(robustness_list)
  } else {
    data.frame()
  }

  sample_stat_df <- dplyr::bind_rows(stat_list)

  list(
    subgraph_list = subgraph_list,
    topology = topology_df,
    Robustness = robustness_df,
    sample_stat = sample_stat_df
  )
}
