#' Build a correlation-based network from a matrix
#'
#' @param mat  Numeric matrix.
#' A numeric matrix with samples in colums and variables in rows
#' @param transfrom.method Character.
#'Data transformation methods applied before correlation analysis.
#' Options include:
#' "none" (raw data),
#' "scale" (z-score standardization),
#' "center" (mean centering only),
#' "log2" (log2 transfrom),
#' "log10" (log10 transfrom),
#' "ln" (natural transfrom ),
#' "rrarefy" (random rarefaction using \code{vegan::rrarefy}),
#' "rrarefy_relative" (rarefy then convert to relative abundance).
#' @param r.threshold Numeric.
#' Correlation coefficient threshold; edges are kept only if |r| >= r.threshold.
#' @param p.threshold Numeric.
#' Significance threshold for correlations; edges are kept only if p < p.threshold.
#' @param method Character.
#' Relationship analysis methods.
#' Options include: "WGCNA", "SpiecEasi", "SPARCC", "cor", and "Hmisc".
#' @param cor.method Character.
#' Correlation analysis method.
#' Options include "pearson", "kendall", and "spearman".
#' @param proc Character.
#' Correlation p-value adjustment methods.
#' Options include:
#' "holm", "hochberg", "hommel", "bonferroni",
#' "BH", "BY", "fdr", and "none".
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#' @param SpiecEasi.method Character.
#' Method used in \code{SpiecEasi} network inference; options include "mb" and "glasso".
#' @param sparcc_R Integer.
#' Number of bootstrap/permutation replicates for SparCC p-values (when \code{method = "SPARCC"}).
#' Default 20.
#' @param node_annotation Data frame.
#' Optional node annotation table, containing metadata such as taxonomy or functional categories.
#' @param top_modules Integer.
#' Number of top-ranked modules to retain for downstream visualization or analysis.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns A graph object representing the correlation-based microbial network.
#' Node/edge attributes include correlation statistics and (optionally) module labels.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' mat <- matrix(stats::rnorm(40 * 20), nrow = 40, ncol = 20)
#' rownames(mat) <- paste0("feature", seq_len(40))
#' colnames(mat) <- paste0("sample",  seq_len(20))
#' obj <- build_graph_from_mat(
#'   mat           = mat,
#'   method        = "cor",
#'   cor.method    = "pearson",
#'   proc          = "none",
#'   r.threshold   = 0.3,
#'   p.threshold   = 0.05,
#'   module.method = "Fast_greedy"
#' )
#' obj
#' }
build_graph_from_mat <- function(mat,
                                 transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy", "rrarefy_relative"),
                                 r.threshold = 0.7,
                                 p.threshold = 0.05,
                                 method = c("WGCNA", "SpiecEasi", "SPARCC", "cor", "Hmisc"),
                                 cor.method = c("pearson", "kendall", "spearman"),
                                 proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                 module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                 SpiecEasi.method = c("mb", "glasso"),
                                 sparcc_R = 20,
                                 node_annotation = NULL,
                                 top_modules = 15,
                                 seed = 1115){

  .ggnv_local_seed(seed)

  # argument check
  if (is.data.frame(mat)){
    mat <- as.matrix(mat)
    }

  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("`mat` must be numeric matrix.", call. = FALSE)
  }

  if (any(!is.finite(mat))) {
    stop("`mat` dont contain any NA/NaN/Inf. please check mat.", call. = FALSE)
  }

  if (is.null(colnames(mat))) {
    stop("`mat` must contains colnames.", call. = FALSE)
  }

  if (anyDuplicated(colnames(mat))) {
    dup <- unique(colnames(mat)[duplicated(colnames(mat))])
    stop(sprintf("`mat` must contain only colname. The duplicated colname: %s", paste(dup, collapse = ", ")), call. = FALSE)
  }

  if (length(top_modules) != 1L || !is.numeric(top_modules) || top_modules < 1) {
    stop("`top_modules` must be a single numeric value >= 1.", call. = FALSE)
  }
  top_modules <- as.integer(top_modules)
  if (length(seed) != 1L || !is.numeric(seed)) {
    stop("`seed` must be a single numeric.", call. = FALSE)
  }
  seed <- as.integer(seed)
  sparcc_R <- as.integer(sparcc_R)[1L]
  if (is.na(sparcc_R) || sparcc_R < 1L) {
    stop("`sparcc_R` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(node_annotation)) {
    if (!is.data.frame(node_annotation)) {
      stop("`annotation` must be a data.frame / tibble.", call. = FALSE)
    }
    if (ncol(node_annotation) < 2) {
      stop("`annotation` requires at least two columns (the first is the name, the rest are the annotation columns to be merged).", call. = FALSE)
    }
  }

  # argument check
  method <- match.arg(method)
  transfrom.method <-  match.arg(transfrom.method)
  cor.method <- match.arg(cor.method)
  proc <- match.arg(proc)
  module.method <- match.arg(module.method)
  SpiecEasi.method <- match.arg(SpiecEasi.method)

  # filter mat
  mat <- mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ID") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sum = sum(dplyr::c_across(where(is.numeric)))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(sum != 0) %>%
    dplyr::select(-sum) %>%
    tibble::column_to_rownames(var = "ID")

  # data transfrom
  mat <- apply_transform_method(mat, transfrom.method)

  # calculate correlation
  # Correct over the n(n-1)/2 unique off-diagonal tests, not all n^2 cells:
  # the duplicated triangle and the uninformative diagonal otherwise distort
  # the ranks that Benjamini-Hochberg depends on. See .ggnv_adjust_p_matrix().
  adjust_p_matrix <- function(p_mat, proc_method) {
    .ggnv_adjust_p_matrix(p_mat, proc_method)
  }

  # WGCNA
  if (method == "WGCNA") {
    # WGCNA for correlation
    occor <- WGCNA::corAndPvalue(t(mat), method = cor.method)
    occor.p <- adjust_p_matrix(occor$p, proc)

    # R and pvalue
    occor.r <- occor$cor
    diag(occor.r) <- 0
    occor.r[is.na(occor.p) | is.na(occor.r) | occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
    occor.r[is.na(occor.r)]=0

    # create igraph object
    g <- igraph::graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  }

  # SpiecEasi selects edges via StARS model selection (sparse inverse
  # covariance), NOT a correlation threshold or p-value. We therefore take the
  # weighted estimate directly -- partial correlation (glasso) or symmetrised
  # beta (mb) -- whose non-zero support IS the StARS-selected edge set and whose
  # values carry signed interaction strength. `r.threshold` / `p.threshold` do
  # not apply to this method.
  if (method == "SpiecEasi") {
    # mat: ASV x samples -> t(mat): samples x taxa for spieceasi_matrix_rcpp
    se_output <- if (SpiecEasi.method == "mb") "beta" else "partial_correlation"
    am <- spieceasi_matrix_rcpp(as.matrix(t(mat)), method = SpiecEasi.method, output = se_output,
                                lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 50))
    rownames(am) <- rownames(mat)
    colnames(am) <- rownames(mat)
    am[is.na(am)] <- 0
    diag(am) <- 0
    g <- igraph::graph_from_adjacency_matrix(am, weighted = TRUE, mode = 'undirected')
  }

  # SparCC (sparcc_matrix_rcpp: filter by r.threshold and p.threshold)
  if (method == "SPARCC") {
    # mat: ASV x samples -> t(mat): samples x taxa for sparcc_matrix_rcpp
    occor.r <- sparcc_matrix_rcpp(as.matrix(t(mat)))
    p_mat <- sparcc_pvalue_rcpp(as.matrix(t(mat)), R = sparcc_R, seed = seed)
    diag(occor.r) <- 0
    occor.r[abs(occor.r) < r.threshold | is.na(p_mat) | p_mat > p.threshold] <- 0
    occor.r[is.na(occor.r)] <- 0
    rownames(occor.r) <- rownames(mat)
    colnames(occor.r) <- rownames(mat)
    occor.r <- Matrix::Matrix(occor.r, sparse = TRUE)
    g <- igraph::graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  }

  # cor
  if (method == "cor") {
    # WGCNA for correlation
    occor <- psych::corr.test(t(mat), method = cor.method, adjust = "none")
    occor.p <- adjust_p_matrix(occor$p, proc)

    # R and pvalue
    occor.r <- occor$r
    diag(occor.r) <- 0
    occor.r[is.na(occor.p) | is.na(occor.r) | occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
    occor.r[is.na(occor.r)]=0

    # create igraph object
    g <- igraph::graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  }

  # Hmisc::rcorr
  if (method == "Hmisc") {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop("`method = \"Hmisc\"` requires the 'Hmisc' package. ",
           "Install it with install.packages(\"Hmisc\").", call. = FALSE)
    }
    if (identical(cor.method, "kendall")) {
      stop("`method = 'Hmisc'` uses `Hmisc::rcorr()` and only supports `cor.method = 'pearson'` or `cor.method = 'spearman'`.", call. = FALSE)
    }

    rcorr_type <- switch(
      cor.method,
      pearson = "pearson",
      spearman = "spearman"
    )

    occor <- Hmisc::rcorr(t(mat), type = rcorr_type)
    occor.r <- occor$r
    occor.r[abs(occor.r) < r.threshold] <- 0

    occor.p <- adjust_p_matrix(occor$P, proc)
    occor.p[occor.p >= p.threshold] <- -1
    occor.p[occor.p < p.threshold & occor.p >= 0] <- 1
    occor.p[occor.p == -1] <- 0
    occor.p[is.na(occor.p)] <- 0

    z <- occor.r * occor.p
    z[is.na(z)] <- 0
    diag(z) <- 0

    g <- igraph::graph_from_adjacency_matrix(z, weighted = TRUE, mode = "undirected")
  }

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
