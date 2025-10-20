#' Build a correlation-based network from a matrix
#'
#' @param mat  Numeric matrix.
#' A numeric matrix with samples in rows and variables in columns.
#'
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
#'
#' @param r.threshold Numeric.
#' Correlation coefficient threshold; edges are kept only if |r| >= r.threshold.
#'
#' @param p.threshold
#' Significance threshold for correlations; edges are kept only if p < p.threshold.
#'
#' @param method Character.
#' Relationship analysis methods.
#' Options include: "WGCNA", "SpiecEasi", and "SPARCC".
#'
#' @param cor.method Character.
#' Correlation analysis method.
#' Options include "pearson", "kendall", and "spearman".
#'
#' @param proc Character.
#' Correlation p-value adjustment methods.
#' Options include:
#' "Bonferroni", "Holm", "Hochberg", "
#' SidakSS", "SidakSD","BH",
#' "BY", "ABH", and "TSBH".
#'
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#'
#' @param SpiecEasi.method Character.
#' Method used in \code{SpiecEasi} network inference; options include "mb" and "glasso".
#'
#' @param node_annotation Data frame.
#' Optional node annotation table, containing metadata such as taxonomy or functional categories.
#'
#' @param top_modules Integer.
#' Number of top-ranked modules to retain for downstream visualization or analysis.
#'
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility..
#'
#' @returns A graph object representing the correlation-based microbial network.
#' Node/edge attributes include correlation statistics and (optionally) module labels.
#'
#'
#' @export
#'
#' @examples NULL
build_graph_from_mat <- function(mat,
                                 transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy", "rrarefy_relative"),
                                 r.threshold = 0.7,
                                 p.threshold = 0.05,
                                 method = c("WGCNA", "SpiecEasi", "SPARCC"),
                                 cor.method = c("pearson", "kendall", "spearman"),
                                 proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD","BH", "BY","ABH","TSBH"),
                                 module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                 SpiecEasi.method = c("mb", "glasso"),
                                 node_annotation = NULL,
                                 top_modules = 15,
                                 seed = 1115){

  set.seed(seed)

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



  allowed_proc <- c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY","ABH","TSBH")
  if (is.null(proc) || length(proc) < 1L) {
    stop("`proc` provides at least one multi-correlation methods.", call. = FALSE)
  }

  if (!all(proc %in% allowed_proc)) {
    bad <- proc[!proc %in% allowed_proc]
    stop(sprintf("`proc` contains unsupported method: %s, Supported: %s",
                 paste(bad, collapse = ", "),
                 paste(allowed_proc, collapse = ", ")),
         call. = FALSE)
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
  transfrom.method <-  match.arg(transfrom.method)
  module.method <- match.arg(module.method)
  cor.method <- match.arg(cor.method)
  proc <- match.arg(proc)
  module.method <- match.arg(module.method)
  SpiecEasi.method <- match.arg(SpiecEasi.method)

  # data transfrom
  mat <- switch (
    transfrom.method,
    none = mat,
    scale = t(scale(t(mat), scale = T, center = T)),
    center = t(scale(t(mat), scale = F, center = T)),
    log2 = log2(mat + 1),
    log10 = log10(mat + 1),
    ln = log(mat + 1),
    rrarefy = t(vegan::rrarefy(t(mat), min(colSums(mat)))),
    rrarefy_relative = t(vegan::rrarefy(t(mat), min(colSums(mat)))) / colSums(t(vegan::rrarefy(t(mat), min(colSums(mat)))))
  )

  # calculate correlation

  # WGCNA
  if (method == "WGCNA") {
    # WGCNA for correlation
    occor <- WGCNA::corAndPvalue(t(mat), method = cor.method)
    mtadj <- multtest::mt.rawp2adjp(unlist(occor$p),proc=proc)
    adpcor <- mtadj$adjp[order(mtadj$index),2]
    occor.p <- matrix(adpcor, dim(t(mat))[2])

    # R and pvalue
    occor.r <- occor$cor
    diag(occor.r) <- 0
    occor.r[occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
    occor.r[is.na(occor.r)]=0

    # create igraph object
    g <- igraph::graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  }

  # SpiecEasi
  if (method == "SpiecEasi") {
    # SpiecEasi for correlation
    SpiecEasi_obj <- SpiecEasi::spiec.easi(as.matrix(t(mat)),
                                           method = SpiecEasi.method,
                                           lambda.min.ratio=1e-2,
                                           nlambda=20,
                                           pulsar.params=list(rep.num=50)
                                           )

    # return adjacency matrix
    am <- SpiecEasi::getRefit(SpiecEasi_obj)

    rownames(am) <- rownames(mat)
    colnames(am) <- rownames(mat)

    # Create igraph objects
    g <- igraph::graph_from_adjacency_matrix(am, weighted = TRUE, mode = 'undirected')
  }

  # SparCC
  if (method == "SparCC") {
    # Sparcc for correlation
    SparCC_obj <- SpiecEasi::sparcc(as.matrix(t(mat)))

    SparCC_graph <- abs(SparCC_obj$Cor) >= r.threshold

    diag(SparCC_graph) <- 0

    rownames(SparCC_graph) <- rownames(mat)
    colnames(SparCC_graph) <- rownames(mat)

    SparCC_graph <- Matrix::Matrix(SparCC_graph, sparse=TRUE)

    # Create igraph objects
    g <- igraph::graph_from_adjacency_matrix(SparCC_graph, weighted = TRUE, mode = 'undirected')
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

  table(igraph::V(g)$modularity2) %>% sort(., decreasing = T)

  # max model length
  max_model <- length(table(igraph::V(g)$modularity2) %>% sort(., decreasing = T))

  if (max_model < top_modules) {

    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:max_model] %>% names()

  }else if (max_model >= top_modules) {

    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:top_modules] %>% names()
  }

  igraph::V(g)$modularity2 <- ifelse(igraph::V(g)$modularity2 %in% modularity_top_15, igraph::V(g)$modularity2, "Others")

  if (is.null(node_annotation)) {
    # create ggraph_obj
    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::mutate(modularity = factor(modularity),
                        modularity2 = factor(modularity2),
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
                        modularity2 = factor(modularity2),
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
