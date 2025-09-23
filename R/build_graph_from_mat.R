#' Build a graph object from raw matrix for correlation analysis
#'
#' @param mat Numeric matrix.
#' A matrix with variables in columns; pairwise correlations are computed.
#' @param r.threshold Numeric.
#' Correlation coefficient threshold; edges are kept only if |r| >= r.threshold.
#' @param p.threshold Numeric.
#' Significance threshold for correlations; edges are kept only if p < p.threshold.
#' @param top_modules Integer.
#' Number of top-ranked modules to retain
#' @param seed Integer.
#' Random seed for reproducibility; if `NULL`, no seed is set.
#' @param cor.method Charecter.
#' Correlation analysis methods contains "pearson", "kendall", "spearman"
#' @param proc
#' Corelation adjust pvalue methods contains "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD","BH", "BY","ABH","TSBH"
#' @param module.method Character
#' Module analysis methods contains "Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"
#' @param annotation Data Frame
#' The annotation file of nodes in network
#'
#' @returns An graph object representing the correlation network.
#' Node/edge attributes include correlation statistics and (optionally) module labels.
#'
#' @export
#'
#' @examples
#' set.seed(1115)
#' mat1 <- matrix(rnorm(10000), nrow = 100)  # 100 samples x 100 variables
#' g <- build_graph_from_mat(
#'   mat = mat1,
#'   r.threshold = 0.7,
#'   p.threshold = 0.05,
#'   method = "cluster_fast_greedy",
#'   top_modules = 15,
#'   seed = 1115
#' )
#' g
build_graph_from_mat <- function(mat,
                                 r.threshold = 0.7,
                                 p.threshold = 0.05,
                                 cor.method = c("pearson", "kendall", "spearman"),
                                 proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD","BH", "BY","ABH","TSBH"),
                                 module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                 annotation = NULL,
                                 top_modules = 15,
                                 seed = 1115){
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

  module.method <- match.arg(module.method)

  allowed_proc <- c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY","ABH","TSBH")
  if (is.null(proc) || length(proc) < 1L) {
    stop("`proc` provides at least one multi-correlation methods.", call. = FALSE)
  }

  if (!all(proc %in% allowed_proc)) {
    bad <- proc[!proc %in% allowed_proc]
    stop(sprintf("`proc` contains unsupported method: %s。Supported: %s",
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
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("`progress` must be a single TRUE/FALSE.", call. = FALSE)
  }
  if (!is.null(annotation)) {
    if (!is.data.frame(annotation)) {
      stop("`annotation` must be a data.frame / tibble.", call. = FALSE)
    }
    if (ncol(annotation) < 2) {
      stop("`annotation` requires at least two columns (the first is the name, the rest are the annotation columns to be merged).", call. = FALSE)
    }
  }

  set.seed(seed)

  # WGCNA 计算相关性
  occor <- WGCNA::corAndPvalue(t(mat), method = 'pearson')
  mtadj <- multtest::mt.rawp2adjp(unlist(occor$p),proc=proc)
  adpcor <- mtadj$adjp[order(mtadj$index),2]
  occor.p <- matrix(adpcor, dim(t(mat))[2])
  # R and pvalue
  occor.r <- occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
  occor.r[is.na(occor.r)]=0

  # 构建igraph对象
  g <- igraph::graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')

  # 删除自相关
  g <- igraph::simplify(g)

  # 删除孤立节点
  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))

  ## 设置网络的weight，为计算模块性做准备
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  # 模块化
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

  if (is.null(annotation)) {
    # 构建ggraph对象
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
    # 构建ggraph对象
    graph_obj <- tidygraph::as_tbl_graph(g) %>%
      tidygraph::mutate(modularity = factor(modularity),
                        modularity2 = factor(modularity2),
                        modularity3 = as.character(modularity2),
                        Modularity = modularity2,
                        Degree = tidygraph::centrality_degree(mode = "out"),
                        Strength = tidygraph::centrality_degree(weights = weight)
                        ) %>%
      tidygraph::arrange(Modularity, desc(Degree)) %>%
      tidygraph::left_join(annotation %>%
                             purrr::set_names(c("name", colnames(annotation)[-1])),
                           by = "name")
  }

  return(graph_obj)
}
