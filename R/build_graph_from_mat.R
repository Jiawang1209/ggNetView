#' Build a graph object from raw matrix for correlation analysis
#'
#' @param mat Numeric matrix.
#' A matrix with variables in columns; pairwise correlations are computed.
#' @param r.threshold Numeric.
#' Correlation coefficient threshold; edges are kept only if |r| >= r.threshold.
#' @param p.threshold Numeric.
#' Significance threshold for correlations; edges are kept only if p < p.threshold.
#' @param method Character string.
#'   Community detection method; one of
#'   \code{"cluster_fast_greedy"}, \code{"cluster_walktrap"},
#'   \code{"cluster_edge_betweenness"}, \code{"cluster_spinglass"}.
#' @param top_modules Integer.
#' Number of top-ranked modules to retain
#' @param seed Integer.
#' Random seed for reproducibility; if `NULL`, no seed is set.
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
                                 method = "cluster_fast_greedy",
                                 top_modules = 15,
                                 seed = 1115){

  set.seed(seed)
  # WGCNA 计算相关性
  occor <- WGCNA::corAndPvalue(t(mat), method = 'pearson')
  mtadj <- multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
  adpcor <- mtadj$adjp[order(mtadj$index),2]
  occor.p <- matrix(adpcor, dim(t(mat))[2])

  # R and pvalue
  occor.r <- occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
  occor.r[is.na(occor.r)]=0

  # 构建igraph对象
  g <- igraph::graph.adjacency(occor.r, weighted = TRUE, mode = 'undirected')

  # 删除自相关
  g <- igraph::simplify(g)

  # 删除孤立节点
  g <- igraph::delete_vertices(g, which(igraph::degree(g)==0))

  ## 设置网络的weight，为计算模块性做准备
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  # 模块化
  if (method == "cluster_fast_greedy") {
    igraph::V(g)$modularity <- igraph::membership(igraph::cluster_fast_greedy(g))
    igraph::V(g)$modularity2 <- as.character(igraph::V(g)$modularity)
  }

  if (method == "cluster_walktrap") {
    igraph::V(g)$modularity <- igraph::membership(igraph::cluster_walktrap(g))
    igraph::V(g)$modularity2 <- as.character(igraph::V(g)$modularity)
  }

  if (method == "cluster_edge_betweenness") {
    igraph::V(g)$modularity <- igraph::membership(igraph::cluster_edge_betweenness(g))
    igraph::V(g)$modularity2 <- as.character(igraph::V(g)$modularity)
  }

  if (method == "cluster_spinglass") {
    igraph::V(g)$modularity <- igraph::membership(igraph::cluster_spinglass(g))
    igraph::V(g)$modularity2 <- as.character(igraph::V(g)$modularity)
  }

  table(igraph::V(g)$modularity2) %>% sort(., decreasing = T)

  # max model length
  max_model <- length(table(igraph::V(g)$modularity2) %>% sort(., decreasing = T))

  if (max_model < top_modules) {

    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:max_model] %>% names()

  }else if (max_model > top_modules) {

    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:top_modules] %>% names()
  }

  igraph::V(g)$modularity2 <- ifelse(igraph::V(g)$modularity2 %in% modularity_top_15, igraph::V(g)$modularity2, "Others")

  # 构建ggraph对象
  graph_obj <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::mutate(modularity = factor(modularity),
                      modularity2 = factor(modularity2),
                      modularity3 = as.character(modularity2),
                      degree = tidygraph::centrality_degree(mode = "out"),
                      strength = tidygraph::centrality_degree(weights = weight)
    ) %>%
    tidygraph::arrange(modularity2, desc(degree))

  return(graph_obj)
}
