#' Create network topology
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#'   The network object to be visualized.
#'
#' @returns data frame of network topolog
#' @export
#'
#' @examples NULL
get_network_topology <- function(graph_obj){
  # 生成igraph对象
  ig <- tidygraph::as.igraph(graph_obj)

  # 计算常见的拓扑属性
  # 节点度
  degree_vals <- igraph::degree(ig, mode = "all")

  # 平均路径长度
  avg_path <- igraph::mean_distance(ig, directed = FALSE)

  # 直径
  dia <- igraph::diameter(ig, directed = FALSE)

  # 聚类系数 (全局/局部)
  clust_global <- igraph::transitivity(ig, type = "global")
  clust_local  <- igraph::transitivity(ig, type = "local")

  # 中介中心性
  betweenness_vals <- igraph::betweenness(ig)

  # 接近中心性
  closeness_vals <- igraph::closeness(ig)

  # 特征向量中心性
  eigen_vals <- igraph::eigen_centrality(ig)$vector

  # 模块度 / 社区划分
  fast_greedy_com <- igraph::cluster_fast_greedy(ig)
  fast_greedy_com_modularity <- igraph::modularity(fast_greedy_com)

  # walktrap_com <- igraph::cluster_walktrap(ig)
  # walktrap_com_modularity <- igraph::modularity(walktrap_com)
  #
  # edge_betweenness_com <- igraph::cluster_edge_betweenness(ig)
  # edge_betweenness_com_modularity <- igraph::modularity(edge_betweenness_com)
  #
  # spinglass_com <- igraph::cluster_spinglass(ig)
  # spinglass_com_modularity <- igraph::modularity(spinglass_com)


  out <- data.frame(
    degree_vals = degree_vals,
    avg_path = avg_path,
    dia = dia,
    clust_global = clust_global,
    clust_local = clust_local,
    betweenness_vals = betweenness_vals,
    closeness_vals = closeness_vals,
    eigen_vals = eigen_vals,
    fast_greedy_com_modularity = fast_greedy_com_modularity
    # walktrap_com_modularity = walktrap_com_modularity,
    # edge_betweenness_com_modularity = edge_betweenness_com_modularity,
    # spinglass_com_modularity = spinglass_com_modularity
  )

  return(out)


}
