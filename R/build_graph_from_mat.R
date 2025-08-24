build_graph_from_mat <- function(matrix,
                                 cor_cutoff = 0.7,
                                 pvalue_cutoff = 0.05,
                                 method = "cluster_fast_greedy",
                                 top_module = 15,
                                 seed = 2025){
  # test
  # matrix = otu_rare_relative
  # cor_cutoff = 0.75
  # pvalue_cutoff = 0.05
  # seed = 2025
  # top_module = 15
  # method = "cluster_fast_greedy"
  
  set.seed(seed)
  # WGCNA 计算相关性
  occor<-WGCNA::corAndPvalue(t(matrix), method = 'pearson')
  mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor, dim(t(matrix))[2])
  
  # R and pvalue
  occor.r<-occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p>pvalue_cutoff|abs(occor.r)<cor_cutoff] = 0
  occor.r[is.na(occor.r)]=0
  
  # 构建igraph对象
  g <-  graph.adjacency(occor.r, weighted = TRUE, mode = 'undirected')
  
  # 删除自相关
  g <- simplify(g)
  
  # 删除孤立节点
  g <- delete_vertices(g, which(degree(g)==0))
  
  ## 设置网络的weight，为计算模块性做准备
  E(g)$correlation <- E(g)$weight
  E(g)$weight <- abs(E(g)$weight)
  
  # 模块化
  if (method == "cluster_fast_greedy") {
    V(g)$modularity <- membership(cluster_fast_greedy(g))
    V(g)$modularity2 <- as.character(V(g)$modularity)
  }
  
  if (method == "cluster_walktrap") {
    V(g)$modularity <- membership(cluster_walktrap(g))
    V(g)$modularity2 <- as.character(V(g)$modularity)
  }
  
  if (method == "cluster_edge_betweenness") {
    V(g)$modularity <- membership(cluster_edge_betweenness(g))
    V(g)$modularity2 <- as.character(V(g)$modularity)
  }
  
  if (method == "cluster_spinglass") {
    V(g)$modularity <- membership(cluster_spinglass(g))
    V(g)$modularity2 <- as.character(V(g)$modularity)
  }
  
  table(V(g)$modularity2) %>% sort(., decreasing = T) 
  
  modularity_top_15 <- V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:top_module] %>% names()
  
  
  V(g)$modularity2 <- ifelse(V(g)$modularity2 %in% modularity_top_15, V(g)$modularity2, "Others")
  
  # 构建ggraph对象
  graph_obj <- as_tbl_graph(g) %>%
    tidygraph::mutate(modularity = factor(modularity),
                      modularity2 = factor(modularity2),
                      modularity3 = as.character(modularity2),
                      degree = centrality_degree(mode = "out"),
                      strength = centrality_degree(weights = weight) 
    ) %>%
    tidygraph::arrange(modularity2, desc(degree)) 
  
  return(graph_obj)
}