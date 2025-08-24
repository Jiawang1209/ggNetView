build_graph_from_df <- function(df,
                                node_annotation = NULL,
                                directed = F,
                                method = "cluster_fast_greedy",
                                top_module = 15,
                                seed = 2025){
  # test
  df = NULL
  node_annotation = NULL
  directed = F
  top_module = 15
  
  
  set.seed(seed)
 
  # 构建igraph对象
  g <- graph_from_data_frame(
    d = df,
    vertices = node_annotation,
    directed = directed
  )
  
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
  
  # max model length
  max_model <- length(table(V(g)$modularity2) %>% sort(., decreasing = T))
  
  if (max_model < top_module) {
    
    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:top_module] %>% names()
    
  }else if (max_model < top_module) {
    
    modularity_top_15 <- V(g)$modularity2 %>% table() %>% sort(., decreasing = T) %>% .[1:top_module] %>% names()
  }
  
  
  
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