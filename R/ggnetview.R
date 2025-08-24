ggNetView <- function(graph_obj, 
                      layout = NULL, 
                      node_add = 7, 
                      r = 1, 
                      center = T, 
                      idx = NULL, 
                      shrink = 1,
                      label = F,
                      add_outer = F
                      ){
  
  # # test
  # graph_obj = graph_obj
  # layout = "gephi"
  # node_add = 7
  # r = 1
  
  # 首先拿到布局函数
  func_name <- paste0("create_layout_", layout)
  lay_func <- get(func_name, envir = .GlobalEnv)  # 从全局环境找对应函数
  
  # 获取布局
  ly1 = lay_func(graph_obj = graph_obj, node_add = node_add, r = r)
  
  # 圆形布局 添加模块化 获取模块
  ly1_1 <- module_layout(graph_obj, layout = ly1, center = center, idx = idx, shrink = shrink)
  
  # 可视化结果
  
  # label = F add_outer = F
  if (isFALSE(label) & isFALSE(add_outer)) {
    
    p1_1 <- ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      geom_edge_link(alpha = 0.2, colour = "grey70") +
      geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                   '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                   '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                   '#bdbdbd'),
                        name = "modularity") +
      coord_equal(clip = "off") +
      theme_void() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }
  
  # label = F add_outer = T
  if (isFALSE(label) & isTRUE(add_outer)) {
    
    maskTable <- generateMask(dims= ly1_1[["layout"]],
                              clusters=ly1_1[["graph_obj"]] %>% activate(nodes) %>% as_tibble() %>% dplyr::pull(modularity3))
    
    p1_1 <- ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      geom_edge_link(alpha = 0.2, colour = "grey70") +
      geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                   '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                   '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                   '#bdbdbd'),
                        name = "modularity") +
      new_scale_fill() + 
      geom_polygon(data=maskTable %>% dplyr::filter(cluster != "Others"), 
                   mapping = aes(x = x, y = y, group=group, fill = group, color = group),
                   linewidth = 1.25, linetype = 2,
                   # color = "#000000",
                   alpha = 0.5,
                   show.legend = F) + 
      scale_color_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                    '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                    '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                    '#bdbdbd'),
                         name = "modularity") + 
      scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                   '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                   '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                   '#bdbdbd'),
                        name = "modularity") +
      
      coord_equal(clip = "off") +
      theme_void() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }
  
  
  
  return(p1_1)
}