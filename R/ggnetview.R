#' Visualize a network with custom layouts
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#'   The network object to be visualized.
#' @param layout Character string.
#'   custom layouts; one of
#'   "gephi", "square", "square2", "petal",
#'   "petal2", "heart_centered","diamond",
#'   "star", "star_concentric","rectangle,
#'   "rightiso_layers"
#'
#' @param node_add  Integer (default = 7).
#'   Number of nodes to add in each layer of the layout.
#' @param r Numeric (default = 1).
#'   Radius increment for concentric or layered layouts.
#' @param center Logical (default = TRUE).
#'   Whether to place a node at the center of the layout.
#' @param idx Optional.
#'   Index of nodes to be emphasized or centered in the layout
#' @param shrink Numeric (default = 1).
#'   Shrinkage factor applied to the center points.
#' @param label Logical (default = FALSE).
#'   Whether to display node labels in the center points
#' @param add_outer  Logical (default = FALSE).
#'   Whether to add an outer circle/border around the layout.
#' @param orientation Character string.
#'   custom orientation; one of
#'   "up","down","left","right"
#' @param angle Integer  (default = 0).
#'  change  orientation angle
#' @param split Numeric (default = 1).
#'   split factor applied to the center points.
#' @param linealpha Integer  (default = 0.25).
#'  change  line alpha
#' @param linecolor   Character  (default = "grey70").
#'  change  line color
#' @param outerwidth Integer  (default = 1.25).
#'  change  outer linewidth
#' @param outerlinetype  Integer  (default = 2).
#'  change  outer linetype
#' @param outeralpha Integer  (default = 0.5).
#'  change  outer alpha
#'
#' @returns A ggplot object representing the network visualization.
#' @export
#'
#' @examples NULL
ggNetView <- function(graph_obj,
                      layout = NULL,
                      node_add = 7,
                      r = 1,
                      center = T,
                      idx = NULL,
                      shrink = 1,
                      split = 1,
                      label = F,
                      linealpha = 0.25,
                      linecolor = "grey70",
                      add_outer = F,
                      outerwidth = 1.25,
                      outerlinetype = 2,
                      outeralpha = 0.5,
                      orientation = c("up","down","left","right"),
                      angle = 0 # 在 orientation 基础上的微调（弧度）
                      ){

  # 首先拿到布局函数
  func_name <- paste0("create_layout_", layout)
  lay_func <- utils::getFromNamespace(func_name, "ggNetView")  # 从全局环境找对应函数

  # 获取布局
  ly1 = lay_func(graph_obj = graph_obj, node_add = node_add, r = r, orientation = orientation, angle = angle)

  # 圆形布局 添加模块化 获取模块
  ly1_1 <- module_layout(graph_obj,
                         layout = ly1,
                         center = center,
                         idx = idx,
                         shrink = shrink,
                         split = split)

  # 可视化结果

  # label = F add_outer = F
  if (isFALSE(label) & isFALSE(add_outer)) {

    p1_1 <- ggraph::ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      ggraph::geom_edge_link(alpha = linealpha, colour = linecolor) +
      ggraph::geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                   '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                   '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                   '#bdbdbd',
                                   '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                   '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                   '#ffff99','#b15928'),
                        name = "modularity") +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }

  # label = T add_outer = F
  if (isTRUE(label) & isFALSE(add_outer)) {

    p1_1 <- ggraph::ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      ggraph::geom_edge_link(alpha = linealpha, colour = linecolor) +
      ggraph::geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      ggrepel::geom_text_repel(data = ly1_1[["graph_ly_final"]] %>%
                                 dplyr::distinct(modularity3, .keep_all = T) %>%
                                 dplyr::filter(modularity3 != "Others"),
                               mapping = aes(x = x, y = y, label = paste0("Module", modularity3))
                                 ) +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                                 name = "modularity") +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }

  # label = F add_outer = T
  if (isFALSE(label) & isTRUE(add_outer)) {

    maskTable <- mascarade::generateMask(dims= ly1_1[["layout"]],
                              clusters=ly1_1[["graph_obj"]] %>%
                                tidygraph::activate(nodes) %>%
                                tidygraph::as_tibble() %>%
                                dplyr::pull(modularity3))

    p1_1 <- ggraph::ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      ggraph::geom_edge_link(alpha = linealpha, colour = linecolor) +
      ggraph::geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                        name = "modularity") +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data=maskTable %>% dplyr::filter(cluster != "Others"),
                   mapping = aes(x = x, y = y, group=group, fill = group, color = group),
                   linewidth = outerwidth,
                   linetype = outerlinetype,
                   alpha = outeralpha,
                   show.legend = F) +
      ggplot2::scale_color_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                             '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                             '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                             '#bdbdbd',
                                             '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                             '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                             '#ffff99','#b15928'),
                         name = "modularity") +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                        name = "modularity") +

      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }

  # label = T add_outer = T
  if (isTRUE(label) & isTRUE(add_outer)) {

    maskTable <- mascarade::generateMask(dims= ly1_1[["layout"]],
                                         clusters=ly1_1[["graph_obj"]] %>%
                                           tidygraph::activate(nodes) %>%
                                           tidygraph::as_tibble() %>%
                                           dplyr::pull(modularity3))

    p1_1 <- ggraph::ggraph(ly1_1[["graph_obj"]], layout = "manual", x = ly1_1[["layout"]]$x, y = ly1_1[["layout"]]$y) +
      ggraph::geom_edge_link(alpha = linealpha, colour = linecolor) +
      ggraph::geom_node_point(aes(fill = modularity2, size = degree), alpha = 0.9, shape = 21) +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                                 name = "modularity") +
      ggrepel::geom_text_repel(data = ly1_1[["graph_ly_final"]] %>%
                                 dplyr::distinct(modularity3, .keep_all = T) %>%
                                 dplyr::filter(modularity3 != "Others"),
                               mapping = aes(x = x, y = y, label = paste0("Module", modularity3))
      ) +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data=maskTable %>% dplyr::filter(cluster != "Others"),
                            mapping = aes(x = x, y = y, group=group, fill = group, color = group),
                            linewidth = outerwidth,
                            linetype = outerlinetype,
                            alpha = outeralpha,
                            show.legend = F) +
      ggplot2::scale_color_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                             '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                             '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                             '#bdbdbd',
                                             '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                             '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                             '#ffff99','#b15928'),
                                  name = "modularity") +
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                                 name = "modularity") +

      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = margin(1,1,1,1,"cm")
      )
  }

  return(p1_1)
}
