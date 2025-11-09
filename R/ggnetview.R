#' Visualize a network with custom layouts
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#' The network object to be visualized.
#' @param layout Character string.
#' Custom layouts; one of "gephi", "square", "square2", "petal",
#' "petal2", "heart_centered","diamond", "star", "star_concentric","rectangle,
#' "rightiso_layers" etc.
#' @param node_add Integer (default = 7).
#' Number of nodes to add in each layer of the layout.
#' @param ring_n Numeric (default = 7)
#' Numbers of ring in rings layout.
#' @param r Numeric (default = 1).
#' Radius increment for concentric or layered layouts.
#' @param center Logical (default = TRUE).
#' Whether to place a node at the center of the layout.
#' @param idx Optional.
#' Index of nodes to be emphasized or centered in the layout
#' @param shrink Numeric (default = 1).
#' Shrinkage factor applied to the center points.
#' @param k_nn Numeric (default = 8).
#' Number of nearest neighbors used to build the local adjacency graph.
#' @param push_others_delta Numeric (default = 0).
#' Radial offset applied to the "Others" module to slightly
#' @param layout.module Character  (default = "random")
#‘ - random : modules are distributed more randomly and independently.
#' - adjacent : modules are positioned close to each other, minimizing inter-module gaps.
#'
#' @param shape Integer  (default = 21).
#' The point shape likely in ggplot2.
#' @param pointalpha Integer  (default = 1).
#' The point alpha
#' @param pointsize Vector (default =  c(1,10))
#' The point size rang.
#' @param pointstroke Integer  (default = 0.3).
#' @param group.by Character (default = "Modularity").
#' Change variable for group
#' @param jitter Logical (default = FALSE).
#' Whether to apply jitter to points.
#' @param jitter_sd  Integer  (default = 0.1).
#' The standard deviation of the jitter applied when `jitter = TRUE`.
#' @param mapping_line  Logical (default = FALSE).
#' Whether to mapping line in ggNetView.
#' @param linealpha  Integer  (default = 0.25).
#' Change  line alpha.
#' @param linecolor Character  (default = "grey70").
#' Change  line color.
#' @param label Logical (default = FALSE).
#' Whether to display node labels in the center points.
#' @param labelsize Integer  (default = 10).
#' Change Module label size.
#' @param labelsegmentsize Integer  (default = 1).
#' Change  label segment size.
#' @param labelsegmentalpha Integer  (default = 1).
#' Change  label segment alpha.
#' @param add_outer Logical (default = FALSE).
#' Whether to add an outer circle/border around the layout.
#' @param outerwidth Integer  (default = 1.25).
#' Change  outer linewidth.
#' @param outerlinetype Integer  (default = 2).
#' Change  outer linetype.
#' @param outeralpha Integer  (default = 0.5).
#' Change  outer alpha.
#' @param nodelabsize Integer  (default = 5).
#' Change  node label size.
#' @param remove Logical (default = FALSE).
#' Delect nodes that are not modules.
#' @param orientation Character string.
#' Custom orientation; one of "up","down","left","right".
#' @param angle Integer  (default = 0).
#' Change  orientation angle.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.

#'
#' @returns A ggplot object representing the network visualization.
#' @export
#'
#' @examples NULL
ggNetView <- function(graph_obj,
                      layout = NULL,
                      node_add = 7,
                      ring_n = NULL,
                      r = 1,
                      center = TRUE,
                      idx = NULL,
                      shrink = 1,
                      k_nn = 8,
                      push_others_delta = 0,
                      layout.module = c("random", "adjacent"),
                      shape = 21,
                      pointalpha = 1,
                      pointsize = c(1,10),
                      pointstroke = 0.3,
                      group.by = "Modularity",
                      jitter = FALSE,
                      jitter_sd = 0.1,
                      mapping_line = FALSE,
                      linealpha = 0.25,
                      linecolor = "grey70",
                      label = FALSE,
                      labelsize = 10,
                      labelsegmentsize = 1,
                      labelsegmentalpha = 1,
                      add_outer = FALSE,
                      outerwidth = 1.25,
                      outerlinetype = 2,
                      outeralpha = 0.5,
                      nodelabsize = 5,
                      remove = FALSE,
                      orientation = "up",
                      angle = 0, # 在 orientation 基础上的微调（弧度）
                      seed = 1115
                      ){

  set.seed(seed)
  # find layout function
  func_name <- paste0("create_layout_", layout)

  # find layout functions from ggNetView package
  lay_func <- utils::getFromNamespace(func_name, "ggNetView")

  # get real data ly1
  if (func_name == "create_layout_rings") {
    ly1 = lay_func(graph_obj = graph_obj,
                   r = r,
                   ring_n = ring_n,
                   orientation = orientation,
                   angle = angle)
  }else{
    ly1 = lay_func(graph_obj = graph_obj,
                   node_add = node_add,
                   r = r,
                   orientation = orientation,
                   angle = angle)
  }

  # Normal layout
  # 只有当我们需要模块的时候，我们才要获取模块的布局
  # 并且只有模块化之后，我们才有必要去remove无模块的节点
  if (group.by != "pie") {
    # 圆形布局 添加模块化 获取模块
    if (layout.module == "random") {
      ly1_1 <- module_layout(graph_obj,
                             layout = ly1,
                             center = center,
                             idx = idx,
                             shrink = shrink# ,
                             # seed = seed
      )
    }

    if (layout.module == "adjacent") {
      ly1_1 <- module_layout3(graph_obj,
                              layout = ly1,
                              center = center,
                              k_nn = k_nn,
                              push_others_delta = push_others_delta,
                              shrink = shrink# ,
                              # seed = seed
      )
    }


    module_number <- ly1_1$graph_ly_final$Modularity %>% as.character() %>% unique() %>% length()

    # module info
    if (module_number == 1) {
      module_info <-  ly1_1$graph_ly_final$Modularity %>% as.character() %>% unique()

      ly1_1[["graph_ly_final"]] <- ly1_1[["graph_ly_final"]] %>%
        dplyr::mutate(Modularity = as.character(Modularity)) %>%
        dplyr::mutate(Modularity = factor(Modularity))

      ly1_1[["graph_obj"]] <- ly1_1[["graph_obj"]] %>%
        tidygraph::mutate(Modularity = as.character(Modularity)) %>%
        tidygraph::mutate(Modularity = factor(Modularity, ordered = T))

    }else{
      module_info <- levels(ly1_1$graph_ly_final$Modularity)
      module_info <- module_info[module_info!="Others"]
    }

    # remove Others
    if (isFALSE(remove)) {
      ly1_1 <- ly1_1
    }else{
      ly1_1[["graph_ly_final"]] <- ly1_1[["graph_ly_final"]] %>%
        dplyr::filter(as.character(Modularity) %in% module_info)

      ly1_1[["graph_obj"]] <- ly1_1[["graph_obj"]] %>%
        tidygraph::filter(name %in% ly1_1[["graph_ly_final"]]$name)


      ly1_1[["layout"]] <- dplyr::select(ly1_1[["graph_ly_final"]], x, y)

      ly1_1[["ggplot_data"]] <- get_location(ly1_1[["graph_ly_final"]],
                                             ly1_1[["graph_obj"]])
    }

    .build_label_location <- function(){
      # compute label location
      xr <<- range(ly1_1[["layout"]]$x)
      yr <<- range(ly1_1[["layout"]]$y)
      x_mid <<- stats::median(ly1_1[["layout"]]$x)
      dx <<- diff(xr) * 0.12
      pad <<- dx * 1.2

      # label df location
      lab_df <<- ly1_1[["graph_ly_final"]] %>%
        dplyr::distinct(modularity3, .keep_all = T) %>%
        dplyr::filter(modularity3 != "Others") %>%
        dplyr::mutate(side = ifelse(x < x_mid, "left", "right")) %>%
        dplyr::group_by(side) %>%
        dplyr::arrange(y, .by_group = TRUE) %>%
        dplyr::mutate(
          y_rank   = dplyr::row_number(),
          y_target = scales::rescale(y_rank, to = yr),                    # 把 rank 均匀映射到全局 y 范围
          x_anchor = dplyr::if_else(side == "left", xr[1] - dx, xr[2] + dx),
          nudge_x  = x_anchor - x,                                # 横向把标签推到两侧锚点
          nudge_y  = y_target - y,                                # 纵向把标签均匀拉开
          hjust    = dplyr::if_else(side == "left", 1, 0)
        ) %>%
        dplyr::ungroup()
    }


    # outlier df location
    .build_mask_table <- function(){
      maskTable <- mascarade::generateMask(dims= ly1_1[["layout"]],
                              clusters=ly1_1[["graph_obj"]] %>%
                                tidygraph::activate(nodes) %>%
                                tidygraph::as_tibble() %>%
                                dplyr::pull(modularity3)
      )

      return(maskTable)
    }

    ####----Plot----####
    # base plot
    p1_1 <- ggplot2::ggplot()

    # line parameter
    if (isFALSE(mapping_line)) {
       p1_1 <- p1_1 +
        ggplot2::geom_segment(data = ly1_1[["ggplot_data"]][[2]],
                              mapping = ggplot2::aes(x = from_x,
                                                     xend = to_x,
                                                     y = from_y,
                                                     yend = to_y),
                              alpha = linealpha,
                              colour = linecolor) +
         # ggplot2::coord_fixed() +
         theme_ggnetview() +
         scale_fill_ggnetview(levels(ly1_1[["graph_ly_final"]]$Modularity))

    }else{
      p1_1 <- p1_1 +
        ggplot2::geom_segment(data = ly1_1[["ggplot_data"]][[2]],
                              mapping = ggplot2::aes(x = from_x,
                                                     xend = to_x,
                                                     y = from_y,
                                                     yend = to_y,
                                                     colour = corr_direction),
                              alpha = linealpha) +
        # ggplot2::coord_fixed() +
        theme_ggnetview() +
        scale_fill_ggnetview(levels(ly1_1[["graph_ly_final"]]$Modularity))
    }

    # point paramers
    if (isFALSE(jitter)) {
      p1_1 <- p1_1 +
        ggplot2::geom_point(data = ly1_1[["ggplot_data"]][[1]],
                            mapping = ggplot2::aes(x = x, y = y, fill = .data[[group.by]], size = Degree),
                            shape = shape,
                            alpha = pointalpha,
                            stroke = pointstroke) +
        ggplot2::scale_size(range = pointsize) +
        ggplot2::coord_fixed() +
        theme_ggnetview() +
        scale_fill_ggnetview(levels(ly1_1[["graph_ly_final"]]$Modularity))
    }else{
      p1_1 <- p1_1 +
        ggplot2::geom_jitter(data = ly1_1[["ggplot_data"]][[1]],
                             mapping = ggplot2::aes(x = x, y = y, fill = .data[[group.by]], size = Degree),
                             shape = shape,
                             alpha = pointalpha,
                             stroke = pointstroke,
                             position = ggplot2::position_jitter(width = jitter_sd, height = jitter_sd, seed = seed)) +
        ggplot2::scale_size(range = pointsize) +
        ggplot2::coord_fixed() +
        theme_ggnetview() +
        scale_fill_ggnetview(levels(ly1_1[["graph_ly_final"]]$Modularity))

    }

    # label = F add_outer = F
    if (isFALSE(label) & isFALSE(add_outer)) {
      p1_1 <- p1_1

    }

    # label = T add_outer = F
    if (isTRUE(label) & isFALSE(add_outer)) {

      .build_label_location()

      p1_1 <- p1_1 +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggrepel::geom_text_repel(data = lab_df,
                                 mapping = ggplot2::aes(x = x,
                                               y = y,
                                               label = paste0("Module", modularity3),
                                               color = .data[[group.by]]),
                                 size = labelsize,
                                 nudge_x = lab_df$nudge_x,
                                 nudge_y = lab_df$nudge_y,
                                 hjust   = lab_df$hjust,
                                 min.segment.length = 0,
                                 segment.size = labelsegmentsize,
                                 segment.alpha = labelsegmentalpha,
                                 max.overlaps = Inf,
                                 box.padding = 0.15,
                                 point.padding = 0.15,
                                 force = 0.05,
                                 show.legend = F
        ) +
        scale_fill_ggnetview(levels(lab_df$Modularity)) +
        scale_color_ggnetview(levels(lab_df$Modularity)) +
        ggplot2::coord_equal(clip = "off",
                             xlim = c(xr[1] - pad, xr[2] + pad),
                             ylim = yr) +
        theme_ggnetview()
    }

    # label = F add_outer = T
    if (isFALSE(label) & isTRUE(add_outer)) {

      maskTable <- .build_mask_table()

      maskTable <- maskTable %>% dplyr::mutate(cluster = factor(cluster, levels = levels(ly1_1[["graph_ly_final"]]$Modularity), ordered = T))

      p1_1 <- p1_1 +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggplot2::geom_polygon(data=maskTable %>%
                                dplyr::filter(cluster != "Others"),
                              mapping = ggplot2::aes(x = x, y = y, group=cluster, fill = cluster, color = cluster),
                              linewidth = outerwidth,
                              linetype = outerlinetype,
                              alpha = outeralpha,
                              show.legend = F) +
        scale_fill_ggnetview(levels(maskTable$cluster)) +
        scale_color_ggnetview(levels(maskTable$cluster)) +
        ggplot2::coord_equal(clip = "off") +
        theme_ggnetview()
    }

    # label = T add_outer = T
    if (isTRUE(label) & isTRUE(add_outer)) {

      .build_label_location()
      maskTable <- .build_mask_table()

      maskTable <- maskTable %>% dplyr::mutate(cluster = factor(cluster, levels = levels(ly1_1[["graph_ly_final"]]$Modularity), ordered = T))

      p1_1 <- p1_1 +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggrepel::geom_text_repel(data = lab_df,
                                 mapping = ggplot2::aes(x = x,
                                               y = y,
                                               label = paste0("Module", modularity3),
                                               color = modularity2),
                                 size = labelsize,
                                 nudge_x = lab_df$nudge_x,
                                 nudge_y = lab_df$nudge_y,
                                 hjust   = lab_df$hjust,
                                 min.segment.length = 0,
                                 segment.size = labelsegmentsize,
                                 segment.alpha = labelsegmentalpha,
                                 max.overlaps = Inf,
                                 box.padding = 0.15,
                                 point.padding = 0.15,
                                 force = 0.05,
                                 show.legend = F
        ) +
        scale_fill_ggnetview(levels(lab_df$Modularity)) +
        scale_color_ggnetview(levels(lab_df$Modularity)) +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggplot2::geom_polygon(data= maskTable %>% dplyr::filter(cluster != "Others"),
                              mapping = ggplot2::aes(x = x, y = y, group=cluster, fill = cluster, color = cluster),
                              linewidth = outerwidth,
                              linetype = outerlinetype,
                              alpha = outeralpha,
                              show.legend = F) +
        scale_fill_ggnetview(levels(maskTable$cluster), na_value = NA) +
        scale_color_ggnetview(levels(maskTable$cluster), na_value = NA) +
        ggplot2::coord_equal(clip = "off",
                             xlim = c(xr[1] - pad, xr[2] + pad),
                             ylim = yr) +
        theme_ggnetview()
    }


  }

  # specific layout dendrogram
  if (func_name == "create_layout_rings") {
    p1_1 <- ggraph::ggraph(ly1)  +
      ggraph::geom_edge_link(alpha = linealpha, colour = linecolor) +
      ggraph::geom_node_point(ggplot2::aes(fill = group, size = Degree), shape = 21, alpha = pointalpha) +
      ggraph::geom_node_text(
        ggplot2::aes(x = 1.1 * x,
            y = 1.1 * y,
            label = name,
            angle = -((-ggraph::node_angle(x, y) + 90) %% 180) + 90),
        size = nodelabsize, hjust = 'outward'
      ) +
      ggplot2::scale_shape_manual(values = 20:25) +
      # scale_edge_color_gradientn(colors = c("#74add1","#abd9e9","#ffffbf","#fdae61","#f46d43"))+
      ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                            '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                            '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                            '#bdbdbd',
                                            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                            '#ffff99','#b15928'),
                                 name = "Group") +
      ggraph::scale_edge_width(range = c(0.1, 1)) +
      ggplot2::coord_equal(clip = "off") +
      theme_ggnetview()

  }

  if (layout == "dendrogram") {
    p1_1 <- ggraph::ggraph(graph_obj,layout = layout, circular = TRUE) +
      ggraph::geom_node_point(ggplot2::aes(size=node_size, color=type),alpha=pointalpha) +
      ggraph::geom_edge_diagonal(ggplot2::aes(color = node1.node), alpha=linealpha) +
      ggraph::scale_edge_color_manual(values = c('#66c2a5','#fc8d62','#a6d854','#e78ac3')) +
      ggplot2::scale_color_manual(values = c('#66c2a5','#fc8d62','#a6d854','#e78ac3')) +
      ggplot2::scale_size(range = c(3,15)) +
      ggraph::geom_node_text(
        ggplot2::aes(
          x = 1.0175 * x,
          y = 1.0175 * y,
          label = node,
          angle = -((-ggraph::node_angle(x, y) + 90) %% 180) + 90,
          filter = leaf,
          color = type
        ),
        size = 2, hjust = 'outward'
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label=node,
            filter = !leaf,
            color = type),
        fontface="bold",
        size=3,
        family="sans"
      ) +
      ggplot2::coord_fixed(clip = "off") +
      theme_ggnetview()

    return(p1_1)
  }

  # specific layout pie
  if (group.by == "pie") {

    ly <- ggraph::create_layout(graph_obj, layout = layout)

    col_index_start = which(colnames(ly) == "name")
    col_index_end = which(colnames(ly) == ".ggraph.orig_index")
    col_index = colnames(ly)[(col_index_start+1) : (col_index_end -1)]

    # 然后开始可视化
    p1_1 <- ggraph::ggraph(ly, layout = "manual", x = ly[["x"]], y = ly[["y"]]) +
      ggraph::geom_edge_link(color = "#6baed6") +
      scatterpie::geom_scatterpie(
        data = ly,
        cols = col_index,
        colour = "#000000",
        pie_scale = 2
      ) +
      ggplot2::scale_fill_manual(values = c('#66c2a5','#fc8d62','#a6d854','#e78ac3')) +
      ggplot2::coord_fixed() +
      theme_ggnetview()

    return(p1_1)
  }

  return(p1_1)
}
