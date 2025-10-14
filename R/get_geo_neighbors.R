get_neighbors <- function(ly,
                          k = 5,
                          idx = NULL,
                          coord = NULL,
                          seed = 1115,
                          tol = 1e-12) {
  # ly: data.frame，至少包含 x, y 两列
  # k:  需要的邻居个数（包含中心点本身）
  # idx: 指定中心点行号（可选）
  # coord: 指定中心点坐标 c(x, y)（可选，不一定在 ly 中）
  # seed: 随机中心时用于可复现
  # tol:  判断零距离的容差（浮点误差保护）

  stopifnot(all(c("x","y") %in% names(ly)))
  stopifnot(is.numeric(ly$x), is.numeric(ly$y))
  n <- nrow(ly)
  if (n < 1) stop("ly requires at least one point.")

  set.seed(seed)

  # 1) 确定中心 (fx, fy)，以及与中心重合的点们
  # 如果 coord 不是空的
  if (!is.null(coord)) {
    fx <- coord[1]
    fy <- coord[2]
    self_ids <- which(abs(ly$x - fx) <= tol & abs(ly$y - fy) <= tol)
    idx_used <- if (length(self_ids) == 1) self_ids else NA_integer_
  } else {
    # 如果id 不是空的
    if (is.null(idx)) idx <- sample.int(n, 1)
    stopifnot(idx >= 1, idx <= n)
    fx <- ly$x[idx]; fy <- ly$y[idx]
    self_ids <- idx
    idx_used <- idx
  }

  # 2) 计算距离
  d <- sqrt((ly$x - fx)^2 + (ly$y - fy)^2)

  # 3) 候选集合：finite 距离
  cand <- which(is.finite(d))
  o <- cand[order(d[cand], cand)]  # 按距离升序（再按行号稳定排序）

  # 4) 先放中心/重合点，再补最近的其他点，确保总数 = k
  self_pick <- intersect(self_ids, o)
  others <- setdiff(o, self_pick)
  need_self <- length(self_pick)
  k_self <- min(need_self, k)
  k_other <- max(0, k - k_self)
  nn_idx <- c(utils::head(self_pick, k_self), utils::head(others, k_other))

  # 5) 输出
  neighbors <- data.frame(
    node = nn_idx,
    x    = ly$x[nn_idx],
    y    = ly$y[nn_idx],
    dist = d[nn_idx]
  )
  focal <- data.frame(
    node = idx_used,  # 如果用 coord 且多点重合，可能 NA
    x = fx,
    y = fy
  )

  if (nrow(neighbors) < k) {
    message(sprintf("Only %d neighbors returned (fewer than the requested k = %d).", nrow(neighbors), k))
  }

  list(focal = focal, neighbors = neighbors)
}

# 没有 split 参数的
# module_layout <- function(graph_obj, layout, center = T, idx = NULL, shrink = 1){
#   # 从graph对象中提取出数据
#   node_df <- graph_obj %>%
#     tidygraph::activate(nodes) %>%
#     tidygraph::as_tibble()
#
#   # 按照模块进行排序
#   node_df %>%
#     dplyr::count(modularity3, name = "size") %>%
#     dplyr::arrange(desc(size)) %>%
#     dplyr::mutate(modularity4 = factor(modularity3,
#                                        levels = c(setdiff(modularity3, "Others"), "Others"),
#                                        ordered = T)) %>%
#     dplyr::arrange(modularity4) %>%
#     dplyr::mutate(modularity4 = as.character(modularity4)) %>%
#     dplyr::pull(modularity4) -> mod_levels
#
#   node_df_sorted <- node_df %>%
#     tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
#     tidygraph::arrange(modularity3, desc(degree))
#
#   node_df_sorted_number <- node_df_sorted %>%
#     dplyr::count(modularity3)
#
#   graph_obj_sort <- graph_obj %>%
#     tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = T)) %>%
#     tidygraph::arrange(modularity3)
#
#   # 是否进行缩放
#   shrink_rings_global <- function(df, shrink = shrink) {
#     anchor <- c(df$x[1], df$y[1])
#     r <- sqrt((df$x - anchor[1])^2 + (df$y - anchor[2])^2)
#     # 第一个点不动
#     r[1] <- 0
#     r_new <- r * shrink
#     theta <- atan2(df$y - anchor[2], df$x - anchor[1])
#     df$x <- anchor[1] + r_new * cos(theta)
#     df$y <- anchor[2] + r_new * sin(theta)
#     df
#   }
#
#
#   if (isTRUE(center)) {
#     coord = c(0,0)
#   }else{
#     coord = NULL
#   }
#
#   neighbors_list <- list()
#
#   for (i in 1:dim(node_df_sorted_number)[1]) {
#     # print(i)
#
#     if (i == 1) {
#       # 最大的模块，我要制定他的原点
#       out <- get_neighbors(ly = layout,
#                            k = node_df_sorted_number$n[i],
#                            coord = coord,
#                            idx = idx)
#       # 真实的坐标
#       out_ly <- out$neighbors
#       neighbors_list[[i]] <- shrink_rings_global(out_ly %>% dplyr::select(2,3),
#                                                  shrink = shrink)
#
#       # ly_sub <- ly %>%
#       #   tibble::rownames_to_column(var = "node") %>%
#       #   dplyr::filter(!node %in% out_ly$node) %>%
#       #   tibble::column_to_rownames(var = "node")
#       ly_sub <- layout[-out_ly$node, , drop = FALSE]
#
#     }else if (i  == dim(node_df_sorted_number)[1]) {
#       neighbors_list[[i]] <- ly_sub
#
#     }else{
#       out <- get_neighbors(ly = ly_sub,
#                            k = node_df_sorted_number$n[i])
#       # 真实的坐标
#       out_ly <- out$neighbors
#       neighbors_list[[i]] <- shrink_rings_global(out_ly %>% dplyr::select(2,3),
#                                                  shrink = shrink)
#       ly_sub <- ly_sub[-out_ly$node, , drop = FALSE]
#     }
#   }
#
#   ly_final <- do.call(rbind, neighbors_list)
#
#   dim(ly_final)
#
#   return(list(layout = ly_final,
#               graph_obj = graph_obj_sort))
# }


# # 新增加了split参数
# module_layout <- function(graph_obj,
#                           layout,
#                           center = TRUE,
#                           idx = NULL,
#                           shrink = 1,
#                           split = 1) {
#
#   # 只能执行一个：split 优先
#   mode <- if (!isTRUE(all.equal(split, 1))) {
#     "split"
#   } else if (!isTRUE(all.equal(shrink, 1))) {
#     "shrink"
#   } else {
#     "none"
#   }
#
#     # 从graph对象中提取出数据
#     node_df <- graph_obj %>%
#       tidygraph::activate(nodes) %>%
#       tidygraph::as_tibble()
#
#     # 按照模块进行排序
#     node_df %>%
#       dplyr::count(modularity3, name = "size") %>%
#       dplyr::arrange(desc(size)) %>%
#       dplyr::mutate(modularity4 = factor(modularity3,
#                                          levels = c(setdiff(modularity3, "Others"), "Others"),
#                                          ordered = T)) %>%
#       dplyr::arrange(modularity4) %>%
#       dplyr::mutate(modularity4 = as.character(modularity4)) %>%
#       dplyr::pull(modularity4) -> mod_levels
#
#     node_df_sorted <- node_df %>%
#       tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
#       tidygraph::arrange(modularity3, desc(degree))
#
#     node_df_sorted_number <- node_df_sorted %>%
#       dplyr::count(modularity3)
#
#     graph_obj_sort <- graph_obj %>%
#       tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = T)) %>%
#       tidygraph::arrange(modularity3)
#
#   # 缩放函数1：绕首点（shrink 用）
#   scale_anchor_first <- function(df_xy, factor) {
#     df <- df_xy; names(df)[1:2] <- c("x","y")
#     anchor <- c(df$x[1], df$y[1])
#     r <- sqrt((df$x - anchor[1])^2 + (df$y - anchor[2])^2)
#     theta <- atan2(df$y - anchor[2], df$x - anchor[1])
#     r[1] <- 0
#     r_new <- r * factor
#     df$x <- anchor[1] + r_new * cos(theta)
#     df$y <- anchor[2] + r_new * sin(theta)
#     df
#   }
#
#   # 缩放函数2：绕质心（split 用）
#   scale_around_centroid <- function(df_xy, factor) {
#     df <- df_xy; names(df)[1:2] <- c("x","y")
#     cx <- mean(df$x); cy <- mean(df$y)
#     dx <- df$x - cx; dy <- df$y - cy
#     df$x <- cx + dx * factor
#     df$y <- cy + dy * factor
#     df
#   }
#
#
#   neighbors_list <- vector("list", length = nrow(node_df_sorted_number))
#   ly_sub <- layout; if (is.matrix(ly_sub)) ly_sub <- as.data.frame(ly_sub)
#   if (!all(c("x","y") %in% names(ly_sub))) names(ly_sub)[1:2] <- c("x","y")
#   coord <- if (isTRUE(center)) c(0,0) else NULL
#
#   for (i in seq_len(nrow(node_df_sorted_number))) {
#     k_i <- node_df_sorted_number$n[i]
#
#     if (i == 1) {
#       out <- get_neighbors(ly = ly_sub, k = k_i, coord = coord, idx = idx)
#       xy  <- out$neighbors %>% dplyr::select(2,3)
#
#       if (mode == "split") {
#         neighbors_list[[i]] <- scale_around_centroid(xy, factor = split)
#       } else if (mode == "shrink") {
#         neighbors_list[[i]] <- scale_anchor_first(xy, factor = shrink)
#       } else {
#         neighbors_list[[i]] <- xy  # 保持不变
#       }
#
#       ly_sub <- ly_sub[-out$neighbors$node, , drop = FALSE]
#
#     } else if (i == nrow(node_df_sorted_number)) {
#       xy <- ly_sub
#
#       if (mode == "split") {
#         neighbors_list[[i]] <- scale_around_centroid(xy, factor = split)
#       } else if (mode == "shrink") {
#         neighbors_list[[i]] <- scale_anchor_first(xy, factor = shrink)
#       } else {
#         neighbors_list[[i]] <- xy
#       }
#
#     } else {
#       out <- get_neighbors(ly = ly_sub, k = k_i)
#       xy  <- out$neighbors %>% dplyr::select(2,3)
#
#       if (mode == "split") {
#         neighbors_list[[i]] <- scale_around_centroid(xy, factor = split)
#       } else if (mode == "shrink") {
#         neighbors_list[[i]] <- scale_anchor_first(xy, factor = shrink)
#       } else {
#         neighbors_list[[i]] <- xy
#       }
#
#       ly_sub <- ly_sub[-out$neighbors$node, , drop = FALSE]
#     }
#   }
#
#   ly_final <- do.call(rbind, neighbors_list)
#   list(layout = ly_final, graph_obj = graph_obj_sort)
# }


module_layout <- function(graph_obj,
                          layout,
                          center = TRUE,
                          idx = NULL,
                          shrink = 1,
                          split = NULL,
                          seed = 1115){

  set.seed(seed)
  # --- 工具函数：围绕子集第一个点做整体缩放（保留你原逻辑） ---
  shrink_rings_global <- function(df, shrink = shrink) {
    anchor <- c(df$x[1], df$y[1])
    r <- sqrt((df$x - anchor[1])^2 + (df$y - anchor[2])^2)
    r[1]  <- 0
    r_new <- r * shrink
    theta <- atan2(df$y - anchor[2], df$x - anchor[1])
    df$x  <- anchor[1] + r_new * cos(theta)
    df$y  <- anchor[2] + r_new * sin(theta)
    df
  }

  # --- 沿径向整体平移（给定“要加/减的半径”） ---
  radial_offset <- function(df, delta){
    r <- sqrt(df$x^2 + df$y^2)
    theta <- atan2(df$y, df$x)
    r_new <- r + delta
    df$x  <- r_new * cos(theta)
    df$y  <- r_new * sin(theta)
    df
  }

  # --- 代表半径：用模块当前坐标的平均半径（也可换成 median 或 max）---
  rep_radius <- function(df){
    mean(sqrt(df$x^2 + df$y^2))
  }

  # 1) 取节点数据
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  # 2) 确定模块顺序（大到小，Others 最后）
  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels

  # 3) 模块内按度数排（仅用于统计数量）
  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))

  # 4) 每模块节点数
  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)

  # 5) 返回的图对象（只按模块顺序排）
  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  # 6) 第一个模块的锚点（get_neighbors 用）
  coord <- if (isTRUE(center)) c(0, 0) else NULL

  neighbors_list <- list()

  # --- 关键：维护“上一圈目标半径” ---
  prev_target_r <- NULL

  # 7) 逐模块分配
  for (i in 1:nrow(node_df_sorted_number)) {

    if (i == 1) {
      out <- get_neighbors(ly = layout,
                           k = node_df_sorted_number$n[i],
                           coord = coord,
                           idx = idx)
      out_ly <- out$neighbors
      coords <- out_ly %>% dplyr::select(x, y)

      if (is.null(split)) {
        # 仅 shrink
        coords <- shrink_rings_global(coords, shrink = shrink)
      } else {
        # “按圈推进”逻辑：
        # 1) 先拿到该模块当前代表半径
        curr_r <- rep_radius(coords)
        # 2) 第一圈的目标半径：上一圈不存在，就以“当前代表半径 + split”为起点
        target_r <- curr_r + split
        # 3) 需要的平移量（半径增量）
        delta <- target_r - curr_r
        coords <- radial_offset(coords, delta = delta)
        # 4) 记录“上一圈目标半径”
        prev_target_r <- target_r
      }

      neighbors_list[[i]] <- coords
      ly_sub <- layout[-out_ly$node, , drop = FALSE]

    } else if (i == nrow(node_df_sorted_number)) {

      if (is.null(split)) {
        neighbors_list[[i]] <- ly_sub
      } else {
        coords <- ly_sub %>% dplyr::select(x, y)
        curr_r <- rep_radius(coords)
        # 目标半径 = 上一圈目标半径 + split
        target_r <- prev_target_r + split
        delta <- target_r - curr_r
        coords <- radial_offset(coords, delta = delta)
        neighbors_list[[i]] <- coords
        prev_target_r <- target_r
      }

    } else {
      out <- get_neighbors(ly = ly_sub,
                           k = node_df_sorted_number$n[i])
      out_ly <- out$neighbors
      coords <- out_ly %>% dplyr::select(x, y)

      if (is.null(split)) {
        coords <- shrink_rings_global(coords, shrink = shrink)
      } else {
        curr_r <- rep_radius(coords)
        target_r <- prev_target_r + split
        delta <- target_r - curr_r
        coords <- radial_offset(coords, delta = delta)
        prev_target_r <- target_r
      }

      neighbors_list[[i]] <- coords
      ly_sub <- ly_sub[-out_ly$node, , drop = FALSE]
    }
  }

  ly_final <- do.call(rbind, neighbors_list)

  # combine

  graph_ly_final <- dplyr::bind_cols(
    ly_final,
    graph_obj_sort %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
  )

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort,
              graph_ly_final = graph_ly_final))
}

