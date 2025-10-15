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


module_layout <- function(graph_obj,
                          layout,
                          center = TRUE,
                          idx = NULL,
                          shrink = 1,
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

      coords <- shrink_rings_global(coords, shrink = shrink)

      neighbors_list[[i]] <- coords
      ly_sub <- layout[-out_ly$node, , drop = FALSE]

    } else if (i == nrow(node_df_sorted_number)) {
        neighbors_list[[i]] <- ly_sub
    } else {
      out <- get_neighbors(ly = ly_sub,
                           k = node_df_sorted_number$n[i])
      out_ly <- out$neighbors
      coords <- out_ly %>% dplyr::select(x, y)

      coords <- shrink_rings_global(coords, shrink = shrink)

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


module_layout2 <- function(graph_obj,
                           layout,
                           center = TRUE,
                           idx = NULL,
                           shrink = 1,
                           arrange_by_radius = TRUE,
                           push_others_out   = TRUE,
                           seed = 1115){

  set.seed(seed)

  # --- 校验 ---
  # stopifnot(tidygraph::is_tbl_graph(graph_obj))
  # stopifnot(is.data.frame(layout), all(c("x","y") %in% names(layout)))
  # stopifnot(is.numeric(shrink), length(shrink) == 1, shrink > 0)

  # === 后面保持你之前的主体逻辑 ===
  shrink_rings_global <- function(df, shrink){
    anchor <- c(df$x[1], df$y[1])
    r      <- sqrt((df$x - anchor[1])^2 + (df$y - anchor[2])^2)
    r[1]   <- 0
    r_new  <- r * shrink
    theta  <- atan2(df$y - anchor[2], df$x - anchor[1])
    df$x   <- anchor[1] + r_new * cos(theta)
    df$y   <- anchor[2] + r_new * sin(theta)
    df
  }
  radial_offset <- function(df, delta){
    r     <- sqrt(df$x^2 + df$y^2)
    theta <- atan2(df$y, df$x)
    r_new <- r + delta
    df$x  <- r_new * cos(theta)
    df$y  <- r_new * sin(theta)
    df
  }

  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels

  node_df_sorted <- node_df %>%
    dplyr::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    dplyr::arrange(modularity3, dplyr::desc(Degree))

  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3, name = "n")

  if (isTRUE(arrange_by_radius)) {
    layout <- layout %>%
      dplyr::mutate(.r = sqrt(x^2 + y^2)) %>%
      dplyr::arrange(.r) %>%
      dplyr::select(-.r)
  }
  ly_sub <- layout

  pieces <- list()
  coord  <- if (isTRUE(center)) c(0,0) else NULL

  for (i in seq_len(nrow(node_df_sorted_number))) {
    n_i   <- node_df_sorted_number$n[i]
    mod_i <- levels(node_df_sorted$modularity3)[i]

    nodes_i <- node_df_sorted %>% dplyr::filter(modularity3 == mod_i)

    if (isTRUE(arrange_by_radius)) {
      if (i == nrow(node_df_sorted_number)) {
        used <- ly_sub
        if (nrow(used) != n_i) stop(sprintf("最后一圈剩余坐标(%d)与模块大小(%d)不一致。", nrow(used), n_i))
      } else {
        used_idx <- seq_len(n_i)
        used     <- ly_sub[used_idx, , drop = FALSE]
        ly_sub   <- ly_sub[-used_idx, , drop = FALSE]
      }
    } else {
      if (i == 1) {
        out    <- get_neighbors(ly = ly_sub, k = n_i, coord = coord, idx = idx)
        used   <- out$neighbors
        ly_sub <- ly_sub[-used$node, , drop = FALSE]
      } else if (i == nrow(node_df_sorted_number)) {
        used <- ly_sub
        if (nrow(used) != n_i) stop(sprintf("最后一圈剩余坐标(%d)与模块大小(%d)不一致。", nrow(used), n_i))
      } else {
        out    <- get_neighbors(ly = ly_sub, k = n_i)
        used   <- out$neighbors
        ly_sub <- ly_sub[-used$node, , drop = FALSE]
      }
    }

    coords  <- used %>% dplyr::select(x, y) %>% shrink_rings_global(shrink = shrink)
    piece_i <- dplyr::bind_cols(coords, nodes_i)
    pieces[[i]] <- piece_i
  }

  graph_ly_final <- dplyr::bind_rows(pieces)

  if (isTRUE(push_others_out) && any(graph_ly_final$modularity3 == "Others")) {
    r_all     <- sqrt(graph_ly_final$x^2 + graph_ly_final$y^2)
    is_mod    <- graph_ly_final$modularity3 != "Others"
    r_mod_max <- if (any(is_mod)) max(r_all[is_mod], na.rm = TRUE) else 0

    others_idx <- which(!is_mod)
    if (length(others_idx) > 0) {
      df_o    <- graph_ly_final[others_idx, c("x","y")]
      curr_r  <- sqrt(df_o$x^2 + df_o$y^2)
      gap_rel <- 0.10
      gap_abs <- max(1e-8, r_mod_max * gap_rel)
      target_r <- pmax(curr_r, r_mod_max + gap_abs)
      delta    <- target_r - curr_r
      df_o2    <- radial_offset(df_o, delta = delta)
      graph_ly_final$x[others_idx] <- df_o2$x
      graph_ly_final$y[others_idx] <- df_o2$y
    }
  }

  ly_final       <- graph_ly_final %>% dplyr::select(x, y)
  # graph_obj_sort <- graph_ly_final %>% dplyr::select(-x, -y)
  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  list(
    layout         = ly_final,
    graph_obj      = graph_obj_sort,
    graph_ly_final = graph_ly_final
  )
}


