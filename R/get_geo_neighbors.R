get_neighbors <- function(ly, k = 5, idx = NULL, coord = NULL, seed = NULL, tol = 1e-12) {
  # ly: data.frame，至少包含 x, y 两列
  # k:  需要的邻居个数（包含中心点本身）
  # idx: 指定中心点行号（可选）
  # coord: 指定中心点坐标 c(x, y)（可选，不一定在 ly 中）
  # seed: 随机中心时用于可复现
  # tol:  判断零距离的容差（浮点误差保护）

  stopifnot(all(c("x","y") %in% names(ly)))
  stopifnot(is.numeric(ly$x), is.numeric(ly$y))
  n <- nrow(ly); if (n < 1) stop("ly requires at least one point.")
  if (!is.null(seed) && is.null(idx) && is.null(coord)) set.seed(seed)

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
  nn_idx <- c(head(self_pick, k_self), head(others, k_other))

  # 5) 输出
  neighbors <- data.frame(
    node = nn_idx,
    x    = ly$x[nn_idx],
    y    = ly$y[nn_idx],
    dist = d[nn_idx]
  )
  focal <- data.frame(
    node = idx_used,  # 如果用 coord 且多点重合，可能 NA
    x = fx, y = fy
  )

  if (nrow(neighbors) < k) {
    message(sprintf("Only %d neighbors returned (fewer than the requested k = %d).", nrow(neighbors), k))
  }

  list(focal = focal, neighbors = neighbors)
}


module_layout <- function(graph_obj, layout, center = T, idx = NULL, shrink = 1){
  # 从graph对象中提取出数据
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  # 按照模块进行排序
  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = T)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels

  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, desc(degree))

  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)

  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = T)) %>%
    tidygraph::arrange(modularity3)

  # 是否进行缩放
  shrink_rings_global <- function(df, shrink = shrink) {
    anchor <- c(df$x[1], df$y[1])
    r <- sqrt((df$x - anchor[1])^2 + (df$y - anchor[2])^2)
    # 第一个点不动
    r[1] <- 0
    r_new <- r * shrink
    theta <- atan2(df$y - anchor[2], df$x - anchor[1])
    df$x <- anchor[1] + r_new * cos(theta)
    df$y <- anchor[2] + r_new * sin(theta)
    df
  }


  if (isTRUE(center)) {
    coord = c(0,0)
  }else{
    coord = NULL
  }

  neighbors_list <- list()

  for (i in 1:dim(node_df_sorted_number)[1]) {
    # print(i)

    if (i == 1) {
      # 最大的模块，我要制定他的原点
      out <- get_neighbors(ly = layout,
                           k = node_df_sorted_number$n[i],
                           coord = coord,
                           idx = idx)
      # 真实的坐标
      out_ly <- out$neighbors
      neighbors_list[[i]] <- shrink_rings_global(out_ly %>% dplyr::select(2,3),
                                                 shrink = shrink)

      # ly_sub <- ly %>%
      #   tibble::rownames_to_column(var = "node") %>%
      #   dplyr::filter(!node %in% out_ly$node) %>%
      #   tibble::column_to_rownames(var = "node")
      ly_sub <- layout[-out_ly$node, , drop = FALSE]

    }else if (i  == dim(node_df_sorted_number)[1]) {
      neighbors_list[[i]] <- ly_sub

    }else{
      out <- get_neighbors(ly = ly_sub,
                           k = node_df_sorted_number$n[i])
      # 真实的坐标
      out_ly <- out$neighbors
      neighbors_list[[i]] <- shrink_rings_global(out_ly %>% dplyr::select(2,3),
                                                 shrink = shrink)
      ly_sub <- ly_sub[-out_ly$node, , drop = FALSE]
    }
  }

  ly_final <- do.call(rbind, neighbors_list)

  dim(ly_final)

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort))
}



