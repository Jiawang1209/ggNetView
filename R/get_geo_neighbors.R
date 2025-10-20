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
    # 如果 coord 是空的
    # 如果id 不是空的
    if (is.null(idx)) idx <- sample.int(n, 1)
    stopifnot(idx >= 1, idx <= n)
    fx <- ly$x[idx]
    fy <- ly$y[idx]
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


# 仅对 get_neighbors 的输出做“局部半径升序”的重排
get_neighbors_hub <- function(ly, k=5, idx=NULL, coord=NULL, seed=1115, tol=1e-12) {
  out <- get_neighbors(ly = ly, k = k, idx = idx, coord = coord, seed = seed, tol = tol)
  nb  <- out$neighbors
  # 以该块点集的质心为“局部中心”
  cx  <- mean(nb$x); cy <- mean(nb$y)
  r   <- sqrt((nb$x - cx)^2 + (nb$y - cy)^2)
  ord <- order(r, nb$dist, nb$node)   # 半径优先，随后备份排序
  out$neighbors <- nb[ord, , drop = FALSE]
  out
}

shrink_rings_global <- function(df, shrink){
  cx <- mean(df$x); cy <- mean(df$y)
  r  <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
  th <- atan2(df$y - cy, df$x - cx)
  r2 <- r * shrink
  df$x <- cx + r2 * cos(th)
  df$y <- cy + r2 * sin(th)
  df
}

radial_offset <- function(df, delta){
  r <- sqrt(df$x^2 + df$y^2)
  th <- atan2(df$y, df$x)
  df$x <- (r + delta) * cos(th)
  df$y <- (r + delta) * sin(th)
  df
}



module_layout <- function(graph_obj,
                          layout,
                          center = TRUE,
                          idx = NULL,
                          shrink = 1,
                          seed = seed){

  set.seed(seed)

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


  # 7) 逐模块分配
  for (i in 1:nrow(node_df_sorted_number)) {

    if (i == 1) {
      out <- get_neighbors_hub(ly = layout,
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

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)

  # get location result

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort,
              graph_ly_final = graph_ly_final,
              ggplot_data = ggplot_data
              ))
}



module_layout2 <- function(graph_obj,
                           layout,
                           center = TRUE,
                           shrink = 0.9,
                           k_nn = 8,
                           push_others_delta = 0.2,
                           seed = seed) {

  stopifnot(all(c("x","y") %in% names(layout)))

  set.seed(seed)

  ## ---------- 1) layout 邻接：kNN 图 ----------
  xy <- as.matrix(layout[, c("x","y")])
  nn  <- FNN::get.knn(xy, k = k_nn)$nn.index
  adj <- lapply(seq_len(nrow(layout)), function(i) unique(stats::na.omit(nn[i,])))

  ## ---------- 2) 节点与模块顺序 ----------
  nodes_tb <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_order <- nodes_tb %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::pull(modularity4)

  nodes_sorted <- nodes_tb %>%
    dplyr::mutate(modularity3 = factor(modularity3, levels = mod_order)) %>%
    dplyr::arrange(modularity3, dplyr::desc(Degree))

  need_per_mod <- nodes_sorted %>%
    dplyr::count(modularity3)

  ## ---------- 3) 区域生长占位 ----------
  nL    <- nrow(layout)
  free  <- rep(TRUE, nL)            # 该 layout 点是否空闲
  claim <- rep(NA_character_, nL)   # 被哪个模块占用
  placed_coords <- setNames(vector("list", nrow(need_per_mod)),
                            as.character(need_per_mod$modularity3))

  # 当前已占区域的“边界空位”
  get_frontier <- function(){
    used_idx <- which(!free)
    if (length(used_idx) == 0) return(integer())
    uniq <- unique(unlist(adj[used_idx], use.names = FALSE))
    uniq[free[uniq]]
  }

  # 第一模块种子：尽量靠中心（取最内圈若干随机）
  pick_seed_first <- function(){
    if (center) {
      r <- sqrt(layout$x^2 + layout$y^2)
      cand <- order(r)[seq_len(max(10, ceiling(0.01 * length(r))))]
      sample(cand, 1)
    } else {
      sample(which(free), 1)
    }
  }

  # 后续模块种子：从 frontier 里按“更靠当前已占质心”的概率抽
  pick_seed_next <- function(){
    fr <- get_frontier()
    if (length(fr) == 0) return(sample(which(free), 1))
    cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
    dd <- (layout$x[fr] - cx)^2 + (layout$y[fr] - cy)^2
    prob <- (max(dd) + 1e-9 - dd); prob <- prob / sum(prob)
    sample(fr, 1, prob = prob)
  }

  # BFS 生长（带轻启发式）：从种子扩张到 quota
  grow_region <- function(seed_idx, quota){
    q <- seed_idx
    region <- integer()
    while (length(region) < quota && length(q) > 0) {
      i <- q[1]; q <- q[-1]
      if (!free[i]) next
      region <- c(region, i); free[i] <<- FALSE
      nb <- adj[[i]]
      nb <- nb[free[nb]]
      if (length(nb) == 0) next

      # 轻微偏好“向当前已占质心靠近”，但仍保留随机
      if (any(!is.na(claim))) {
        cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
      } else {
        cx <- 0; cy <- 0
      }
      sc  <- -((layout$x[nb]-cx)^2 + (layout$y[nb]-cy)^2)
      # ord <- order(sc, runif(length(nb)))
      ord <- order(sc, nb)
      nb  <- nb[ord]
      q   <- c(q, nb)
    }
    region
  }

  # 若 BFS 不足 quota，则从全局 free 里按到 region/seed 最近补齐
  fill_deficit <- function(region, quota, target_idx){
    need <- quota - length(region)
    if (need <= 0) return(region)
    free_idx <- which(free)
    if (length(free_idx) == 0) return(region)

    if (length(region) > 0) {
      tx <- mean(layout$x[region]); ty <- mean(layout$y[region])
    } else {
      tx <- layout$x[target_idx];    ty <- layout$y[target_idx]
    }
    dd <- (layout$x[free_idx] - tx)^2 + (layout$y[free_idx] - ty)^2
    pick <- head(free_idx[order(dd)], min(need, length(free_idx)))
    free[pick] <<- FALSE
    c(region, pick)
  }

  # 逐模块（Others 留到最后）
  for (mi in seq_len(nrow(need_per_mod))) {
    mod    <- as.character(need_per_mod$modularity3[mi])
    k_need <- need_per_mod$n[mi]
    if (mod == "Others") next

    if (sum(free) < k_need) {
      stop(sprintf("剩余 layout 点不足：模块 %s 需要 %d，但只剩 %d。",
                   mod, k_need, sum(free)))
    }

    seed_idx <- if (mi == 1) pick_seed_first() else pick_seed_next()
    region   <- grow_region(seed_idx, k_need)
    if (length(region) < k_need) {
      region <- fill_deficit(region, k_need, target_idx = seed_idx)
    }

    claim[region] <- mod
    coords <- layout[region, c("x","y"), drop = FALSE]
    coords <- shrink_rings_global(coords, shrink = shrink)

    # 模块内：Degree降序 ↔ 半径升序（hub居中）
    cx <- mean(coords$x); cy <- mean(coords$y)
    rr <- sqrt((coords$x - cx)^2 + (coords$y - cy)^2)
    ord_pts <- order(rr, coords$x, coords$y)

    nodes_this <- nodes_sorted %>%
      dplyr::filter(modularity3 == mod)

    placed_coords[[mod]] <- dplyr::bind_cols(coords[ord_pts, , drop=FALSE], nodes_this)
  }

  ## ---------- 4) Others：外缘分配 ----------
  if ("Others" %in% names(placed_coords)) {
    n_oth <- need_per_mod$n[need_per_mod$modularity3 == "Others"]
    if (n_oth > 0) {
      idx_free <- which(free)
      if (length(idx_free) < n_oth) {
        stop(sprintf("Others 需要 %d 个点，但仅剩 %d 个。", n_oth, length(idx_free)))
      }
      r_free <- sqrt(layout$x[idx_free]^2 + layout$y[idx_free]^2)
      take   <- idx_free[order(-r_free)][seq_len(n_oth)]
      free[take] <- FALSE

      coords <- layout[take, c("x","y")]
      coords <- radial_offset(coords, delta = push_others_delta)  # 稍微外移，包边

      nodes_oth <- nodes_sorted %>%
        dplyr::filter(modularity3 == "Others")

      placed_coords[["Others"]] <- dplyr::bind_cols(coords, nodes_oth)
    }
  }

  ## ---------- 5) 拼装返回 ----------
  graph_ly_final <- dplyr::bind_rows(placed_coords[as.character(need_per_mod$modularity3)])
  ly_final <- graph_ly_final[, c("x","y")]

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

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)


  list(
    layout = ly_final,
    graph_obj = graph_obj_sort,
    graph_ly_final = graph_ly_final,
    ggplot_data = ggplot_data
  )
}


module_layout3 <- function(graph_obj,
                           layout,                # data.frame(x, y)
                           center = TRUE,         # 第一个模块优先靠中心放
                           shrink = 0.9,          # 模块内轻度收紧
                           k_nn = 8,              # layout 邻接度（6~10合适）
                           push_others_delta = 0.2, # Others 外移量
                           seed = seed) {

  stopifnot(all(c("x","y") %in% names(layout)))
  set.seed(seed)

  ## ---------- 小工具 ----------
  # shrink_rings_global <- function(df, shrink){
  #   cx <- mean(df$x); cy <- mean(df$y)
  #   r  <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
  #   th <- atan2(df$y - cy, df$x - cx)
  #   r2 <- r * shrink
  #   df$x <- cx + r2 * cos(th)
  #   df$y <- cy + r2 * sin(th)
  #   df
  # }
  # radial_offset <- function(df, delta){
  #   r <- sqrt(df$x^2 + df$y^2)
  #   th <- atan2(df$y, df$x)
  #   df$x <- (r + delta) * cos(th)
  #   df$y <- (r + delta) * sin(th)
  #   df
  # }

  ## ---------- 1) layout 邻接：kNN 图 ----------
  xy <- as.matrix(layout[, c("x","y")])
  nn  <- FNN::get.knn(xy, k = k_nn)$nn.index
  adj <- lapply(seq_len(nrow(layout)), function(i) unique(stats::na.omit(nn[i,])))

  ## ---------- 2) 节点与模块顺序 ----------
  nodes_tb <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_order <- nodes_tb %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::pull(modularity4)

  nodes_sorted <- nodes_tb %>%
    dplyr::mutate(modularity3 = factor(modularity3, levels = mod_order)) %>%
    dplyr::arrange(modularity3, dplyr::desc(Degree))

  need_per_mod <- nodes_sorted %>%
    dplyr::count(modularity3)

  ## ---------- 3) 区域生长占位 ----------
  nL    <- nrow(layout)
  free  <- rep(TRUE, nL)            # 该 layout 点是否空闲
  claim <- rep(NA_character_, nL)   # 被哪个模块占用
  placed_coords <- setNames(vector("list", nrow(need_per_mod)),
                            as.character(need_per_mod$modularity3))

  # 当前已占区域的“边界空位”
  get_frontier <- function(){
    used_idx <- which(!free)
    if (length(used_idx) == 0) return(integer())
    uniq <- unique(unlist(adj[used_idx], use.names = FALSE))
    uniq[free[uniq]]
  }

  # 第一模块种子：尽量靠中心（取最内圈若干随机）
  pick_seed_first <- function(){
    if (center) {
      r <- sqrt(layout$x^2 + layout$y^2)
      cand <- order(r)[seq_len(max(10, ceiling(0.01 * length(r))))]
      sample(cand, 1)
    } else {
      sample(which(free), 1)
    }
  }

  # 后续模块种子：从 frontier 里按“更靠当前已占质心”的概率抽
  pick_seed_next <- function(){
    fr <- get_frontier()
    if (length(fr) == 0) return(sample(which(free), 1))
    cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
    dd <- (layout$x[fr] - cx)^2 + (layout$y[fr] - cy)^2
    prob <- (max(dd) + 1e-9 - dd); prob <- prob / sum(prob)
    sample(fr, 1, prob = prob)
  }

  # BFS 生长（带轻启发式）：从种子扩张到 quota
  grow_region <- function(seed_idx, quota){
    q <- seed_idx
    region <- integer()
    while (length(region) < quota && length(q) > 0) {
      i <- q[1]; q <- q[-1]
      if (!free[i]) next
      region <- c(region, i); free[i] <<- FALSE
      nb <- adj[[i]]
      nb <- nb[free[nb]]
      if (length(nb) == 0) next

      # 轻微偏好“向当前已占质心靠近”，但仍保留随机
      if (any(!is.na(claim))) {
        cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
      } else {
        cx <- 0; cy <- 0
      }
      sc  <- -((layout$x[nb]-cx)^2 + (layout$y[nb]-cy)^2)
      # ord <- order(sc, runif(length(nb)))
      ord <- order(sc, nb)
      nb  <- nb[ord]
      q   <- c(q, nb)
    }
    region
  }

  # 若 BFS 不足 quota，则从全局 free 里按到 region/seed 最近补齐
  fill_deficit <- function(region, quota, target_idx){
    need <- quota - length(region)
    if (need <= 0) return(region)
    free_idx <- which(free)
    if (length(free_idx) == 0) return(region)

    if (length(region) > 0) {
      tx <- mean(layout$x[region]); ty <- mean(layout$y[region])
    } else {
      tx <- layout$x[target_idx];    ty <- layout$y[target_idx]
    }
    dd <- (layout$x[free_idx] - tx)^2 + (layout$y[free_idx] - ty)^2
    pick <- head(free_idx[order(dd)], min(need, length(free_idx)))
    free[pick] <<- FALSE
    c(region, pick)
  }

  # 逐模块（Others 留到最后）
  for (mi in seq_len(nrow(need_per_mod))) {
    mod    <- as.character(need_per_mod$modularity3[mi])
    k_need <- need_per_mod$n[mi]
    if (mod == "Others") next

    if (sum(free) < k_need) {
      stop(sprintf("剩余 layout 点不足：模块 %s 需要 %d，但只剩 %d。",
                   mod, k_need, sum(free)))
    }

    seed_idx <- if (mi == 1) pick_seed_first() else pick_seed_next()
    region   <- grow_region(seed_idx, k_need)
    if (length(region) < k_need) {
      region <- fill_deficit(region, k_need, target_idx = seed_idx)
    }

    claim[region] <- mod
    coords <- layout[region, c("x","y"), drop = FALSE]
    coords <- shrink_rings_global(coords, shrink = shrink)

    # 模块内：Degree降序 ↔ 半径升序（hub居中）
    cx <- mean(coords$x); cy <- mean(coords$y)
    rr <- sqrt((coords$x - cx)^2 + (coords$y - cy)^2)
    ord_pts <- order(rr, coords$x, coords$y)

    nodes_this <- nodes_sorted %>%
      dplyr::filter(modularity3 == mod)

    placed_coords[[mod]] <- dplyr::bind_cols(coords[ord_pts, , drop=FALSE], nodes_this)
  }

  ## ---------- 4) Others：外缘分配（里大外小；不改坐标，只改映射） ----------
  if ("Others" %in% names(placed_coords)) {
    n_oth <- need_per_mod$n[need_per_mod$modularity3 == "Others"]
    if (n_oth > 0) {
      idx_free <- which(free)
      if (length(idx_free) < n_oth) {
        stop(sprintf("Others 需要 %d 个点，但仅剩 %d 个。", n_oth, length(idx_free)))
      }
      # 仍然选“最外圈”的 n_oth 个点，并可整体外移
      r_free <- sqrt(layout$x[idx_free]^2 + layout$y[idx_free]^2)
      take   <- idx_free[order(-r_free)][seq_len(n_oth)]
      free[take] <- FALSE
      coords <- layout[take, c("x","y"), drop = FALSE]
      coords <- radial_offset(coords, delta = push_others_delta)  # 可为 0

      # 关键：只重排“节点→坐标”的对应关系
      # - 节点按 Degree 降序（高在前）
      nodes_oth <- nodes_sorted %>%
        dplyr::filter(modularity3 == "Others") %>%
        dplyr::arrange(dplyr::desc(Degree))

      # - 坐标按半径升序（越靠内的坐标放在前）
      rr   <- sqrt(coords$x^2 + coords$y^2)
      ordc <- order(rr, coords$x, coords$y)
      coords <- coords[ordc, , drop = FALSE]

      # 绑定：高 Degree ↔ 更靠内的坐标
      placed_coords[["Others"]] <- dplyr::bind_cols(coords, nodes_oth)
    }
  }

  ## ---------- 5) 拼装返回 ----------
  graph_ly_final <- dplyr::bind_rows(placed_coords[as.character(need_per_mod$modularity3)])
  ly_final <- graph_ly_final[, c("x","y")]

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

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)


  list(
    layout = ly_final,
    graph_obj = graph_obj_sort,
    graph_ly_final = graph_ly_final,
    ggplot_data = ggplot_data
  )
}
