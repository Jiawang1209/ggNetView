# 扇区环切布局（支持模块间 gap；hub 聚到各自模块中心）
module_layout_sector_by_ring <- function(
    graph_obj,
    layout,
    order = c("size_desc","custom"),
    modules_order = NULL,
    width_mode = c("proportional","equal"),
    others_last = TRUE,
    theta0_deg = 0,
    # hub 布局相关
    hub_mode = c("module_centroid","global_inner","none"),
    hub_alpha = 0.6,      # 半径靠模块中心的权重
    hub_beta  = 0.4,      # 角度靠模块中轴的权重
    r_center_stat = c("median","mean"),
    # 新增：模块间缝隙（度）
    sector_gap_deg = 6,   # 每个扇面之间留多少度；0 = 无缝
    seed = NULL
){
  stopifnot(all(c("x","y") %in% names(layout)))
  if (!is.null(seed)) set.seed(seed)
  order <- match.arg(order); width_mode <- match.arg(width_mode)
  hub_mode <- match.arg(hub_mode); r_center_stat <- match.arg(r_center_stat)

  to_rad <- function(d) d/180*pi
  wrap_theta <- function(th){ th2 <- th %% (2*pi); th2[th2 < 0] <- th2[th2 < 0] + 2*pi; th2 }
  ang <- function(x,y) wrap_theta(atan2(y,x))
  rad <- function(x,y) sqrt(x^2 + y^2)

  # ---------- 节点与模块信息 ----------
  node_df <- graph_obj %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
  stopifnot("modularity3" %in% names(node_df))
  mod_sizes <- node_df %>% dplyr::count(modularity3, name = "n") %>% dplyr::arrange(dplyr::desc(n))
  mods <- mod_sizes$modularity3

  if (order == "custom") {
    stopifnot(!is.null(modules_order))
    miss <- setdiff(mods, modules_order)
    if (length(miss)) stop("modules_order 缺少模块: ", paste(miss, collapse = ", "))
    mods_order <- modules_order
  } else {
    mods_order <- if (others_last && "Others" %in% mods) c(setdiff(mods, "Others"), "Others") else mods
  }

  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mods_order, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  sizes_vec <- setNames(mod_sizes$n[match(mods_order, mod_sizes$modularity3)], mods_order)
  total_n   <- sum(sizes_vec)

  # ---------- 识别同心环 ----------
  L <- layout
  r_all <- rad(L$x, L$y); theta_all <- ang(L$x, L$y)
  tol <- max(1e-9, stats::median(diff(sort(unique(r_all))))/100, na.rm = TRUE)
  r_levels <- sort(unique(round(r_all / tol) * tol))
  ring_id  <- match(round(r_all / tol) * tol, r_levels)
  rings    <- split(seq_len(nrow(L)), ring_id)
  rings    <- lapply(rings, function(idx){ idx[order(theta_all[idx], r_all[idx])] })
  cap      <- vapply(rings, length, integer(1))
  num_rings <- length(rings)

  # ---------- 环级配额（Hamilton + 剩余法） ----------
  remaining <- sizes_vec
  per_ring_take <- vector("list", num_rings)
  for (j in seq_len(num_rings)) {
    if (sum(remaining) == 0L) { per_ring_take[[j]] <- setNames(rep(0L, length(mods_order)), mods_order); next }
    Mj <- cap[j]
    share <- remaining / sum(remaining)
    quota <- share * Mj
    base  <- floor(quota)
    base[base > remaining] <- remaining[base > remaining]
    rest <- Mj - sum(base)
    if (rest > 0) {
      frac <- quota - floor(quota); frac[base >= remaining] <- -Inf
      add_order <- order(frac, decreasing = TRUE); i <- 1
      while (rest > 0 && i <= length(add_order)) {
        m <- add_order[i]
        if (base[m] < remaining[m]) { base[m] <- base[m] + 1L; rest <- rest - 1L }
        i <- i + 1L
      }
    }
    per_ring_take[[j]] <- setNames(as.integer(base), mods_order)
    remaining <- remaining - base
  }
  if (any(remaining != 0L)) stop("环级分配后全局剩余未清零，请检查。")

  # ---------- 环内切扇区（统一起点） ----------
  theta0 <- wrap_theta(to_rad(theta0_deg))
  ring_theta <- lapply(rings, function(idx) theta_all[idx])
  start_idx  <- vapply(ring_theta, function(th){ which.min(abs(atan2(sin(th - theta0), cos(th - theta0)))) }, integer(1))
  rings_rot  <- mapply(function(idx, s){ if (length(idx)==0) idx else c(idx[s:length(idx)], idx[1:(s-1)]) },
                       rings, start_idx, SIMPLIFY = FALSE)

  ring_assign <- vector("list", num_rings)
  for (j in seq_len(num_rings)) {
    idx_order <- rings_rot[[j]]
    takes <- per_ring_take[[j]]
    assign_j <- list(); cursor <- 1L
    for (m in mods_order) {
      k <- takes[[m]]
      assign_j[[m]] <- if (k > 0) idx_order[cursor:(cursor + k - 1L)] else integer(0)
      cursor <- cursor + k
    }
    ring_assign[[j]] <- assign_j
  }

  # ---------- 模块中轴角（全局一致，用于 gap 与 hub） ----------
  module_mid_angle <- setNames(numeric(length(mods_order)), mods_order)
  for (m in mods_order) {
    pos_m <- unlist(lapply(ring_assign, `[[`, m), use.names = FALSE)
    module_mid_angle[m] <- if (length(pos_m)) {
      atan2(mean(sin(theta_all[pos_m])), mean(cos(theta_all[pos_m])))
    } else NA_real_
  }

  # ---------- 节点→位置（hub 聚到各自模块中心） ----------
  g_nodes <- graph_obj_sort %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
  node_idx_by_mod <- split(seq_len(nrow(g_nodes)), g_nodes$modularity3)
  layout_row_for_node <- integer(nrow(g_nodes))

  wsum <- hub_alpha + hub_beta
  if (wsum <= 0) { hub_alpha <- 1; hub_beta <- 0 } else { hub_alpha <- hub_alpha/wsum; hub_beta <- hub_beta/wsum }

  for (m in mods_order) {
    nodes_m <- node_idx_by_mod[[m]]; if (length(nodes_m) == 0) next
    pos_m   <- unlist(lapply(ring_assign, `[[`, m), use.names = FALSE)
    if (length(pos_m) != length(nodes_m))
      stop(sprintf("模块 %s 的位置数(%d)与节点数(%d)不一致。", m, length(pos_m), length(nodes_m)))

    r_m  <- r_all[pos_m]; th_m <- theta_all[pos_m]
    th0 <- if (is.finite(module_mid_angle[m])) module_mid_angle[m] else 0
    r0  <- if (r_center_stat == "median") stats::median(r_m) else mean(r_m)

    r_rng <- max(1e-12, max(r_m) - min(r_m))
    ang_diff <- abs(atan2(sin(th_m - th0), cos(th_m - th0))) / pi   # 0..1
    rad_diff <- abs(r_m - r0) / r_rng                               # 0..1
    score_pos <- hub_alpha * rad_diff + hub_beta * ang_diff
    pos_order <- order(score_pos, rad_diff, ang_diff)

    if (hub_mode == "none") {
      node_order <- order(nodes_m)
    } else {
      if ("degree" %in% names(g_nodes)) node_order <- order(-g_nodes$degree[nodes_m], nodes_m) else node_order <- order(nodes_m)
    }

    layout_row_for_node[ nodes_m[ node_order ] ] <- pos_m[ pos_order ]
  }

  # ---------- 应用扇区 gap：按模块中轴线压缩角度 ----------
  K <- length(mods_order)
  gap_rad <- to_rad(sector_gap_deg)
  if (gap_rad < 0) stop("sector_gap_deg 不能为负。")
  total_gap <- K * gap_rad
  if (total_gap >= 2*pi - 1e-9)
    stop(sprintf("sector_gap_deg 过大：K*gap >= 360°。目前 K=%d, gap=%.2f°", K, sector_gap_deg))
  s_comp <- 1 - total_gap/(2*pi)   # (0,1] 压缩比例

  new_theta_pos <- theta_all
  for (j in seq_len(num_rings)) {
    for (m in mods_order) {
      pos <- ring_assign[[j]][[m]]; if (length(pos) == 0) next
      th  <- theta_all[pos]
      mid <- module_mid_angle[m]
      if (!is.finite(mid)) next
      d   <- atan2(sin(th - mid), cos(th - mid))  # 最短角差
      new_theta_pos[pos] <- mid + s_comp * d      # 向中轴压缩 => 留出模块间 gap
    }
  }

  x_pos <- r_all * cos(new_theta_pos)
  y_pos <- r_all * sin(new_theta_pos)
  ly_final <- data.frame(
    x = x_pos[ layout_row_for_node ],
    y = y_pos[ layout_row_for_node ]
  )

  list(layout = ly_final, graph_obj = graph_obj_sort)
}

res <- module_layout_sector_by_ring(
  graph_obj = graph_obj,
  layout    = ly,
  theta0_deg = 0,
  sector_gap_deg = 8,          # ← 模块之间的缝隙；常见 4~10°
  hub_mode   = "module_centroid",
  hub_alpha  = 0.65, hub_beta = 0.35,
  r_center_stat = "median",
  seed = 1115
)

####----test2----####
#' 同心圆扇区布局（融合“重要度排序 / 小模块合并 / 稀疏对齐 / 模块gap / 分层压缩 / hub就地”）
#' @param graph_obj   tidygraph::tbl_graph（节点需有 modularity3；若无 degree 会自动计算度）
#' @param layout      data.frame(x,y) 同心圆坐标池（行数=节点数）
#' @param order       c("auto","custom")；auto: 按 sum(degree) 再按 n；custom: 用 modules_order
#' @param modules_order 手动模块顺序（含 Others），当 order="custom" 时生效
#' @param others_last  TRUE/FALSE，是否把 Others 固定到最后
#' @param min_size     小模块合并阈值（<阈值并到 Others；0 关闭），默认 0
#' @param theta0_deg   起始参考角（度，0=正上方，顺时针）
#' @param auto_theta0  TRUE 时自动把 theta0 对齐到“全局最稀疏角域”的中心（更均衡）
#' @param sector_gap_deg 模块之间的缝隙角度（度；全局等宽）
#' @param hub_mode     c("module_centroid","global_inner","none")
#' @param hub_alpha    位置评分中“半径靠模块中心”的权重
#' @param hub_beta     位置评分中“角度靠中轴”的权重（两者会归一到和=1）
#' @param r_center_stat c("median","mean") 模块中心半径统计量
#' @param layer_comp_inner 内环角度压缩系数（0~1，越小越紧）
#' @param layer_comp_outer 外环角度压缩系数（0~1，外圈通常更松）
#' @param jitter_deg   环内同扇段轻微角度抖动（度；0 关闭）
#' @param seed         随机种子（用于并列打破/抖动）
#' @return list(layout = data.frame(x,y), graph_obj = 排序后的图对象)
module_layout_sector_hybrid <- function(
    graph_obj,
    layout,
    order = c("auto","custom"),
    modules_order = NULL,
    others_last = TRUE,
    min_size = 0,
    theta0_deg = 0,
    auto_theta0 = TRUE,
    sector_gap_deg = 6,
    hub_mode = c("module_centroid","global_inner","none"),
    hub_alpha = 0.6,
    hub_beta  = 0.4,
    r_center_stat = c("median","mean"),
    layer_comp_inner = 0.9,
    layer_comp_outer = 0.75,
    jitter_deg = 0,
    seed = NULL
){
  stopifnot(all(c("x","y") %in% names(layout)))
  if (!is.null(seed)) set.seed(seed)
  order <- match.arg(order)
  hub_mode <- match.arg(hub_mode)
  r_center_stat <- match.arg(r_center_stat)

  to_rad <- function(d) d/180*pi
  wrap_theta <- function(th){ th2 <- th %% (2*pi); th2[th2 < 0] <- th2[th2 < 0] + 2*pi; th2 }
  ang <- function(x,y) wrap_theta(atan2(y,x))
  rad <- function(x,y) sqrt(x^2 + y^2)

  # ---------- 准备节点&度 ----------
  g0 <- graph_obj
  # 若没有 degree 就算一个（无权度）
  if (!"degree" %in% names(tidygraph::as_tibble(tidygraph::activate(g0, "nodes")))) {
    g0 <- g0 %>% tidygraph::mutate(degree = tidygraph::centrality_degree(mode = "all"))
  }

  node_df <- g0 %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
  stopifnot("modularity3" %in% names(node_df))

  # 小模块合并到 Others
  if (min_size > 0) {
    tabn <- dplyr::count(node_df, modularity3, name = "n")
    small <- tabn$modularity3[ tabn$n < min_size & tabn$modularity3 != "Others" ]
    if (length(small)) {
      g0 <- g0 %>% tidygraph::mutate(
        modularity3 = ifelse(modularity3 %in% small, "Others", modularity3)
      )
      node_df <- g0 %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
    }
  }

  # 模块统计：sum_degree 优先，再 n
  mods_stat <- node_df %>%
    dplyr::group_by(modularity3) %>%
    dplyr::summarise(
      n = dplyr::n(),
      sum_degree = sum(.data$degree, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_degree), dplyr::desc(n))

  mods_all <- mods_stat$modularity3
  if (order == "custom") {
    stopifnot(!is.null(modules_order))
    miss <- setdiff(mods_all, modules_order)
    if (length(miss)) stop("modules_order 缺少模块: ", paste(miss, collapse = ", "))
    mods_order <- modules_order
  } else {
    mods_order <- mods_stat$modularity3
  }
  if (others_last && "Others" %in% mods_order) {
    mods_order <- c(setdiff(mods_order, "Others"), "Others")
  }

  # 规范后的图对象顺序
  graph_obj_sort <- g0 %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mods_order, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  # 模块规模向量
  sizes_vec <- setNames(mods_stat$n[match(mods_order, mods_stat$modularity3)], mods_order)
  total_n   <- sum(sizes_vec)

  # ---------- 同心环解析 ----------
  L <- layout
  r_all <- rad(L$x, L$y); theta_all <- ang(L$x, L$y)
  tol <- max(1e-9, stats::median(diff(sort(unique(r_all))))/100, na.rm = TRUE)
  r_levels <- sort(unique(round(r_all / tol) * tol))
  ring_id  <- match(round(r_all / tol) * tol, r_levels)
  rings    <- split(seq_len(nrow(L)), ring_id)
  rings    <- lapply(rings, function(idx){ idx[order(theta_all[idx], r_all[idx])] })
  cap      <- vapply(rings, length, integer(1))
  num_rings <- length(rings)

  # ---------- 自动选择稀疏起点（可选） ----------
  theta0 <- to_rad(theta0_deg)
  if (auto_theta0) {
    # 角度直方图：找到最稀疏的 bin 中心作为 theta0
    bins <- 180
    brks <- seq(0, 2*pi, length.out = bins+1)
    h <- hist(theta_all, breaks = brks, plot = FALSE)$counts
    k <- which.min(h); mid <- (brks[k] + brks[k+1])/2
    theta0 <- mid
  }
  theta0 <- wrap_theta(theta0)

  # ---------- 环级配额（Hamilton + 剩余） ----------
  remaining <- sizes_vec
  per_ring_take <- vector("list", num_rings)
  for (j in seq_len(num_rings)) {
    if (sum(remaining) == 0L) { per_ring_take[[j]] <- setNames(rep(0L, length(mods_order)), mods_order); next }
    Mj <- cap[j]
    share <- remaining / sum(remaining)
    quota <- share * Mj
    base  <- floor(quota)
    base[base > remaining] <- remaining[base > remaining]
    rest <- Mj - sum(base)
    if (rest > 0) {
      frac <- quota - floor(quota); frac[base >= remaining] <- -Inf
      add_order <- order(frac, decreasing = TRUE); i <- 1
      while (rest > 0 && i <= length(add_order)) {
        m <- add_order[i]
        if (base[m] < remaining[m]) { base[m] <- base[m] + 1L; rest <- rest - 1L }
        i <- i + 1L
      }
    }
    per_ring_take[[j]] <- setNames(as.integer(base), mods_order)
    remaining <- remaining - base
  }
  if (any(remaining != 0L)) stop("环级分配后全局剩余未清零，请检查。")

  # ---------- 统一起点并切扇段 ----------
  ring_theta <- lapply(rings, function(idx) theta_all[idx])
  start_idx  <- vapply(ring_theta, function(th){ which.min(abs(atan2(sin(th - theta0), cos(th - theta0)))) }, integer(1))
  rings_rot  <- mapply(function(idx, s){ if (length(idx)==0) idx else c(idx[s:length(idx)], idx[1:(s-1)]) },
                       rings, start_idx, SIMPLIFY = FALSE)

  ring_assign <- vector("list", num_rings)
  for (j in seq_len(num_rings)) {
    idx_order <- rings_rot[[j]]
    takes <- per_ring_take[[j]]
    assign_j <- list(); cursor <- 1L
    for (m in mods_order) {
      k <- takes[[m]]
      assign_j[[m]] <- if (k > 0) idx_order[cursor:(cursor + k - 1L)] else integer(0)
      cursor <- cursor + k
    }
    ring_assign[[j]] <- assign_j
  }

  # ---------- 模块中轴角（供 gap/定位/分层压缩） ----------
  module_mid_angle <- setNames(numeric(length(mods_order)), mods_order)
  for (m in mods_order) {
    pos_m <- unlist(lapply(ring_assign, `[[`, m), use.names = FALSE)
    module_mid_angle[m] <- if (length(pos_m)) {
      atan2(mean(sin(theta_all[pos_m])), mean(cos(theta_all[pos_m])))
    } else NA_real_
  }

  # ---------- 节点 → 位置（让 hub 落在模块中心带） ----------
  g_nodes <- graph_obj_sort %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
  node_idx_by_mod <- split(seq_len(nrow(g_nodes)), g_nodes$modularity3)
  layout_row_for_node <- integer(nrow(g_nodes))

  wsum <- hub_alpha + hub_beta
  if (wsum <= 0) { hub_alpha <- 1; hub_beta <- 0 } else { hub_alpha <- hub_alpha/wsum; hub_beta <- hub_beta/wsum }

  for (m in mods_order) {
    nodes_m <- node_idx_by_mod[[m]]; if (length(nodes_m) == 0) next
    pos_m   <- unlist(lapply(ring_assign, `[[`, m), use.names = FALSE)
    if (length(pos_m) != length(nodes_m))
      stop(sprintf("模块 %s 的位置数(%d)与节点数(%d)不一致。", m, length(pos_m), length(nodes_m)))

    r_m  <- r_all[pos_m]; th_m <- theta_all[pos_m]
    th0 <- if (is.finite(module_mid_angle[m])) module_mid_angle[m] else 0
    r0  <- if (r_center_stat == "median") stats::median(r_m) else mean(r_m)

    r_rng <- max(1e-12, max(r_m) - min(r_m))
    ang_diff <- abs(atan2(sin(th_m - th0), cos(th_m - th0))) / pi   # 0..1
    rad_diff <- abs(r_m - r0) / r_rng                               # 0..1
    score_pos <- hub_alpha * rad_diff + hub_beta * ang_diff
    pos_order <- order(score_pos, rad_diff, ang_diff)

    if (hub_mode == "none") {
      node_order <- order(nodes_m)
    } else if (hub_mode == "global_inner") {
      if ("degree" %in% names(g_nodes)) node_order <- order(-g_nodes$degree[nodes_m], nodes_m) else node_order <- order(nodes_m)
    } else { # module_centroid
      if ("degree" %in% names(g_nodes)) node_order <- order(-g_nodes$degree[nodes_m], nodes_m) else node_order <- order(nodes_m)
    }

    layout_row_for_node[ nodes_m[ node_order ] ] <- pos_m[ pos_order ]
  }

  # ---------- 应用模块 gap（角度向中轴压缩，留出缝隙） ----------
  K <- length(mods_order)
  gap_rad <- to_rad(sector_gap_deg)
  if (gap_rad < 0) stop("sector_gap_deg 不能为负。")
  total_gap <- K * gap_rad
  if (total_gap >= 2*pi - 1e-9)
    stop(sprintf("sector_gap_deg 过大：K*gap >= 360°。目前 K=%d, gap=%.2f°", K, sector_gap_deg))
  s_comp_gap <- 1 - total_gap/(2*pi)   # 全局压缩比例 (0,1]

  theta_gap <- theta_all
  for (j in seq_len(num_rings)) {
    for (m in mods_order) {
      pos <- ring_assign[[j]][[m]]; if (length(pos) == 0) next
      th  <- theta_all[pos]
      mid <- module_mid_angle[m]
      if (!is.finite(mid)) next
      d   <- atan2(sin(th - mid), cos(th - mid))
      theta_gap[pos] <- mid + s_comp_gap * d
    }
  }

  # ---------- 分层角度压缩（内紧外松的“花瓣感”） ----------
  comp_inner <- layer_comp_inner; comp_outer <- layer_comp_outer
  comp_inner <- max(min(comp_inner, 1), 0); comp_outer <- max(min(comp_outer, 1), 0)
  if (num_rings > 1 && (abs(comp_inner - comp_outer) > 1e-9)) {
    theta_layer <- theta_gap
    for (j in seq_len(num_rings)) {
      ring_pos <- rings[[j]]; if (!length(ring_pos)) next
      # 线性插值：内环用 comp_inner，外环用 comp_outer
      t <- if (num_rings == 1) 0 else (j-1)/(num_rings-1)
      s_comp_layer <- comp_inner + t * (comp_outer - comp_inner)
      for (m in mods_order) {
        pos <- intersect(ring_pos, unlist(lapply(list(ring_assign[[j]][[m]]), identity), use.names = FALSE))
        if (!length(pos)) next
        th  <- theta_gap[pos]
        mid <- module_mid_angle[m]; if (!is.finite(mid)) next
        d   <- atan2(sin(th - mid), cos(th - mid))
        theta_layer[pos] <- mid + s_comp_layer * d
      }
    }
    theta_fin <- theta_layer
  } else {
    theta_fin <- theta_gap
  }

  # ---------- 轻微角度抖动（可选，缓解重合） ----------
  if (jitter_deg > 0) {
    jit <- rnorm(length(theta_fin), mean = 0, sd = to_rad(jitter_deg))
    theta_fin <- wrap_theta(theta_fin + jit)
  }

  # ---------- 输出最终布局（与 graph_obj_sort 顺序一一对应） ----------
  x_new <- r_all * cos(theta_fin); y_new <- r_all * sin(theta_fin)
  ly_final <- data.frame(
    x = x_new[ layout_row_for_node ],
    y = y_new[ layout_row_for_node ]
  )

  list(layout = ly_final, graph_obj = graph_obj_sort)
}

res <- module_layout_sector_hybrid(
  graph_obj = graph_obj,
  layout    = ly,
  # 顺序&合并
  order = "auto", others_last = TRUE, min_size = 60,
  # 稀疏角起点 + 扇区缝隙
  theta0_deg = 0, auto_theta0 = TRUE, sector_gap_deg = 8,
  # hub落位：模块中心
  hub_mode = "module_centroid", hub_alpha = 0.65, hub_beta = 0.35, r_center_stat = "median",
  # 分层压缩（内紧外松）+ 轻抖动
  layer_comp_inner = 0.9, layer_comp_outer = 0.75, jitter_deg = 0.5,
  seed = 1115
)




maskTable <- mascarade::generateMask(dims= res[["layout"]],
                                     clusters=res[["graph_obj"]] %>%
                                       tidygraph::activate(nodes) %>%
                                       tidygraph::as_tibble() %>%
                                       dplyr::pull(modularity3))

ggraph(res$graph_obj, layout = "manual", x = res$layout$x, y = res$layout$y) +
  geom_edge_link(alpha = 0.15, colour = "grey70") +
  geom_node_point(aes(fill = modularity3, size = degree),
                  shape = 21, stroke = 0.2, alpha = 0.9) +
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
                        linewidth = 1.25,
                        linetype = 2,
                        alpha = 0.5,
                        show.legend = F) +
  ggplot2::scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                        '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                        '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                        '#bdbdbd',
                                        '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                        '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                        '#ffff99','#b15928'),
                             name = "modularity") +
  ggplot2::scale_color_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                        '#fdb462','#b3de69','#fccde5','#cab2d6','#bc80bd',
                                        '#ccebc5','#ffed6f','#a6cee3','#b2df8a', '#fb9a99',
                                        '#bdbdbd',
                                        '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                                        '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a',
                                        '#ffff99','#b15928'),
                             name = "modularity") +
  theme_void()
