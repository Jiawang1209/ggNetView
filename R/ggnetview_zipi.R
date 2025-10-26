ggnetview_zipi <- function(nodes_bulk, z_bulk_mat, modularity_col, degree_col) {
  stopifnot(is.matrix(z_bulk_mat) || is.data.frame(z_bulk_mat))
  z_bulk_mat <- as.matrix(z_bulk_mat)

  # —— 对齐：行名必须覆盖 nodes_bulk 的行名
  if (!all(rownames(nodes_bulk) %in% rownames(z_bulk_mat))) {
    stop("rownames(nodes_bulk) 必须是 z_bulk_mat 的子集并对齐。")
  }
  # 重排矩阵顺序以匹配 nodes_bulk
  z_bulk_mat <- z_bulk_mat[rownames(nodes_bulk), rownames(nodes_bulk), drop = FALSE]

  # —— 二值化 & 处理对角线（可选：把对角线先置 0，再不做 -1）
  A <- (abs(z_bulk_mat) > 0) * 1L
  diag(A) <- 1L  # 确保自连为 1，便于后续 -1

  mod  <- nodes_bulk[[modularity_col]]
  deg  <- nodes_bulk[[degree_col]]
  ids  <- rownames(nodes_bulk)

  # 安全性
  if (any(is.na(mod))) stop("模块列存在 NA。")
  if (any(is.na(deg))) stop("度列存在 NA。")

  # —— 计算 within-module degree z
  # 按模块拆分索引
  split_idx <- split(seq_along(ids), f = factor(mod, levels = unique(mod)))
  z_vec <- numeric(length(ids)); names(z_vec) <- ids

  for (lev in names(split_idx)) {
    idx <- split_idx[[lev]]
    if (length(idx) <= 1) {
      z_vec[idx] <- 0
      next
    }
    Aii <- A[idx, idx, drop = FALSE]
    k_in <- rowSums(Aii) - 1L
    sd_k <- stats::sd(k_in)
    if (sd_k == 0) z_vec[idx] <- 0 else z_vec[idx] <- (k_in - mean(k_in)) / sd_k
  }

  # —— 计算参与系数 P
  # k_is：每个节点对每个模块的边数
  modules <- names(split_idx)
  kis_mat <- sapply(modules, function(lev) {
    idx <- split_idx[[lev]]
    rowSums(A[, idx, drop = FALSE])
  })
  # 去掉自身对角线：属于该模块的节点，k_is 减 1
  for (j in seq_along(modules)) {
    idx <- split_idx[[modules[j]]]
    kis_mat[idx, j] <- kis_mat[idx, j] - 1L
  }
  kis_mat[kis_mat < 0] <- 0  # 理论上不会小于0，保险

  sum_kis2 <- rowSums(kis_mat^2)
  k_tot    <- as.numeric(deg)
  P <- numeric(length(k_tot))
  P[k_tot == 0] <- 0
  nz <- (k_tot > 0)
  P[nz] <- 1 - (sum_kis2[nz] / (k_tot[nz]^2))
  names(P) <- ids

  # —— 组织输出
  out <- data.frame(
    nodes_id = ids,
    within_module_connectivities = z_vec[ids],
    among_module_connectivities  = P[ids],
    row.names = NULL,
    check.names = FALSE
  )
  nodes_bulk$nodes_id <- ids
  merge(out, nodes_bulk, by = "nodes_id")
}
