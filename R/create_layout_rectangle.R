create_layout_rectangle <- function(graph_obj, node_add = 7, r = 0.1, ratio = 1.5){
  # graph 的节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  
  # 每圈放多少个点（与原先一致：第1圈1个，其后第 i 圈增 node_add*(i-1)）
  ring_counts <- (function(n, node_add){
    counts <- 1
    total  <- 1
    i <- 2
    while (total < n) {
      add <- node_add * (i - 1)
      if (total + add <= n) {
        counts <- c(counts, add)
        total  <- total + add
      } else {
        counts <- c(counts, n - total)
        total  <- n
      }
      i <- i + 1
    }
    counts
  })(n, node_add)
  
  layout_df_info <- data.frame(
    number_circle = seq_along(ring_counts),
    number = ring_counts
  )
  
  # 把 u∈[0,1) 等距映射到矩形四条边（按周长比例分配）
  # a = 半宽, b = 半高；四边长度：上/下各 2a，左/右各 2b
  rectangle_param_to_xy <- function(u, a, b){
    u <- u %% 1
    
    fa <- a / (a + b)           # 每条水平边占的周长比例 (2a / 4(a+b))
    fb <- b / (a + b)           # 每条垂直边占的周长比例 (2b / 4(a+b))
    
    s0 <- fa            # 顶边区间 [0, s0)
    s1 <- fa + fb       # 右边区间 [s0, s1)
    s2 <- fa + fb + fa  # 底边区间 [s1, s2)
    # 左边区间 [s2, 1)
    
    x <- numeric(length(u))
    y <- numeric(length(u))
    
    # 顶边：(-a, b) -> (a, b)
    idx0 <- u < s0
    if (any(idx0)) {
      t0 <- u[idx0] / fa
      x[idx0] <- -a + 2 * a * t0
      y[idx0] <-  b
    }
    
    # 右边：(a, b) -> (a, -b)
    idx1 <- (u >= s0) & (u < s1)
    if (any(idx1)) {
      t1 <- (u[idx1] - s0) / fb
      x[idx1] <-  a
      y[idx1] <-  b - 2 * b * t1
    }
    
    # 底边：(a, -b) -> (-a, -b)
    idx2 <- (u >= s1) & (u < s2)
    if (any(idx2)) {
      t2 <- (u[idx2] - s1) / fa
      x[idx2] <-  a - 2 * a * t2
      y[idx2] <- -b
    }
    
    # 左边：(-a, -b) -> (-a, b)
    idx3 <- (u >= s2)
    if (any(idx3)) {
      t3 <- (u[idx3] - s2) / fb
      x[idx3] <- -a
      y[idx3] <- -b + 2 * b * t3
    }
    
    data.frame(x = x, y = y)
  }
  
  # 逐圈生成坐标
  ly <- data.frame(x = 0, y = 0)    # 第1圈中心点
  offset_accum <- 0                 # 沿周长方向的错位，避免重叠
  
  for (index in 2:nrow(layout_df_info)) {
    m <- layout_df_info$number[index]
    # 当前圈的半宽/半高：让宽 = ratio * 高
    half_height <- (index - 1) * r
    half_width  <- ratio * half_height
    
    # 累积错位：每圈平移 1/(2m) 的周长比例
    offset_accum <- (offset_accum + 0.5 / m) %% 1
    u <- ((0:(m - 1)) / m + offset_accum) %% 1
    
    coords <- rectangle_param_to_xy(u, a = half_width, b = half_height)
    ly <- dplyr::bind_rows(ly, coords)
  }
  
  ly
}
