# 内圈同心圆 + 外圈花瓣（可选平滑过渡），严格对称、无偏移
create_layout_petal2 <- function(graph_obj, 
                                 node_add = 7, r = 0.1,
                                 petals = 6, amp = 0.35,
                                 inner_rings = 2,          # 有多少“内圈”保持纯圆形（不含中心点）
                                 transition_rings = 0      # 从圆形过渡到花瓣需要的圈数（0=立即切换）
){
  # 节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  
  # 每圈点数分配：第1圈1个，其后第 i 圈增 node_add*(i-1)
  ring_counts <- (function(n, node_add){
    counts <- 1; total <- 1; i <- 2
    while (total < n) {
      add <- node_add * (i - 1)
      if (total + add <= n) { counts <- c(counts, add); total <- total + add
      } else { counts <- c(counts, n - total); total <- n }
      i <- i + 1
    }
    counts
  })(n, node_add)
  
  layout_df_info <- data.frame(
    number_circle = seq_along(ring_counts),
    number = ring_counts
  )
  
  # 结果容器：中心点
  ly <- data.frame(x = 0, y = 0)
  
  # 从第二圈开始
  for (index in 2:nrow(layout_df_info)) {
    m <- layout_df_info$number[index]
    R_base <- (index - 1) * r               # 本圈基准半径
    theta  <- 2 * pi * (0:(m - 1)) / m      # 等角度采样（无偏移，保证对称）
    
    # 决定这一圈的“有效幅度”：内圈=0；外圈逐步从0过渡到 amp
    k <- index - 1                           # 逻辑上的圈号（1=第一圈外圈）
    if (k <= inner_rings) {
      amp_eff <- 0
    } else if (transition_rings <= 0) {
      amp_eff <- amp
    } else {
      t <- (k - inner_rings) / transition_rings   # 0→1
      t <- max(0, min(1, t))
      amp_eff <- amp * t
    }
    
    radius <- R_base * (1 + amp_eff * cos(petals * theta))
    
    # 生成坐标（内圈就是纯圆，外圈是花瓣）
    x <- radius * cos(theta)
    y <- radius * sin(theta)
    
    ly <- dplyr::bind_rows(ly, data.frame(x = x, y = y))
  }
  
  return(ly)
}
