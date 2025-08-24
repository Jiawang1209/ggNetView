create_layout_diamond <- function(graph_obj, node_add = 7, r = 0.1){
  # test
  graph_obj = graph_obj
  node_add = 7
  r = 0.1
  # 获取节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  # 节点个数
  n <- dim(node_df)[1]
  
  # # 每圈点数分配：与原先 circle_layout 一致
  # ring_counts <- function(n, node_add){
  #   counts <- 1
  #   total  <- 1
  #   i <- 2
  #   while (total < n) {
  #     add <- node_add * (i - 1)
  #     if (total + add <= n) {
  #       counts <- c(counts, add)
  #       total  <- total + add
  #     } else {
  #       counts <- c(counts, n - total)
  #       total  <- n
  #     }
  #     i <- i + 1
  #   }
  #   return(counts)
  # }
  
  
  # 获取到layout的具体信息 每一层，有多少个点
  # layout_df_info <- data.frame(
  #   number_circle = seq_along(ring_counts(n, node_add)),
  #   number = ring_counts(n, node_add)
  # )
  
  # 这里有一个立即调用函数的版本， 不需要再次调用函数，直接可以运行
  ring_counts <- (function(n, node_add){
    counts <- 1        # 初始中心点 1 个
    total  <- 1        # 已放置的总数
    i <- 2             # 从第 2 圈开始计算
    while (total < n) {
      add <- node_add * (i - 1)   # 第 i 圈计划放多少
      if (total + add <= n) {
        counts <- c(counts, add)  # 放得下就追加
        total  <- total + add
      } else {
        counts <- c(counts, n - total) # 不够一整圈就把剩余全放进去
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
  
  layout_df_info
  
  # 将 [0,1) 上的参数 u 均匀映射到菱形四条边（每条边长度占 1/4）
  diamond_param_to_xy <- function(u, radius){
    # u 可以是向量
    u <- u %% 1
    seg <- floor(u * 4)            # 0,1,2,3 四段
    t   <- (u * 4) - seg           # 段内 0-1
    # 分段计算
    x <- numeric(length(u))
    y <- numeric(length(u))
    # 段0： (0,r) -> (r,0)
    idx0 <- seg == 0
    x[idx0] <-  radius * t[idx0]
    y[idx0] <-  radius * (1 - t[idx0])
    # 段1： (r,0) -> (0,-r)
    idx1 <- seg == 1
    x[idx1] <-  radius * (1 - t[idx1])
    y[idx1] <- -radius * t[idx1]
    # 段2： (0,-r) -> (-r,0)
    idx2 <- seg == 2
    x[idx2] <- -radius * t[idx2]
    y[idx2] <- -radius * (1 - t[idx2])
    # 段3： (-r,0) -> (0,r)
    idx3 <- seg == 3
    x[idx3] <- -radius * (1 - t[idx3])
    y[idx3] <-  radius * t[idx3]
    data.frame(x = x, y = y)
  }
  
  # 逐圈生成坐标
  ly <- data.frame(x = 0, y = 0)  # 第一圈：中心点
  # 用一个逐圈累加的偏移，打散竖直/水平重合
  offset_accum <- 0
  
  for (index in 2:nrow(layout_df_info)) {
    m <- layout_df_info$number[index]
    radius <- (index - 1) * r
    # 均匀等距的参数；再加一个小偏移避免重合
    # 这里的偏移策略：每圈沿周长错开 1/(2m)，并与上一圈累计，类似你在圆形里的做法
    offset_accum <- (offset_accum + 0.5 / m) %% 1
    u <- ( (0:(m-1)) / m + offset_accum ) %% 1
    coords <- diamond_param_to_xy(u, radius)
    ly <- dplyr::bind_rows(ly, coords)
  }
  
  return(ly)
}
