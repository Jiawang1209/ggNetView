# 同心五角星布局（无重复坐标；沿边等距中点采样）
create_layout_star <- function(
    graph_obj,
    node_add = 7,            # 每圈名义新增点数（与原逻辑一致）
    r = 0.1,                 # 外半径步长（相邻两圈外半径差）
    inner_ratio = 0.45,      # 内半径/外半径（0.35~0.55 常用，越小越尖）
    orientation = c("up","down","left","right"),
    angle = 0                # 在 orientation 基础上的微调（弧度）
){
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle
  
  # 节点数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  if (n <= 0) return(data.frame(x=numeric(), y=numeric()))
  
  # 与你原函数一致的配额：1, 7, 14, ... 直到放满 n
  ring_counts <- (function(n, node_add){
    counts <- 1L; total <- 1L; i <- 2L
    while (total < n) {
      add <- node_add * (i - 1L)
      if (total + add <= n) { counts <- c(counts, add); total <- total + add
      } else { counts <- c(counts, n - total); total <- n }
      i <- i + 1L
    }
    counts
  })(n, node_add)
  
  info <- data.frame(ring = seq_along(ring_counts), m = ring_counts)
  
  # —— 构造“朝上”的五角星 10 顶点（外内外内…交替），给定外半径 Ro 与内半径 Ri —— #
  star_vertices <- function(Ro, Ri){
    # 让第一个外顶点朝正上（pi/2）
    ang0 <- pi/2
    angs <- ang0 + (0:9) * (pi/5)  # 10 个等角
    rad  <- ifelse((0:9) %% 2 == 0, Ro, Ri)
    x <- rad * cos(angs)
    y <- rad * sin(angs)
    cbind(x, y)
  }
  
  # —— 按周长位置 s（0..P）插值到五角星边上（线性插值；半开分段避免重复） —— #
  map_s_to_star <- function(s, Ro, Ri){
    V  <- star_vertices(Ro, Ri)        # 10 x 2
    V2 <- V[c(2:10, 1), , drop = FALSE]# 邻边终点
    seglen <- sqrt(rowSums((V2 - V)^2))
    P <- sum(seglen)
    if (P == 0) return(cbind(0,0))
    
    # 归一化到 [0, P)
    s <- s %% P
    cum <- c(0, cumsum(seglen))        # 长度 11
    # 找到 s 落在哪一段：cum[j] <= s < cum[j+1]
    idx <- findInterval(s, cum, rightmost.closed = FALSE)  # 1..10
    idx[idx < 1]  <- 1L
    idx[idx > 10] <- 10L
    
    s0 <- cum[idx]
    s1 <- cum[idx + 1L]
    w  <- (s - s0) / pmax(s1 - s0, .Machine$double.eps)
    
    P0 <- V [idx, , drop = FALSE]
    P1 <- V2[idx, , drop = FALSE]
    xy <- P0 + (P1 - P0) * w
    xy
  }
  
  out <- list()
  # 第 1 圈：中心点
  out[[1]] <- data.frame(x = 0, y = 0, ring = 1L, idx = 1L)
  
  # 后续各圈：在五角星边上做“等分段中点”采样，避开顶点
  if (nrow(info) >= 2){
    for (row in 2:nrow(info)) {
      ring_i <- info$ring[row]
      m      <- info$m[row]
      Ro     <- (ring_i - 1) * r
      Ri     <- Ro * inner_ratio
      if (m <= 0 || Ro <= 0) next
      
      # 用该圈的“满额 m”把周长分成 m 段，取中点，避免采到 10 个顶点
      # （若最后一圈是截断，这里的 m 已经是截断后的）
      # 先算总周长
      Vtmp  <- star_vertices(Ro, Ri)
      seg   <- sqrt(rowSums((Vtmp[c(2:10,1),]-Vtmp)^2))
      P     <- sum(seg)
      s_mid <- ((0:(m-1)) + 0.5) / m * P
      
      xy <- map_s_to_star(s_mid, Ro, Ri)
      
      # 统一旋转
      if (theta_shift != 0) {
        Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                       sin(theta_shift),  cos(theta_shift)), nrow = 2)
        xy <- t(Rm %*% t(xy))
      }
      
      out[[ring_i]] <- data.frame(x = xy[,1], y = xy[,2],
                                  ring = ring_i, idx = seq_len(m))
    }
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::select(x, y)
}
