# 等腰直角三角形同心布局（质心居中到原点；第一个点 = (0,0) = 总中心）
# 层配额：1, 3, 6, 9, ...（第 k 层放 3*(k-1) 个，均匀落在三条边的中点位置上）
create_layout_rightiso_layers <- function(
    graph_obj,
    r = 0.2,                                  # 层间尺度步长（直角边长度 L 的增量）
    orientation = c("up","down","left","right"),
    angle = 0                                 # 在 orientation 基础上的微调（弧度）
){
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle
  
  # 1) 节点总数
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  n <- nrow(node_df)
  if (n <= 0) return(data.frame(x=numeric(), y=numeric(), ring=integer(), idx=integer()))
  
  # 2) 层配额：1, 3, 6, 9, ... 截断到 n
  ring_sizes <- c(1L)
  total <- 1L; k <- 2L
  while (total < n) {
    add <- 3L * (k - 1L)
    take <- min(add, n - total)
    ring_sizes <- c(ring_sizes, take)
    total <- total + take
    k <- k + 1L
  }
  
  # 3) 给定直角边长 L，构造“质心在原点”的等腰直角三角形三个顶点
  #    先用标准顶点 A(0,0), B(-L/√2,-L/√2), C(L/√2,-L/√2)，其质心 G=(A+B+C)/3
  #    然后把三个顶点整体平移 -G，使质心到原点
  verts_centered <- function(L){
    if (L == 0) {
      A <- c(0,0); B <- c(0,0); C <- c(0,0)
      return(list(A=A,B=B,C=C, lengths=c(0,0,0), P=0))
    }
    A <- c(0, 0)
    B <- c(-L / sqrt(2), -L / sqrt(2))
    C <- c( L / sqrt(2), -L / sqrt(2))
    G <- (A + B + C) / 3  # 质心
    A <- A - G; B <- B - G; C <- C - G
    
    len_AB <- sqrt(sum((A-B)^2))  # = L
    len_BC <- sqrt(sum((B-C)^2))  # = sqrt(2)*L
    len_CA <- sqrt(sum((C-A)^2))  # = L
    P <- len_AB + len_BC + len_CA
    list(A=A,B=B,C=C, lengths=c(len_AB,len_BC,len_CA), P=P)
  }
  
  # 4) 周长参数 s -> (x,y) 映射（在“质心居中”的三角形上）
  map_s_to_xy_centered <- function(s, L){
    vc <- verts_centered(L)
    A <- vc$A; B <- vc$B; C <- vc$C
    len_AB <- vc$lengths[1]; len_BC <- vc$lengths[2]; len_CA <- vc$lengths[3]
    P <- vc$P
    if (P == 0) return(cbind(0,0))
    
    s <- s %% P
    xy <- matrix(NA_real_, nrow = length(s), ncol = 2)
    
    # 段1：A -> B（半开）
    idx1 <- s < len_AB
    if (any(idx1)) {
      t <- s[idx1] / len_AB
      xy[idx1,1] <- (1 - t) * A[1] + t * B[1]
      xy[idx1,2] <- (1 - t) * A[2] + t * B[2]
    }
    # 段2：B -> C（半开）
    idx2 <- (s >= len_AB) & (s < len_AB + len_BC)
    if (any(idx2)) {
      s2 <- s[idx2] - len_AB
      t  <- s2 / len_BC
      xy[idx2,1] <- B[1] + t * (C[1] - B[1])
      xy[idx2,2] <- B[2] + t * (C[2] - B[2])
    }
    # 段3：C -> A（半开）
    idx3 <- s >= (len_AB + len_BC)
    if (any(idx3)) {
      s3 <- s[idx3] - (len_AB + len_BC)
      t  <- s3 / len_CA
      xy[idx3,1] <- C[1] + t * (A[1] - C[1])
      xy[idx3,2] <- C[2] + t * (A[2] - C[2])
    }
    
    xy
  }
  
  out <- list()
  # 第1层：把“中心点”放在整体质心处——这里就是原点 (0,0)
  out[[1]] <- data.frame(x = 0, y = 0, ring = 1L, idx = 1L)
  
  # 5) 逐层：等分周长、取各段中点（避开顶点），无重复
  if (length(ring_sizes) >= 2) {
    for (ring in 2:length(ring_sizes)) {
      m <- ring_sizes[ring]
      L <- (ring - 1) * r
      if (m <= 0) next
      
      M_full <- 3L * (ring - 1L)          # 满额容量
      vc <- verts_centered(L)
      P  <- vc$P
      
      # 中点采样：((0…M_full-1)+0.5)/M_full * P
      s_mid <- ((0:(M_full - 1L)) + 0.5) / M_full * P
      s_use <- s_mid[seq_len(m)]
      
      xy <- map_s_to_xy_centered(s_use, L)
      
      # 旋转（先居中，后统一旋转，避免分段边界误判）
      if (theta_shift != 0) {
        Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                       sin(theta_shift),  cos(theta_shift)), nrow = 2)
        xy <- t(Rm %*% t(xy))
      }
      
      out[[ring]] <- data.frame(x = xy[,1], y = xy[,2],
                                ring = ring, idx = seq_len(m))
    }
  }
  
  layout_df <- dplyr::bind_rows(out) %>%
    dplyr::select(x, y)
  
  # （可选）数值去重校验
  # eps <- 1e-12
  # key <- paste(round(layout_df$x/eps), round(layout_df$y/eps))
  # stopifnot(!any(duplicated(key)))
  
  layout_df
}
