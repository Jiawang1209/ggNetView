create_layout_bipartite_gephi_layout <- function(
    graph_obj,
    r = 1,
    node_add = 7,
    scale = scale,
    orientation = c("up","down","left","right"),
    angle = 0
){

  # 既然，要画二分网络，那么你肯定要有二分网络
  # graph_obj = graph_sub[[3]]
  # 旋转角度
  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + angle

  # set radius
  radius = r * 10

  # 获取节点
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  Modularity_name <- node_df$Modularity  %>% droplevels() %>% levels() %>% as.character()

  mod_levels <- node_df$Modularity %>% droplevels() %>% levels() %>% as.character()

  module_list <- node_df %>% dplyr::group_split(Modularity)

  # 每一个分组的节点数
  n_vec <- purrr::map(module_list, ~dim(.x)[1]) %>% unlist()

  if (length(n_vec) < 2) {
    stop("Bipartite layout requires at least 2 modules in `Modularity`.")
  }
  if (length(n_vec) > 2) {
    message("`Bipartite layout more than 2 modules detected.")
  }

  # 我们统计一下每一个分组当中，应该有多少个点
  circle_layout <- function(n, node_add){
    # 初始化
    counts <- 1
    total <- 1
    i <- 2
    while (total < n) {
      add <- node_add * (i-1)
      if (total + add <= n) {
        counts <- c(counts, add)
        total <- total + add
      }else{
        # 最后一圈
        counts <- c(counts, n-total)
        total <- n
      }
      i <- i + 1
    }
    return(counts)
  }

  n_vec

  n_vec_node <- purrr::map(n_vec, function(x){
    data.frame(
      number_circle = seq_along(circle_layout(n = x, node_add = node_add)),
      number_node = circle_layout(n = x, node_add = node_add)
    )
  })

  ly_list <- list()

  for (i in seq_along(n_vec_node)) {
    if (i == 1) {
      ly <- data.frame(x = -1*radius, y = 0)
      offset <- 0
      prev_n <- n_vec_node[[i]]$number_node

      # 第二圈 到 最后一圈
      for (index in 2:(dim(n_vec_node[[i]])[1])) {
        if (index == 2) {
          # index = 2
          # l <- seq(0, 2*pi, length.out = prev_n[index])
          l <- 2* pi * (0:(prev_n[index]-1)) / prev_n[index]
        }else{
          # 第三圈开始 错开半个身位置
          offset <- pi/prev_n[index] %% (2*pi) + offset
          # l <- offset + seq(0, 2*pi, length.out = prev_n[index])
          l <- offset + (2* pi * (0:(prev_n[index]-1)) / prev_n[index])
        }

        x <- -1*radius + sin(l) * (index-1)*r
        y <- 0 + cos(l) * (index-1)*r
        ly_tmp <- data.frame(x = x,
                             y = y)

        ly <- dplyr::bind_rows(ly,ly_tmp)

        ly_list[[i]] <- ly

      }
    }

    if (i == 2) {
      ly <- data.frame(x = 1*radius, y = 0)
      offset <- 0
      prev_n <- n_vec_node[[i]]$number_node

      # 第二圈 到 最后一圈
      for (index in 2:(dim(n_vec_node[[i]])[1])) {
        if (index == 2) {
          # index = 2
          # l <- seq(0, 2*pi, length.out = prev_n[index])
          l <- 2* pi * (0:(prev_n[index]-1)) / prev_n[index]
        }else{
          # 第三圈开始 错开半个身位置
          offset <- pi/prev_n[index] %% (2*pi) + offset
          # l <- offset + seq(0, 2*pi, length.out = prev_n[index])
          l <- offset + (2* pi * (0:(prev_n[index]-1)) / prev_n[index])
        }

        x <- 1*radius + sin(l) * (index-1)*r
        y <- 0 + cos(l) * (index-1)*r
        ly_tmp <- data.frame(x = x,
                             y = y)

        ly <- dplyr::bind_rows(ly,ly_tmp)

        ly_list[[i]] <- ly

      }
    }
  }



  ly <- do.call(rbind, ly_list)

  # ---- 统一旋转（绕原点）----
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(ly[, c("x","y")])
    ly[, c("x","y")] <- t(Rm %*% t(xy))
  }

  return(ly)


}
