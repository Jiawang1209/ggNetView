# library(tidyverse)
# library(ggNetView)
# library(ggnewscale)
# library(ggbump)
# data("Envdf_4st")
# data("Spedf")


gglink_heatmaps4 <- function(
    env,
    spec,
    env_select = NULL,
    spec_select = NULL,
    spec_layout = "ringle",
    relation_method = c("correlation", "mantel"),
    cor.method = c("pearson", "kendall", "spearman"),
    cor.use = c("everything", "all", "complete", "pairwise", "na"),
    mantel.method = c("mantel", "mantel.partial", "mantelhaen.test", "mantel.correlog"),
    mantel.method2 = c("pearson", "kendall", "spearman"),
    mantel.alternative = c("two.sided", "less", "greater"),
    drop_nonsig = FALSE,
    shape = 22,
    distance = 3,
    orientation = c("top_right", "bottom_right", "top_left","bottom_left"),
    r = 6
){

  # argument test
  relation_method <- match.arg(relation_method)
  cor.method      <- match.arg(cor.method)
  cor.use         <- match.arg(cor.use)
  mantel.method   <- match.arg(mantel.method)
  mantel.method2  <- match.arg(mantel.method2)
  mantel.alternative <- match.arg(mantel.alternative)
  orientation     <- match.arg(orientation, several.ok = TRUE)


  # test
  env = Envdf_4st
  spec = Spedf
  orientation = c("top_right", "bottom_right", "top_left","bottom_left")
  # orientation = c("top_right","bottom_left")
  # orientation = "top_right"
  # orientation = "bottom_right"
  cor.method = "pearson"
  cor.use = "pairwise"
  distance = 4
  method = "correlation"
  r = 6



  radius = r
  # if env_select = NULL & spec_select = NULL
  # 说明是最简单的方式 1个点，1个矩阵

  # 如果env_select 含有多个，则出现不同位置的热图

  # 如果 spec_select 出现多个，则出现不同位置的点，这里建议使用环状布局

  ####----split data----####
  # 1 个点， 然后4个环境数据
  spec_select = list(Spec01 = 1:8)

  # different env
  # env_select = list(Env01 = 1:14,
  #                   Env02 = 15:26, # 15:28
  #                   Env03 = 29:38, # 29:42
  #                   Env04 = 43:50 # 43:56
  # )

  # equal env
  env_select = list(Env01 = 1:14,
                    Env02 = 15:28,
                    Env03 = 29:42,
                    Env04 = 43:56)

  # split data
  env_list <- purrr::map(env_select, ~ Envdf_4st[, .x, drop = FALSE])
  spec_list <- purrr::map(spec_select, ~ Spedf[, .x, drop = FALSE])

  # 这里需要统计一下
  k_vec  <- purrr::map_int(env_list, ncol)
  k_ref  <- max(k_vec)

  length_dist <- max(k_vec)  + 0.5 * radius

  # 真实不够的，我们需要使用gap去补充
  k_gap <- length_dist - k_vec

  # 然后分别对env_list的元素，自身内部做相关性分析
  # purrr::map(env_list, function(x){
  #   cor_out_self <- psych::corr.test(x)
  # })

  ####----环境因子自身的相关性----####
  env_cor_self_list <- list()

  # 分别进行分析
  for (i in seq_along(orientation)) {
    # top_right
    if (orientation[i] == "top_right") {
      # 这个是右上角的
      # 单独一个做相关性
      cor_out_self <- psych::corr.test(env_list[[i]], use = cor.use, method = cor.method)

      # correlation
      cor_self_r <- cor_out_self$r %>% as.data.frame()
      cor_self_r[upper.tri(cor_self_r)] <- NA

      # pvalue
      cor_self_p <- cor_out_self$p %>% as.data.frame()
      cor_self_p[upper.tri(cor_self_p)] <- NA

      # combine
      cor_self_r <- cor_self_r %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Correlation") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)
        ) %>%
        stats::na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        stats::na.omit() %>%
        dplyr::mutate(p_signif = dplyr::case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue >= 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6))

      env_cor_self_list[[i]] <- cor_self_r_p

    }

    # bottom_right
    if (orientation[i] == "bottom_right") {
      # 单独一个做相关性
      cor_out_self <- psych::corr.test(env_list[[i]], use = cor.use, method = cor.method)

      # 右下角的
      # correlation
      cor_self_r <- cor_out_self$r %>% as.data.frame()
      cor_self_r[upper.tri(cor_self_r)] <- NA

      # pvalue
      cor_self_p <- cor_out_self$p %>% as.data.frame()
      cor_self_p[upper.tri(cor_self_p)] <- NA
      # combine
      cor_self_r <- cor_self_r %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Correlation") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)
        ) %>%
        stats::na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        stats::na.omit() %>%
        dplyr::mutate(p_signif = dplyr::case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue >= 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6)
      )

      env_cor_self_list[[i]] <- cor_self_r_p
    }

    # top_left
    if (orientation[i] == "top_left") {
      # 单独一个做相关性
      cor_out_self <- psych::corr.test(env_list[[i]], use = cor.use, method = cor.method)

      # correlation
      cor_self_r <- cor_out_self$r %>% as.data.frame()
      cor_self_r[lower.tri(cor_self_r)] <- NA

      # pvalue
      cor_self_p <- cor_out_self$p %>% as.data.frame()
      cor_self_p[lower.tri(cor_self_p)] <- NA

      # combine
      cor_self_r <- cor_self_r %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Correlation") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)
        ) %>%
        stats::na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        stats::na.omit() %>%
        dplyr::mutate(p_signif = dplyr::case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue >= 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6)
      )


      env_cor_self_list[[i]] <- cor_self_r_p
    }

    # bottom_left
    if (orientation[i] == "bottom_left") {
      # 单独一个做相关性
      cor_out_self <- psych::corr.test(env_list[[i]], use = cor.use, method = cor.method)

      # correlation
      cor_self_r <- cor_out_self$r %>% as.data.frame()
      cor_self_r[lower.tri(cor_self_r)] <- NA

      # pvalue
      cor_self_p <- cor_out_self$p %>% as.data.frame()
      cor_self_p[lower.tri(cor_self_p)] <- NA

      # combine
      cor_self_r <- cor_self_r %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Correlation") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)
        ) %>%
        stats::na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        stats::na.omit() %>%
        dplyr::mutate(p_signif = dplyr::case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue >= 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6)
      )


      env_cor_self_list[[i]] <- cor_self_r_p
    }
  }

  env_cor_self_list


  # rename list based on orientation
  names(env_cor_self_list) <- orientation

  # rename k_gap based on orientation
  names(k_gap) <- orientation


  ####----核心物种与环境因子之间的关系----####

  if (method == "correlation") {
    cor_spec_env_list <- list()

    for (p in seq_along(env_list)) {
      for (j in seq_along(spec_list)) {
        # correlation
        cor_env_list_tmp <- psych::corr.test(spec_list[[j]], env_list[[p]])

        cor_env_list_tmp_r <- cor_env_list_tmp$r %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "ID") %>%
          tidyr::pivot_longer(cols = -ID,
                              names_to = "Type",
                              values_to = "Correlation")

        cor_env_list_tmp_p <- cor_env_list_tmp$p %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "ID") %>%
          tidyr::pivot_longer(cols = -ID,
                              names_to = "Type",
                              values_to = "Pvalue") %>%
          dplyr::mutate(p_signif = dplyr::case_when(
            Pvalue > 0.05 ~ "",
            Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
            Pvalue < 0.01 & Pvalue >= 0.001 ~ "**",
            Pvalue < 0.001 ~ "***"
          ))

        cor_env_list_tmp_r_p <- cbind(cor_env_list_tmp_r,
                                      cor_env_list_tmp_p %>% dplyr::select(3,4))


        cor_spec_env_list[[p]] <- cor_env_list_tmp_r_p

      }

    }

    cor_spec_env_list_out <- do.call(rbind, cor_spec_env_list)

  }

  names(cor_spec_env_list) <- orientation

  cor_spec_env_list_out

  # core location layout

  # 查看一下里面有多少个变量, 然后将其均等分
  n_points <- cor_spec_env_list_out$ID %>% unique() %>% length()

  # 计算每一个点的角度
  angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points+1)]
  center_x <- 0
  center_y <- 0

  # 计算坐标
  x <- center_x + radius * cos(angles)
  y <- center_y + radius * sin(angles)

  # 中间点的物理位置
  cor_spec_env <- data.frame(
    ID = cor_spec_env_list_out$ID %>% unique() %>% sort(),
    x = x,
    y = y
  )

  k_vec
  k_gap
  length_dist

  # get targets informations
  # 针对每一个 orientation 设置 x_to y_to
  .make_targets <- function(df, ori, k_gap, length_dist){
    df %>%
      dplyr::mutate(
        ID = as.character(ID),
        Type = as.character(Type)
      ) %>%
      dplyr::filter(ID == Type) %>%
      dplyr::transmute(ID,
                       x_to = dplyr::case_when(
                         ori == "top_right" ~ ID2 + k_gap[ori],
                         ori == "bottom_right" ~ ID2 + k_gap[ori],
                         ori == "top_left" ~ ID2 - length_dist,
                         ori == "bottom_left" ~ID2 - length_dist
                       ),
                       y_to = dplyr::case_when(
                         ori == "top_right" ~ Type2+k_gap[ori]-1,
                         ori == "bottom_right" ~ Type2-length_dist+1,
                         ori == "top_left" ~ Type2+k_gap[ori]-1,
                         ori == "bottom_left" ~Type2-length_dist+1
                       ))
  }

  xy_targets <- purrr::imap_dfr(
    .x = env_cor_self_list[orientation],
    .f = .make_targets,
    k_gap = k_gap,
    length_dist = length_dist
  )

  xy_targets

  cor_spec_env_location <- cor_spec_env_list_out %>%
    dplyr::mutate(ID = as.character(ID), Type = as.character(Type)) %>%
    dplyr::left_join(cor_spec_env, by = "ID") %>%
    dplyr::left_join(xy_targets, by = c("Type" = "ID"))

  .offset_env <- function(df, ori, k_gap, length_dist){
    stopifnot(ori %in% c("top_right","bottom_right","top_left","bottom_left"))
    df <- df %>% dplyr::mutate(ID = as.character(ID), Type = as.character(Type))

    # tile 坐标（四象限通用规则）
    x_tile <- if (ori %in% c("top_right","bottom_right")) df$ID2 + k_gap[[ori]] else df$ID2 - length_dist
    y_tile <- if (ori %in% c("top_right","top_left"))      df$Type2 + k_gap[[ori]] else df$Type2 - length_dist

    tile <- df %>% dplyr::mutate(x_tile = x_tile, y_tile = y_tile, orientation = ori)

    # 对角点（你原来用于标注主对角）
    diag_df <- df %>% dplyr::filter(ID == Type)
    x_diag  <- if (ori %in% c("top_right","bottom_right")) diag_df$ID2 + k_gap[[ori]] else diag_df$ID2 - length_dist
    y_diag  <- if (ori %in% c("top_right","top_left"))     diag_df$Type2 + k_gap[[ori]] - 1 else diag_df$Type2 - length_dist + 1
    diag    <- diag_df %>% dplyr::transmute(ID, x_diag, y_diag, orientation = ori)

    # 轴标签位置
    y_id_lab   <- if (ori %in% c("top_right","top_left")) length_dist + 1 else 0 - length_dist
    x_type_lab <- if (ori %in% c("top_right","bottom_right")) length_dist + 1 else 0 - length_dist
    hjust_type <- if (ori %in% c("top_right","bottom_right")) "left" else "right"

    id_lab <- df %>%
      dplyr::distinct(ID, .keep_all = TRUE) %>%
      dplyr::transmute(ID,
                       x_id = if (ori %in% c("top_right","bottom_right")) ID2 + k_gap[[ori]] else ID2 - length_dist,
                       y_id = y_id_lab, orientation = ori)

    type_lab <- df %>%
      dplyr::distinct(Type, .keep_all = TRUE) %>%
      dplyr::transmute(Type,
                       x_type = x_type_lab,
                       y_type = if (ori %in% c("top_right","top_left")) Type2 + k_gap[[ori]] else Type2 - length_dist,
                       hjust_type = hjust_type, orientation = ori)

    list(tile = tile, diag = diag, id_lab = id_lab, type_lab = type_lab)
  }

  .add_quadrant_layers <- function(p,
                                   pack,
                                   idx,
                                   scale_name = "Env",
                                   low_pal  = c("#4d9221", "#8073ac", "#4393c3", "#66bd63"),
                                   high_pal = c("#c51b7d", "#e08214", "#d6604d", "#f46d43")){
    tile     <- pack$tile
    diag     <- pack$diag
    id_lab   <- pack$id_lab
    type_lab <- pack$type_lab

    p +
      geom_tile(data = tile, aes(x = x_tile, y = y_tile, fill = Correlation)) +
      geom_text(data = tile, aes(x = x_tile, y = y_tile, label = p_signif), size = 5) +
      geom_text(data = id_lab,   aes(x = x_id,   y = y_id,   label = ID)) +
      geom_text(data = type_lab, aes(x = x_type, y = y_type, label = Type),
                hjust = type_lab$hjust_type[1]) +
      geom_point(data = diag, aes(x = x_diag, y = y_diag),
                 shape = 21, fill = "#de77ae", size = 4) +
      scale_fill_gradient2(
        low = low_pal[idx], mid = "#ffffff", high = high_pal[idx],
        midpoint = 0, name = paste0(scale_name, " ", idx),
        guide = guide_colorbar(order = idx)
      )
  }

  # 先为每个方位算好偏移后的数据包
  packs <- purrr::imap(env_cor_self_list, ~ .offset_env(.x, .y, k_gap, length_dist))

  p0 <- ggplot()

  for (i in seq_along(packs)) {
    if (i > 1) p0 <- p0 + ggnewscale::new_scale_fill()
    p0 <- .add_quadrant_layers(p0, packs[[i]], idx = i, scale_name = "Env")
  }

  p0
  # 统一叠加连线 & 中圈节点
  p1 <- p0 +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = cor_spec_env_location,
      aes(x = x, y = y, xend = x_to, yend = y_to, color = Correlation, linetype = p_signif),
      alpha = 0.5
    ) +
    scale_color_gradient(low = "#fdbb84", high = "#d7301f") +
    geom_point(
      data = dplyr::distinct(cor_spec_env_location, ID, .keep_all = TRUE),
      aes(x = x, y = y), shape = 21, fill = "#41b6c4", size = 8.5
    ) +
    geom_text(
      data = dplyr::distinct(cor_spec_env_location, ID, .keep_all = TRUE),
      aes(x = x, y = y, label = ID), size = 5
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(10,10,10,10),
      aspect.ratio = 1,
      legend.position = "top"
    )

  p1


  p2 <- p0 +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = cor_spec_env_location,
      aes(x = x, y = y, xend = x_to, yend = y_to, color = Correlation, linetype = p_signif),
      alpha = 0.5
    ) +
    scale_color_gradient(low = "#fdbb84", high = "#d7301f") +
    geom_point(
      data = dplyr::distinct(cor_spec_env_location, ID, .keep_all = TRUE),
      aes(x = x, y = y), shape = 21, fill = "#41b6c4", size = 8.5
    ) +
    geom_text(
      data = dplyr::distinct(cor_spec_env_location, ID, .keep_all = TRUE),
      aes(x = x, y = y, label = ID), size = 5
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(10,10,10,10),
      aspect.ratio = 1,
      legend.position = "top"
    )

  p2 <- p0 +
    new_scale_fill() +
    geom_curve(data = cor_spec_env_location,
               mapping = aes(x = x, y = y, xend = x_to, yend = y_to, color = Correlation, linetype = p_signif),
               alpha = 0.5,
               curvature = 0.25
    ) +
    scale_color_gradient(low = "#fdbb84", high = "#d7301f") +
    geom_point(data = cor_spec_env_location %>% dplyr::distinct(ID, .keep_all = T),
               mapping = aes(x = x, y = y, fill = ID),
               shape = 21,
               fill = "#41b6c4",
               size = 8.5) +
    geom_text(data = cor_spec_env_location %>% dplyr::distinct(ID, .keep_all = T),
              mapping = aes(x =x, y = y, label = ID),
              size = 5) +
    # geom_line(data = cor_spec_env_location %>% dplyr::distinct(ID, .keep_all = T) %>% dplyr::select(ID, x, y),
    #           mapping = aes(x = x, y = y, group = 1),
    #           linetype = 1,
    #           linewidth = 1.5,
    #           color = "#41b6c4") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(1,1,1,1,"cm"),
      aspect.ratio = 1,
      legend.position = "top"
    )

  p2

  return(list(p1, p2))

}


