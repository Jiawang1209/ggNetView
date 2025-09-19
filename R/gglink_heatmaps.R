# library(tidyverse)
# library(ggNetView)
# library(ggnewscale)
# library(ggbump)
# data("Envdf_4st")
# data("Spedf")


gglink_heatmaps <- function(
    env,
    spec,
    env_select = NULL,
    spec_select = NULL,
    drop_nonsig = FALSE,
    method = c("correlation", "mantel"),
    cor.method = c("pearson", "kendall", "spearman"),
    cor.use = c("everything", "all", "compelete", "pairwise", "na"),
    mantel.alternative = c("two.sided", "less", "greater"),
    shape = 22,
    distance = 3,
    orientation = c("up","down","left","right", "top_right", "bottom_right", "top_left","bottom_left"),
    r = 6
    ){

  # test
  env = Envdf_4st
  spec = Spedf
  orientation = c("top_right", "bottom_right", "top_left","bottom_left")
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

  env_select = list(Env01 = 1:14,
                    Env02 = 15:28,
                    Env03 = 29:42,
                    Env04 = 43:56)
  # split data
  env_list <- purrr::map(env_select, ~ Envdf_4st[, .x, drop = FALSE])
  spec_list <- purrr::map(spec_select, ~ Spedf[, .x, drop = FALSE])

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
        na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        na.omit() %>%
        dplyr::mutate(p_signif = case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue > 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6))

      env_cor_self_list[[i]] <- cor_self_r_p

      # ggplot(data = cor_self_r_p) +
      #   geom_tile(aes(x = ID2, y = Type2, fill = Correlation)) +
      #   geom_text(aes(x = ID2, y = Type2, label = p_signif), size = 5) +
      #   geom_text(data = cor_self_r_p %>% dplyr::distinct(Type, .keep_all = T),
      #             mapping = aes(x = 15, y = Type2, label = Type), hjust = "left") +
      #   geom_text(data = cor_self_r_p %>% dplyr::distinct(ID, .keep_all = T),
      #             mapping = aes(x = ID2, y = 15, label = ID)) +
      #   geom_point(data = cor_self_r_p %>% dplyr::filter(ID == Type),
      #              mapping = aes(x = ID2, y = Type2-1),
      #              shape = 21,
      #              fill = "#de77ae",
      #              size = 4)



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
        na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        na.omit() %>%
        dplyr::mutate(p_signif = case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue > 0.001 ~ "**",
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
        na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = unique(Type), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        na.omit() %>%
        dplyr::mutate(p_signif = case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue > 0.001 ~ "**",
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
        na.omit()

      cor_self_p <- cor_self_p %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID,
                            names_to = "Type",
                            values_to = "Pvalue") %>%
        dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = T),
                      Type = factor(Type, levels = rev(unique(Type)), ordered = T),
                      ID2 = as.numeric(ID),
                      Type2 = as.numeric(Type)) %>%
        na.omit() %>%
        dplyr::mutate(p_signif = case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
          Pvalue < 0.01 & Pvalue > 0.001 ~ "**",
          Pvalue < 0.001 ~ "***"
        ))

      cor_self_r_p <- cbind(cor_self_r %>% dplyr::select(1,2,4,5,3),
                            cor_self_p %>% dplyr::select(3,6)
      )


      env_cor_self_list[[i]] <- cor_self_r_p
    }
  }

  ####----核心物种与环境因子之间的关系----####
  # if (method == "mantel") {
  #   mantal_spec_env <- vegan::mantel(vegan::vegdist(env_list[[1]]),
  #                                    vegan::vegdist(spec_list[[1]]))
  # }
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
          dplyr::mutate(p_signif = case_when(
            Pvalue > 0.05 ~ "",
            Pvalue > 0.01 & Pvalue < 0.05 ~ "*",
            Pvalue < 0.01 & Pvalue > 0.001 ~ "**",
            Pvalue < 0.001 ~ "***"
          ))

        cor_env_list_tmp_r_p <- cbind(cor_env_list_tmp_r,
                                      cor_env_list_tmp_p %>% dplyr::select(3,4))


        cor_spec_env_list[[p]] <- cor_env_list_tmp_r_p

      }

    }

    cor_spec_env_list_out <- do.call(rbind, cor_spec_env_list)

  }

  cor_spec_env_list_out

  # 查看一下里面有多少个变量, 然后将其均等分
  n_points <- cor_spec_env_list_out$ID %>% unique() %>% length()

  # 计算每一个点的角度
  angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points+1)]
  center_x <- 0
  center_y <- 0

  # 计算坐标
  x <- center_x + radius * cos(angles)
  y <- center_y + radius * sin(angles)

  cor_spec_env <- data.frame(
    ID = cor_spec_env_list_out$ID %>% unique() %>% sort(),
    x = x,
    y = y
  )

  cor_spec_env_location <- cor_spec_env_list_out %>%
    dplyr::left_join(cor_spec_env, by = c("ID" = "ID")) %>%
    dplyr::left_join(rbind(
      env_cor_self_list[[1]] %>%
        dplyr::filter(ID == Type) %>%
        dplyr::select(ID, ID2, Type2) %>%
        dplyr::mutate(ID2 = ID2 + 2,
                      Type2 = Type2+2-1) %>%
        purrr::set_names(c("ID", "x_to", "y_to")),
      env_cor_self_list[[2]] %>%
        dplyr::filter(ID == Type) %>%
        dplyr::select(ID, ID2, Type2) %>%
        dplyr::mutate(ID2 = ID2 + 2,
                      Type2 = Type2-17+1) %>%
        purrr::set_names(c("ID", "x_to", "y_to")),
      env_cor_self_list[[3]] %>%
        dplyr::filter(ID == Type) %>%
        dplyr::select(ID, ID2, Type2) %>%
        dplyr::mutate(ID2 = ID2 - 17,
                      Type2 = Type2+2-1) %>%
        purrr::set_names(c("ID", "x_to", "y_to")),
      env_cor_self_list[[4]] %>%
        dplyr::filter(ID == Type) %>%
        dplyr::select(ID, ID2, Type2) %>%
        dplyr::mutate(ID2 = ID2 - 17,
                      Type2 = Type2-17+1) %>%
        purrr::set_names(c("ID", "x_to", "y_to"))),
      by = c("Type" = "ID"))


  # # test 右上
  # ggplot(data = env_cor_self_list[[1]]) +
  #   geom_tile(aes(x = ID2 + 2, y = Type2 + 2, fill = Correlation)) +
  #   geom_text(aes(x = ID2 + 2, y = 15 + 2, label = ID)) +
  #   geom_text(aes(x = 15 + 2, y = Type2 + 2, label = Type), hjust = "left") +
  #   coord_cartesian(clip = "off") +
  #   theme_bw()
  #
  # # test 右下
  # ggplot(data = env_cor_self_list[[2]]) +
  #   geom_tile(aes(x = ID2 + 2, y = Type2 - 17, fill = Correlation)) +
  #   geom_text(aes(x = ID2 + 2, y = 0 - 17, label = ID)) +
  #   geom_text(aes(x = 15 + 2, y = Type2 - 17, label = Type), hjust = "left") +
  #   coord_cartesian(clip = "off") +
  #   theme_bw()
  #
  # # test 左上
  # ggplot(data = env_cor_self_list[[3]]) +
  #   geom_tile(aes(x = ID2 - 17, y = Type2 + 2, fill = Correlation)) +
  #   geom_text(aes(x = ID2- 17, y = 15 + 2, label = ID)) +
  #   geom_text(aes(x = 0 - 17, y = Type2 + 2, label = Type), hjust = "right") +
  #   coord_cartesian(clip = "off") +
  #   theme_bw()
  #
  # # test 左下
  # ggplot(data = env_cor_self_list[[4]]) +
  #   geom_tile(aes(x = ID2 - 17, y = Type2 - 17, fill = Correlation)) +
  #   geom_text(aes(x = ID2 - 17, y = 0 - 17, label = ID)) +
  #   geom_text(aes(x = 0 - 17, y = Type2 - 17, label = Type), hjust = "right") +
  #   coord_cartesian(clip = "off") +
  #   theme_bw()


  ####----Plot----####
  p1 <- ggplot(data = env_cor_self_list[[1]]) +
    geom_tile(data = env_cor_self_list[[1]], mapping = aes(x = ID2 + 2, y = Type2 + 2, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[1]], mapping = aes(x = ID2 + 2, y = Type2 + 2, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[1]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2 + 2, y = 15 + 2, label = ID)) +
    geom_text(data = env_cor_self_list[[1]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 15 + 2, y = Type2 + 2, label = Type), hjust = "left") +
    geom_point(data = env_cor_self_list[[1]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2+2, y = Type2+2-1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#4d9221", mid = "#ffffff", high = "#c51b7d", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[2]], mapping = aes(x = ID2 + 2, y = Type2 - 17, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[2]], mapping = aes(x = ID2 + 2, y = Type2 - 17, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[2]] %>% dplyr::distinct(ID, .keep_all = T),, mapping = aes(x = ID2 + 2, y = 0 - 17, label = ID)) +
    geom_text(data = env_cor_self_list[[2]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 15 + 2, y = Type2 - 17, label = Type), hjust = "left") +
    geom_point(data = env_cor_self_list[[2]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2+2, y = Type2-17+1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#8073ac", mid = "#ffffff", high = "#e08214", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[3]], mapping = aes(x = ID2 - 17, y = Type2 + 2, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[3]], mapping = aes(x = ID2 - 17, y = Type2 + 2, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[3]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2- 17, y = 15 + 2, label = ID)) +
    geom_text(data = env_cor_self_list[[3]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 0 - 17, y = Type2 + 2, label = Type), hjust = "right") +
    geom_point(data = env_cor_self_list[[3]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2-17, y = Type2+2-1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#4393c3", mid = "#ffffff", high = "#d6604d", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[4]], mapping = aes(x = ID2 - 17, y = Type2 - 17, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[4]], mapping = aes(x = ID2 - 17, y = Type2 - 17, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[4]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2 - 17, y = 0 - 17, label = ID)) +
    geom_text(data = env_cor_self_list[[4]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 0 - 17, y = Type2 - 17, label = Type), hjust = "right") +
    scale_fill_gradient2(low = "#66bd63", mid = "#ffffff", high = "#f46d43", midpoint = 0) +
    geom_point(data = env_cor_self_list[[4]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2-17, y = Type2-17+1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    new_scale_fill() +
    geom_segment(data = cor_spec_env_location,
                 mapping = aes(x = x, y = y, xend = x_to, yend = y_to, color = Correlation, linetype = p_signif),
                 alpha = 0.5
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


  p2 <- ggplot(data = env_cor_self_list[[1]]) +
    geom_tile(data = env_cor_self_list[[1]], mapping = aes(x = ID2 + 2, y = Type2 + 2, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[1]], mapping = aes(x = ID2 + 2, y = Type2 + 2, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[1]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2 + 2, y = 15 + 2, label = ID)) +
    geom_text(data = env_cor_self_list[[1]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 15 + 2, y = Type2 + 2, label = Type), hjust = "left") +
    geom_point(data = env_cor_self_list[[1]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2+2, y = Type2+2-1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#4d9221", mid = "#ffffff", high = "#c51b7d", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[2]], mapping = aes(x = ID2 + 2, y = Type2 - 17, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[2]], mapping = aes(x = ID2 + 2, y = Type2 - 17, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[2]] %>% dplyr::distinct(ID, .keep_all = T),, mapping = aes(x = ID2 + 2, y = 0 - 17, label = ID)) +
    geom_text(data = env_cor_self_list[[2]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 15 + 2, y = Type2 - 17, label = Type), hjust = "left") +
    geom_point(data = env_cor_self_list[[2]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2+2, y = Type2-17+1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#8073ac", mid = "#ffffff", high = "#e08214", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[3]], mapping = aes(x = ID2 - 17, y = Type2 + 2, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[3]], mapping = aes(x = ID2 - 17, y = Type2 + 2, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[3]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2- 17, y = 15 + 2, label = ID)) +
    geom_text(data = env_cor_self_list[[3]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 0 - 17, y = Type2 + 2, label = Type), hjust = "right") +
    geom_point(data = env_cor_self_list[[3]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2-17, y = Type2+2-1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
    scale_fill_gradient2(low = "#4393c3", mid = "#ffffff", high = "#d6604d", midpoint = 0) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = env_cor_self_list[[4]], mapping = aes(x = ID2 - 17, y = Type2 - 17, fill = Correlation)) +
    geom_text(data = env_cor_self_list[[4]], mapping = aes(x = ID2 - 17, y = Type2 - 17, label = p_signif), size = 5) +
    geom_text(data = env_cor_self_list[[4]] %>% dplyr::distinct(ID, .keep_all = T), mapping = aes(x = ID2 - 17, y = 0 - 17, label = ID)) +
    geom_text(data = env_cor_self_list[[4]] %>% dplyr::distinct(Type, .keep_all = T), mapping = aes(x = 0 - 17, y = Type2 - 17, label = Type), hjust = "right") +
    scale_fill_gradient2(low = "#66bd63", mid = "#ffffff", high = "#f46d43", midpoint = 0) +
    geom_point(data = env_cor_self_list[[4]] %>% dplyr::filter(ID == Type),
               mapping = aes(x = ID2-17, y = Type2-17+1),
               shape = 21,
               fill = "#de77ae",
               size = 4) +
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

# ggsave(filename = "/Users/liuyue/Desktop/Github_repos/R&Python_Wechat_Official_Account_Platform/20250909_网络可视化R包ggNetView_添加核心变量与多个环境因子的关系网络图/Output/Out.pdf",
#        height = 13,
#        width = 13,
#        plot = p1)
#
#
# ggsave(filename = "/Users/liuyue/Desktop/Github_repos/R&Python_Wechat_Official_Account_Platform/20250909_网络可视化R包ggNetView_添加核心变量与多个环境因子的关系网络图/Output/Out_2.pdf",
#        height = 13,
#        width = 13,
#        plot = p2)


}


