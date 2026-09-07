#' @noRd
cor_test2 <- function(Environment, Experiment){
  # Environment self
  cor_out_self <- psych::corr.test(Environment, adjust = "none")

  cor_out_self_r <- cor_out_self$r %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Value")

  cor_out_self_p <- cor_out_self$p %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Pvalue")

  cor_out_self_stat <- cor_out_self_r %>%
    dplyr::left_join(cor_out_self_p, by = c("Sample", "Experiment"))

  # Environment and Experiment Correlation
  cor_out <- psych::corr.test(Environment, Experiment)

  cor_out_r <- cor_out$r %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Value")

  cor_out_p <- cor_out$p %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Pvalue")

  Sample_n <- cor_out_r$Sample %>% unique() %>% length()
  Sample_n_2 <- dim(Experiment)[2]

  cor_out_stat <- cor_out_r %>%
    dplyr::left_join(cor_out_p, by = c("Sample", "Experiment")) %>%
    dplyr::mutate(p_value = dplyr::case_when(
      Pvalue > 0.05                    ~ "P > 0.05",
      Pvalue > 0.01  & Pvalue <= 0.05  ~ "0.01 < P <= 0.05",
      Pvalue <= 0.01 & Pvalue >= 0.001 ~ "0.001 <= P <= 0.01",
      Pvalue < 0.001                   ~ "P < 0.001"
    )) %>%
    dplyr::mutate(
      p_start1 = rep(1:Sample_n, each = Sample_n_2),
      p_end1 = rep(rev(0:(Sample_n-1)), each = Sample_n_2),
      # One anchor position per Experiment variable (odd positions 1,3,5,...).
      # The historical `seq(1, Sample_n, 2)[1:Sample_n_2]` ran out of values
      # (-> NA coordinates -> hubs and their link segments silently dropped)
      # whenever ncol(Experiment) > ceiling(ncol(Environment)/2); using
      # length.out keeps the same values in the old range and extends beyond.
      start2 = rep(seq(1, by = 2, length.out = Sample_n_2), times = Sample_n) - 4,
      end2 = rev(rep(seq(1, by = 2, length.out = Sample_n_2), times = Sample_n))
    )

  ####----Plot----####
  cor_self <- cor_out_self$r %>% as.data.frame()
  cor_self[upper.tri(cor_self)] <- NA

  cor_self <- cor_self %>%
    tibble::rownames_to_column(var = "ID") %>%
    tidyr::pivot_longer(cols = -ID,
                        names_to = "Type",
                        values_to = "Value") %>%
    dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = TRUE),
                  Type = factor(Type, levels = rev(unique(Type)), ordered = TRUE),
                  ID2 = as.numeric(ID),
                  Type2 = as.numeric(Type)
    ) %>%
    stats::na.omit()

  cor_self_p <- cor_out_self$p %>% as.data.frame()
  cor_self_p[upper.tri(cor_self_p)] <- NA
  cor_self_p <- cor_self_p %>%
    tibble::rownames_to_column(var = "ID") %>%
    tidyr::pivot_longer(cols = -ID,
                        names_to = "Type",
                        values_to = "Pvalue") %>%
    dplyr::mutate(ID = factor(ID, levels = unique(ID), ordered = TRUE),
                  Type = factor(Type, levels = rev(unique(Type)), ordered = TRUE),
                  ID2 = as.numeric(ID),
                  Type2 = as.numeric(Type)) %>%
    stats::na.omit() %>%
    dplyr::mutate(p_value = dplyr::case_when(
      Pvalue > 0.05                    ~ "",
      Pvalue > 0.01  & Pvalue <= 0.05  ~ "*",
      Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
      Pvalue < 0.001                   ~ "***"
    ))

  return(list(cor_self, cor_self_p, cor_out_stat))
}


#' Read the `n_points` attribute attached by `create_layout2()`; returns 0L if
#' the attribute is missing (R-version-agnostic alternative to `%||%`).
#' @noRd
.n_points_attr <- function(ly) {
  n <- attr(ly, "n_points")
  if (is.null(n) || !is.finite(n)) 0L else as.integer(n)
}


create_layout2 <- function(graph, stat_out, hub_names = NULL, hub_n = NULL, r = 10) {

  nodes <- graph %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  if (is.null(hub_names)) {
    if (is.null(hub_n)) {
      hub_names <- nodes$node
    } else {
      deg_df <- graph %>%
        tidygraph::mutate(degree = tidygraph::centrality_degree(mode = "out")) %>%
        tidygraph::as_tibble() %>%
        tidygraph::arrange(dplyr::desc(degree))

      hub_n <- min(hub_n, nrow(deg_df))
      hub_names <- deg_df$node[seq_len(hub_n)]
    }
  }

  # hub
  hub_names

  # non hub
  non_hub_names <- nodes$node[!nodes$node %in% hub_names]


  n_points <- length(non_hub_names)
  radius <- r
  center_x <- -12
  center_y <- 4


  angles <- seq(0, 2*pi, length.out = n_points + 1)[-(n_points+1)]


  x <- center_x + radius * cos(angles)
  y <- center_y + radius * sin(angles)



  circle_df <- data.frame(
    id = seq_len(n_points),
    x = x,
    y = y
  )

  hub_df <- stat_out[[3]] %>%
    dplyr::distinct(Experiment, .keep_all = TRUE) %>%
    dplyr::select(start2, end2) %>%
    purrr::set_names(c("x", "y"))



  layout_manual <- ggraph::create_layout(graph, layout = "circle")

  layout_manual_2 <- rbind(layout_manual[1:n_points,] %>%
                             dplyr::mutate(x = circle_df$x,
                                           y = circle_df$y),
                           layout_manual[-c(1:n_points),] %>%
                             dplyr::mutate(x = hub_df$x,
                                           y = hub_df$y)
  )

  ly <- ggraph::create_layout(
    graph,
    layout = "manual",
    x = layout_manual_2$x,
    y = layout_manual_2$y
  )

  # carry the non-hub / hub split point so callers (e.g. gglink_heatmap_triple)
  # don't need to re-derive it via a hard-coded slice index.
  attr(ly, "n_points") <- n_points
  attr(ly, "n_hub")    <- length(hub_names)

  return(ly)
}
