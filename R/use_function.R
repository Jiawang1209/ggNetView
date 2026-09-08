#' Correlation backend for `gglink_heatmap_triple()`
#'
#' Returns a three-element list:
#' 1. `cor_self`   - lower-triangle Environment x Environment correlations (heatmap tiles)
#' 2. `cor_self_p` - the same cells as significance stars
#' 3. `cor_out_stat` - Environment x Experiment correlations plus the segment anchor
#'    coordinates consumed by `create_layout2()`
#'
#' @param cor.method,cor.use Passed to [psych::corr.test()].
#' @param env_p_adjust Multiple-testing correction for the Environment x
#'   Environment block. `psych::corr.test()` returns raw p-values below the
#'   diagonal of `$p` and corrected ones above it; `$p.adj` is a bare vector for
#'   the single-matrix call and therefore cannot be reused, so the matrix is
#'   transposed to bring the corrected half down when a correction is requested.
#' @param link_p_adjust Multiple-testing correction for the Environment x
#'   Experiment block. Here `$p` stays raw whatever `adjust` is set to and the
#'   corrected values live in `$p.adj`, so the correct matrix is picked explicitly.
#' @param sig_breaks Three strictly decreasing p-value cut points.
#' @noRd
cor_test2 <- function(Environment,
                      Experiment,
                      cor.method = "pearson",
                      cor.use = "pairwise",
                      env_p_adjust = "none",
                      link_p_adjust = "none",
                      sig_breaks = c(0.05, 0.01, 0.001)){

  if (!is.numeric(sig_breaks) || length(sig_breaks) != 3L || anyNA(sig_breaks)) {
    stop("`sig_breaks` must be a numeric vector of length 3.", call. = FALSE)
  }
  if (any(sig_breaks <= 0 | sig_breaks > 1) || any(diff(sig_breaks) >= 0)) {
    stop("`sig_breaks` must be strictly decreasing values in (0, 1], ",
         "e.g. c(0.05, 0.01, 0.001).", call. = FALSE)
  }
  b1 <- sig_breaks[1]; b2 <- sig_breaks[2]; b3 <- sig_breaks[3]
  .fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)
  link_labels <- c(paste0("P > ", .fmt(b1)),
                   paste0(.fmt(b2), " < P <= ", .fmt(b1)),
                   paste0(.fmt(b3), " <= P <= ", .fmt(b2)),
                   paste0("P < ", .fmt(b3)))

  # Environment self
  cor_out_self <- psych::corr.test(Environment,
                                   use = cor.use,
                                   method = cor.method,
                                   adjust = env_p_adjust)
  p_self_mat <- if (identical(env_p_adjust, "none")) {
    cor_out_self$p
  } else {
    t(cor_out_self$p)
  }

  # Environment and Experiment Correlation
  cor_out <- psych::corr.test(Environment,
                              Experiment,
                              use = cor.use,
                              method = cor.method,
                              adjust = link_p_adjust)
  p_link_mat <- if (identical(link_p_adjust, "none")) {
    cor_out$p
  } else {
    matrix(cor_out$p.adj, nrow = nrow(cor_out$p), dimnames = dimnames(cor_out$p))
  }

  cor_out_r <- cor_out$r %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Value")

  cor_out_p <- p_link_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Experiment", values_to = "Pvalue")

  Sample_n <- cor_out_r$Sample %>% unique() %>% length()
  Sample_n_2 <- dim(Experiment)[2]

  cor_out_stat <- cor_out_r %>%
    dplyr::left_join(cor_out_p, by = c("Sample", "Experiment")) %>%
    dplyr::mutate(p_value = dplyr::case_when(
      Pvalue > b1                 ~ link_labels[1],
      Pvalue > b2 & Pvalue <= b1  ~ link_labels[2],
      Pvalue <= b2 & Pvalue >= b3 ~ link_labels[3],
      Pvalue < b3                 ~ link_labels[4]
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

  cor_self_p <- p_self_mat %>% as.data.frame()
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
      Pvalue > b1                 ~ "",
      Pvalue > b2 & Pvalue <= b1  ~ "*",
      Pvalue <= b2 & Pvalue >= b3 ~ "**",
      Pvalue < b3                 ~ "***"
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
