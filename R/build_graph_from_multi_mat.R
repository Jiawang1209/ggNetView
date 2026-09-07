#' Build a graph object from multiple numeric matrices
#'
#' @param mat1 Numeric matrix.
#' A numeric matrix with samples in colums and variables in rows
#' @param mat2 Numeric matrix.
#' A numeric matrix with samples in colums and variables in rows
#' @param ... Additional numeric matrices with samples in columns and variables in rows.
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#' @param node_annotation Data frame.
#' Optional node annotation table, containing metadata such as taxonomy or functional categories.
#' @param directed Logical (default: \code{FALSE}).
#' Whether edges between nodes are directed.
#' @param top_modules  Integer.
#' Number of top-ranked modules to retain for downstream visualization or analysis.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#'
#' @returns A graph object representing the correlation-based multi numeric matrix.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' mat1 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
#' mat2 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
#' mat3 <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
#' rownames(mat1) <- paste0("A", seq_len(10))
#' rownames(mat2) <- paste0("B", seq_len(10))
#' rownames(mat3) <- paste0("C", seq_len(10))
#' colnames(mat1) <- colnames(mat2) <- colnames(mat3) <-
#'   paste0("sample", seq_len(20))
#' obj <- build_graph_from_multi_mat(
#'   mat1          = mat1,
#'   mat2          = mat2,
#'   mat3,
#'   module.method = "Fast_greedy"
#' )
#' obj
#' }
build_graph_from_multi_mat <- function(mat1,
                                       mat2,
                                       ...,
                                       module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                       node_annotation = NULL,
                                       directed = FALSE,
                                       top_modules = 15,
                                       seed = 1115) {
  mats <- c(list(mat1, mat2), list(...))
  if (length(mats) < 2) {
    stop("At least two matrices are required.")
  }

  dfs <- lapply(mats, function(mat) {
    mat %>% t() %>% as.data.frame()
  })

  common_samples <- Reduce(intersect, lapply(dfs, rownames))
  if (length(common_samples) == 0) {
    stop("No shared samples across matrices.")
  }
  dfs <- lapply(dfs, function(df) df[common_samples, , drop = FALSE])

  all_names <- unlist(lapply(dfs, colnames), use.names = FALSE)
  if (any(duplicated(all_names))) {
    stop("Duplicated feature names across matrices are not supported.")
  }

  df_all <- do.call(cbind, dfs)

  cor_out_odata <- psych::corr.test(df_all, adjust = "none")

  cor_out_odata_r <- cor_out_odata$r %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "from") %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "Correlation")

  cor_out_odata_p <- cor_out_odata$p %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "from") %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "Pvalue") %>%
    dplyr::mutate(signif = dplyr::case_when(
      Pvalue > 0.05                    ~ "",
      Pvalue > 0.01  & Pvalue <= 0.05  ~ "*",
      Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
      Pvalue < 0.001                   ~ "***"
    ))

  mat_names <- names(mats)
  if (is.null(mat_names) || any(mat_names == "")) {
    mat_names <- paste0("mat", seq_along(mats))
  }
  var_group <- purrr::map2_dfr(mat_names, dfs, function(nm, df) {
    tibble::tibble(var = colnames(df), group = nm)
  })

  cor_out <- cor_out_odata_r %>%
    dplyr::left_join(cor_out_odata_p, by = c("from", "to")) %>%
    dplyr::left_join(var_group, by = c("from" = "var")) %>%
    dplyr::rename(from_group = group) %>%
    dplyr::left_join(var_group, by = c("to" = "var")) %>%
    dplyr::rename(to_group = group) %>%
    dplyr::filter(from_group != to_group) %>%
    purrr::set_names(c("from", "to", "Correlation", "Pvalue", "Signif", "from_group", "to_group")) %>%
    dplyr::mutate(Correlated = dplyr::case_when(
      Correlation > 0 ~ "Positive",
      Correlation < 0 ~ "Negative"
    )) %>%
    dplyr::mutate(Signif2 = dplyr::case_when(
      Pvalue > 0.05                    ~ "P > 0.05",
      Pvalue > 0.01  & Pvalue <= 0.05  ~ "0.01 < P <= 0.05",
      Pvalue <= 0.01 & Pvalue >= 0.001 ~ "0.001 <= P <= 0.01",
      Pvalue < 0.001                   ~ "P < 0.001"
    ))

  df <- cor_out %>%
    dplyr::select(from, to, Correlation) %>%
    dplyr::rename(weight = Correlation)

  g <- igraph::graph_from_data_frame(
    d = df,
    vertices = node_annotation,
    directed = directed
  )


  g <- igraph::simplify(g)


  g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))


  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  igraph::E(g)$corr_direction <- ifelse(igraph::E(g)$correlation > 0, "Positive", "Negative")


  module.method <- match.arg(module.method)
  # honour the documented `seed`: community detection (e.g. Spinglass) is
  # stochastic, so seed before it to keep module assignment reproducible.
  .ggnv_local_seed(seed)
  membership_vec <- switch(
    module.method,
    Fast_greedy = igraph::membership(igraph::cluster_fast_greedy(g)),
    Walktrap = igraph::membership(igraph::cluster_walktrap(g)),
    Edge_betweenness = igraph::membership(igraph::cluster_edge_betweenness(g)),
    Spinglass = igraph::membership(igraph::cluster_spinglass(g))
  )

  igraph::V(g)$modularity  <- membership_vec
  igraph::V(g)$modularity2 <- as.character(membership_vec)

  table(igraph::V(g)$modularity2) %>% sort(., decreasing = TRUE)

  # max model length
  max_model <- length(table(igraph::V(g)$modularity2) %>% sort(., decreasing = TRUE))

  if (max_model < top_modules) {
    message(paste("The max module in network is", max_model, "we use the", max_model, " modules for next analysis"))
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(max_model)] %>% names()
  } else if (max_model >= top_modules) {
    modularity_top_15 <- igraph::V(g)$modularity2 %>% table() %>% sort(., decreasing = TRUE) %>% .[seq_len(top_modules)] %>% names()
  }

  igraph::V(g)$modularity2 <- ifelse(igraph::V(g)$modularity2 %in% modularity_top_15, igraph::V(g)$modularity2, "Others")


  graph_obj <- tidygraph::as_tbl_graph(g) %>%
    tidygraph::mutate(modularity = factor(modularity),
                      modularity2 = factor(modularity2),
                      modularity3 = as.character(modularity2),
                      Modularity = modularity2,
                      Degree = tidygraph::centrality_degree(mode = "out"),
                      Segree = tidygraph::centrality_degree(mode = "out"),
                      Strength = tidygraph::centrality_degree(weights = weight)
    ) %>%
    tidygraph::arrange(modularity2, desc(Degree))

  return(graph_obj)
}
