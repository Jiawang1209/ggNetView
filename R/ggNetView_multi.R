#' Visualize network with custom layouts in different samples
#'
#' @param mat Numeric matrix.
#' A numeric matrix with samples in rows and variables in columns.
#' @param group_info DataFrame
#' The group information contains: Sample and Group
#' @param transfrom.method Character.
#'Data transformation methods applied before correlation analysis.
#' Options include:
#' "none" (raw data),
#' "scale" (z-score standardization),
#' "center" (mean centering only),
#' "log2" (log2 transfrom),
#' "log10" (log10 transfrom),
#' "ln" (natural transfrom ),
#' "rrarefy" (random rarefaction using \code{vegan::rrarefy}),
#' "rrarefy_relative" (rarefy then convert to relative abundance).
#' @param r.threshold Numeric.
#' Correlation coefficient threshold; edges are kept only if |r| >= r.threshold.
#' @param p.threshold p.threshold
#' Significance threshold for correlations; edges are kept only if p < p.threshold.
#' @param method Character.
#' Relationship analysis methods.
#' Options include: "WGCNA", "SpiecEasi", "SPARCC" and "cor".
#' @param cor.method Character.
#' Correlation analysis method.
#' Options include "pearson", "kendall", and "spearman".
#' @param proc Character.
#' Correlation p-value adjustment methods.
#' Options include:
#' "holm", "hochberg", "hommel", "bonferroni",
#' "BH", "BY", "fdr", and "none".
#' @param module.method Character.
#' Network community detection (module identification) method.
#' Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and "Spinglass".
#' @param SpiecEasi.method Character.
#' Method used in \code{SpiecEasi} network inference; options include "mb" and "glasso".
#' @param sparcc_R Integer.
#' Number of bootstrap/permutation replicates for SparCC p-values (when \code{method = "SPARCC"}).
#' Default 20.
#' @param node_annotation Data frame.
#' Optional node annotation table, containing metadata such as taxonomy or functional categories.
#' @param top_modules Integer.
#' Number of top-ranked modules to retain for downstream visualization or analysis.
#' @param layout Character string naming the layout passed to
#'   \code{ggNetView()} (e.g. "gephi", "fr", "circle", "square").
#' @param ... Additional arguments passed to \code{\link{ggNetView}()}
#'   (node_*, edge_*, module_label_*, module_outline_*, network_outline_*,
#'   layout geometry, ...). Deprecated pre-0.2.0 names (e.g. \code{fill.by},
#'   \code{pointsize}) are still accepted with a lifecycle warning.
#' @param layout_nrow Integer (default = NULL).
#' Number of layout rows passed to \code{ggNetView} when using consensus-module grid layouts.
#' @param layout_ncol Integer (default = NULL).
#' Number of layout columns passed to \code{ggNetView} when using consensus-module grid layouts.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#' @param nrow Integer (default = NULL).
#' Number of rows in the combined patchwork plot.
#' @param ncol Integer (default = NULL).
#' Number of columns in the combined patchwork plot.
#'
#' @returns  A ggplot object representing the network visualization.
#' @export
#'
#' @examples
#' \dontrun{
#' # `mat` is a numeric matrix (features x samples) and
#' # `group_info` is a data frame with columns Sample and Group.
#' p <- ggNetView_multi(
#'   mat        = mat,
#'   group_info = group_info,
#'   method     = "cor",
#'   layout     = "fr"
#' )
#' }
ggNetView_multi <- function(mat,
                            group_info,
                            transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
                                                 "rrarefy_relative"),
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
                            cor.method = c("pearson", "kendall", "spearman"),
                            proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                            module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                            SpiecEasi.method = c("mb", "glasso"),
                            sparcc_R = 20,
                            node_annotation = NULL,
                            top_modules = 15,
                            layout = NULL,
                            ...,
                            layout_nrow = NULL,
                            layout_ncol = NULL,
                            seed = 1115,
                            nrow = NULL,
                            ncol = NULL
                            ){

  method <- match.arg(method)
  p_list <- list()

  for (g in unique(group_info$Group)) {
    group_info_sub <- group_info %>%
      dplyr::filter(Group %in% g)

    mat_sub <- mat %>%
      as.data.frame() %>%
      dplyr::select(all_of(group_info_sub$Sample)) %>%
      tibble::rownames_to_column(var = "ID") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(sum = sum(dplyr::c_across(where(is.numeric)))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(sum != 0) %>%
      dplyr::select(-sum) %>%
      tibble::column_to_rownames(var = "ID")

    graph <- build_graph_from_mat(
      mat = mat_sub,
      transfrom.method = transfrom.method,
      r.threshold = r.threshold,
      p.threshold = p.threshold,
      method = method,
      cor.method = cor.method,
      proc = proc,
      module.method = module.method,
      SpiecEasi.method = SpiecEasi.method,
      sparcc_R = sparcc_R,
      node_annotation = node_annotation,
      top_modules = top_modules,
      seed = seed
    )

    gv_args <- .ggnv_rename_args(list(...), fn = "ggNetView_multi",
                                 env = environment(), user_env = parent.frame())
    gv_args <- utils::modifyList(
      list(graph_obj = graph,
           layout = layout,
           nrow = layout_nrow,
           ncol = layout_ncol,
           seed = seed),
      gv_args
    )
    p <- do.call(ggNetView, gv_args)

    p_list[[g]] <- p

  }

  p_out <- patchwork::wrap_plots(p_list,
                                 # guides = "collect",
                                 nrow = nrow,
                                 ncol = ncol)

  return(p_out)

}
