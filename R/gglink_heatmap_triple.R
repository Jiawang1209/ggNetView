#' Visualize multi-orientation environmental-species correlation heatmaps2
#'
#' @param Environment character or data.frame
#' File path or data frame of environment data.
#' @param Experiment character or data.frame
#' File path or data frame of experiment data.
#' @param edge character or data.frame
#' File path or data frame of edge data. Must contain columns \code{from} and
#' \code{to}; an optional numeric \code{weight} column is mapped to edge
#' colour/width (defaults to \code{1} when absent).
#' @param node character or data.frame
#' File path or data frame of node data. Must contain a \code{node} column
#' listing every node referenced by \code{edge}; node names matching columns
#' of \code{Experiment} become the hub nodes anchored on the central heatmap.
#' An optional \code{annotation} column drives node fill/shape (when absent
#' it is derived automatically: \code{"Experiment"} for hub nodes,
#' \code{"Environment"} otherwise).
#' @param sample_col Character (default = "Sample")
#' Column name used as sample ID when input is a data frame or file.
#' @param delim Character (default = ",")
#' Delimiter for reading input files.
#' @param hub_n Integer (default = NULL)
#' If \code{NULL} (recommended), hubs are the \code{Experiment} variables
#' present in \code{node}. If an integer, the \code{hub_n} highest
#' out-degree nodes are used instead (they must then correspond one-to-one
#' to the Experiment variables, and \code{node} rows must list circle nodes
#' first).
#' @param r numeric (default = 6)
#' Radius of the outer node circle.
#' @param cor.method Character (default = "pearson")
#' Correlation method passed to \code{psych::corr.test()}, used for both the
#' Environment x Environment heatmap and the Environment x Experiment links.
#' One of \code{"pearson"}, \code{"kendall"}, \code{"spearman"}.
#' @param cor.use Character (default = "pairwise")
#' Missing-value handling passed to \code{psych::corr.test()}; same vocabulary as
#' \code{gglink_heatmaps()}. Note the default differs from that function
#' (\code{"everything"}) because \code{psych::corr.test()} itself defaults to
#' \code{"pairwise"}, which is what this plot has always used.
#' @param env_p_adjust Character (default = "none")
#' Multiple-testing correction for the Environment x Environment correlations
#' (the significance stars on the triangular heatmap). Any method accepted by
#' \code{stats::p.adjust()}, or \code{"none"}.
#' @param link_p_adjust Character (default = "none")
#' Multiple-testing correction for the Environment x Experiment correlations
#' (the linetype of the link segments). Any method accepted by
#' \code{stats::p.adjust()}, or \code{"none"}.
#' Note that \code{psych::corr.test()} defaults to \code{"holm"} here but still
#' reports raw p-values in \code{$p}, so the previous hard-coded call was in
#' effect uncorrected; \code{"none"} keeps that behaviour.
#' @param sig_breaks Numeric vector of length 3 (default = c(0.05, 0.01, 0.001))
#' Strictly decreasing p-value cut points shared by the heatmap stars
#' (\code{""} / \code{"*"} / \code{"**"} / \code{"***"}) and the link-segment
#' linetype legend.
#'
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' # Environment / Experiment: samples in rows (with a Sample column),
#' # variables in columns. Edges connect Experiment variables (hubs) to
#' # any other nodes.
#' p <- gglink_heatmap_triple(
#'   Environment = env_df,   # Sample + environmental variables
#'   Experiment  = exp_df,   # Sample + experiment variables (become hubs)
#'   edge        = data.frame(from = c("ExpA", "ExpB"),
#'                            to   = c("pH", "TN"),
#'                            weight = c(0.8, 0.5)),
#'   node        = data.frame(node = c("pH", "TN", "ExpA", "ExpB"))
#' )
#' }
gglink_heatmap_triple <- function(
    Environment,
    Experiment,
    edge,
    node,
    sample_col = "Sample",
    delim = ",",
    hub_n = NULL,
    r = 6,
    cor.method = c("pearson", "kendall", "spearman"),
    cor.use = c("pairwise", "everything", "all", "complete", "na"),
    env_p_adjust = "none",
    link_p_adjust = "none",
    sig_breaks = c(0.05, 0.01, 0.001)
){

  cor.method <- match.arg(cor.method)
  cor.use <- match.arg(cor.use)
  p_adjust_choices <- c(stats::p.adjust.methods, "none")
  if (!is.character(env_p_adjust) || length(env_p_adjust) != 1 ||
      !env_p_adjust %in% p_adjust_choices) {
    stop("`env_p_adjust` must be one of: ",
         paste(unique(p_adjust_choices), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.character(link_p_adjust) || length(link_p_adjust) != 1 ||
      !link_p_adjust %in% p_adjust_choices) {
    stop("`link_p_adjust` must be one of: ",
         paste(unique(p_adjust_choices), collapse = ", "), ".", call. = FALSE)
  }

  .read_table <- function(x) {
    if (is.character(x)) {
      readr::read_delim(file = x, delim = delim)
    } else if (is.data.frame(x)) {
      x
    } else {
      stop("Inputs must be file paths or data frames.")
    }
  }

  # Environment Data
  Environment <- .read_table(Environment) %>%
    tibble::as_tibble()
  if (!sample_col %in% colnames(Environment)) {
    stop("`sample_col` not found in Environment.")
  }
  Environment <- Environment %>%
    tibble::column_to_rownames(var = sample_col)

  # Experiment Data
  Experiment <- .read_table(Experiment) %>%
    tibble::as_tibble()
  if (!sample_col %in% colnames(Experiment)) {
    stop("`sample_col` not found in Experiment.")
  }
  Experiment <- Experiment %>%
    tibble::column_to_rownames(var = sample_col)

  # edge Data
  edge <- .read_table(edge) %>%
    tibble::as_tibble()
  if (!all(c("from", "to") %in% colnames(edge))) {
    stop("`edge` must contain columns: from, to.")
  }

  # node Data
  node <- .read_table(node) %>%
    tibble::as_tibble()
  if (!"node" %in% colnames(node)) {
    stop("`node` must contain column: node.")
  }

  # ---- input hardening (0.2.0) --------------------------------------------
  # The plot maps edge weight and node annotation; provide sensible defaults
  # instead of failing at render time with obscure ggplot2 errors.
  if (!"weight" %in% colnames(edge)) {
    edge$weight <- 1
  }
  exp_vars <- colnames(Experiment)
  if (!"annotation" %in% colnames(node)) {
    node$annotation <- ifelse(node$node %in% exp_vars, "Experiment", "Environment")
  }

  # Hub nodes are anchored onto the Experiment rows of the central heatmap,
  # so there must be exactly one hub per Experiment variable appearing in the
  # node table.  Default (`hub_n = NULL`): hubs are the Experiment variables
  # themselves.  The historical `hub_n = NULL -> every node is a hub` default
  # made the outer circle empty and crashed in create_layout2().
  hub_names <- NULL
  if (is.null(hub_n)) {
    hub_names <- intersect(node$node, exp_vars)
    if (length(hub_names) == 0L) {
      stop("None of `node$node` matches a column of `Experiment`. ",
           "Hub nodes must be Experiment variables; add them to `node`, ",
           "or select hubs by degree via `hub_n`.", call. = FALSE)
    }
    if (length(hub_names) == nrow(node)) {
      stop("All nodes are Experiment variables, so no node is left for the ",
           "outer circle. Add non-Experiment nodes to `node`.", call. = FALSE)
    }
    # create_layout2() assigns circle coordinates to the first rows and hub
    # coordinates to the remaining rows; hub anchor coordinates follow the
    # Experiment column order.  Reorder accordingly so users don't have to.
    node <- dplyr::bind_rows(
      node[!node$node %in% hub_names, , drop = FALSE],
      node[match(intersect(exp_vars, node$node), node$node), , drop = FALSE]
    )
    hub_names <- intersect(exp_vars, node$node)
  }

  # Correlation
  stat_out <- cor_test2(Environment,
                        Experiment,
                        cor.method = cor.method,
                        cor.use = cor.use,
                        env_p_adjust = env_p_adjust,
                        link_p_adjust = link_p_adjust,
                        sig_breaks = sig_breaks)

  n_hub_slots <- stat_out[[3]] %>% dplyr::distinct(Experiment) %>% nrow()
  n_hubs <- if (!is.null(hub_names)) length(hub_names) else min(hub_n, nrow(node))
  if (n_hubs != n_hub_slots) {
    stop(sprintf(paste0(
      "Number of hub nodes (%d) must equal the number of Experiment ",
      "variables in the correlation table (%d): each hub is anchored onto ",
      "one Experiment row of the central heatmap. Include all Experiment ",
      "variables in `node` (recommended, with `hub_n = NULL`), or pass a ",
      "matching `hub_n`."), n_hubs, n_hub_slots), call. = FALSE)
  }
  if (nrow(node) - n_hubs < 1L) {
    stop("At least one non-hub node is required for the outer circle.",
         call. = FALSE)
  }

  graph_obj <- tidygraph::tbl_graph(nodes = node, edges = edge)

  # layout
  layout_manual <- create_layout2(graph_obj,
                                 stat_out = stat_out,
                                 hub_names = hub_names,
                                 hub_n = hub_n,
                                 r = r)

  hm_df <- stat_out[[1]]
  id_lab <- hm_df %>%
    dplyr::distinct(ID, ID2, .keep_all = TRUE) %>%
    dplyr::mutate(
      x_lab = ID2,
      y_lab = max(Type2, na.rm = TRUE) + 1
    )
  type_lab <- hm_df %>%
    dplyr::distinct(Type, Type2, .keep_all = TRUE) %>%
    dplyr::mutate(
      x_lab = max(ID2, na.rm = TRUE) + 1,
      y_lab = Type2
    )

  p <- ggraph::ggraph(layout_manual)  +
    ggraph::geom_edge_link(aes(color = weight, width = weight)) +
    ggraph::scale_edge_color_gradientn(colors = c("#74add1","#abd9e9","#ffffbf","#fdae61","#f46d43"),
                               guide = ggraph::guide_edge_colorbar(direction = "horizontal",
                                                           title.position = "top")) +
    ggraph::scale_edge_width(range = c(0.1, 1),
                     guide = ggplot2::guide_legend(nrow = 2,
                                          direction = "horizontal",
                                          title.position = "top"))  +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(data = hm_df,
              aes(x = ID2, y = Type2), fill = "white", color = "#000000", linewidth = 0.5, inherit.aes = FALSE) +
    ggplot2::geom_text(data = id_lab,
              aes(x = x_lab, y = y_lab, label = ID), inherit.aes = FALSE) +
    ggplot2::geom_text(data = type_lab,
              aes(x = x_lab, y = y_lab, label = Type), hjust = "left") +
    ggplot2::geom_point(data = hm_df,
               aes(x = ID2, y = Type2, fill = Value, size = abs(Value)), shape = 21, color = "black") +
    ggplot2::geom_text(data = stat_out[[2]],
              aes(x = ID2, y = Type2, label = p_value),
              size = 7.5) +
    ggplot2::scale_fill_gradient(low = "#edf8b1", high = "#1d91c0", name = "Env Cor",
                        guide = ggplot2::guide_colorbar(direction = "horizontal",
                                               title.position = "top")) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_size(range = c(6,16),
               guide = ggplot2::guide_legend(direction = "horizontal",
                                    title.position = "top"),
               name = "Env Cor Size") +
    ggplot2::xlab('') +
    ggplot2::ylab('')  +
    ggplot2::geom_segment(data = stat_out[[3]] ,
                 mapping = aes(x = p_start1, y = p_end1, xend = start2, yend = end2,
                               linetype = p_value,
                               color = Value,
                               linewidth = abs(Value))) +
    ggplot2::scale_color_gradient(low = "#e0f3db", high = "#4eb3d3", name = "Correlation",
                         guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                title.position = "top")) +
    ggplot2::scale_linewidth(range = c(1, 2.5),
                    name = "Cor",
                    guide = ggplot2::guide_legend(direction = "horizontal",
                                         nrow = 2,
                                         title.position = "top")) +
    ggplot2::scale_linetype(name = "PValue",
                   guide = ggplot2::guide_legend(direction = "horizontal",
                                        nrow = 2,
                                        title.position = "top")) +
    ggplot2::geom_point(data = stat_out[[3]],
               mapping = aes(x = p_start1, y = p_end1), fill = "#9e9ac8", size = 5, shape = 21) +
   ggnewscale::new_scale_fill() +
    ggraph::geom_node_point(data = layout_manual, aes(fill = annotation, shape = annotation), size = 13.5) +
    ggplot2::scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"),
                      guide = ggplot2::guide_legend(title.position = "top",
                                           nrow = 2)) +
    ggplot2::scale_shape_manual(values = c(21,21:25),
                       guide = ggplot2::guide_legend(title.position = "top",
                                            nrow = 2)) +
    # `n_points` is the number of non-hub nodes (placed on the outer circle);
    # `create_layout2()` attaches it as an attribute so we don't hard-code the split.
    ggraph::geom_node_text(
      data = layout_manual %>%
        tidygraph::slice(-seq_len(.n_points_attr(layout_manual))),
      aes(label = node)
    ) +
    ggraph::geom_node_text(
      data = layout_manual %>%
        tidygraph::slice(seq_len(.n_points_attr(layout_manual))),
      aes(x =  x + 0.25,
          y =  y,
          label = node,
      ),
      color = "#000000",
      size = 3,
      hjust = 'outward'
    ) +
    ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::coord_equal(clip = "off") +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, family = "bold"),
      plot.margin = ggplot2::margin(1,1,1,1,"cm"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(hjust = 0.5),
      legend.ticks = ggplot2::element_line(color = "#000000"),
      legend.frame = ggplot2::element_rect(color = "#000000")
    )

  return(p)
}
