#' Visualize multiple networks and show the connections between them to highlight key members within the networks
#'
#' Provide either \code{(mat, group_info)} or \code{graph_obj_list} as input.
#' When \code{graph_obj_list} is provided, \code{mat} and \code{group_info} are ignored.
#'
#' @param mat Numeric matrix. Required when \code{graph_obj_list} is \code{NULL}.
#' A numeric matrix with variables (e.g. genes, taxa) in rows and samples in columns.
#' @param group_info DataFrame. Required when \code{graph_obj_list} is \code{NULL}.
#' The group information contains: Sample and Group.
#' If Sample or Group contain underscores (\code{_}), they are automatically replaced with hyphens (\code{-}) to avoid parsing issues; \code{mat} column names, \code{order}, and \code{comparisons_groups} are updated accordingly.
#' @param graph_obj_list Named list of \code{tbl_graph} objects. Alternative to \code{mat} and \code{group_info}.
#' Each element is a graph object (e.g. from \code{build_graph_from_mat} or \code{build_graph_from_df}).
#' List names define group names (e.g. \code{list(WT = g1, KO = g2)}).
#' When provided, \code{mat} and \code{group_info} are ignored. Enables custom pre-built networks.
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
#' @param p.threshold Numeric.
#' #' Significance threshold for correlations; edges are kept only if p < p.threshold.
#' @param method Character.
#' Relationship analysis methods.
#' Options include: "WGCNA", "SpiecEasi", "SPARCC" and "cor".
#' @param cor.method Character.
#' Correlation analysis method.
#' Options include "pearson", "kendall", and "spearman".
#' @param proc  Character.
#' Correlation p-value adjustment methods.
#' Options include:
#' "holm", "hochberg", "hommel", "bonferroni",
#' "BH", "BY", "fdr", and "none".
#' @param sig_by Character (default \code{"pvalue"}). Which column decides
#'   whether a between-module comparison is kept (< 0.05): the raw
#'   \code{"pvalue"} (default; preserves historical behaviour) or the
#'   BH-adjusted \code{"FDR"} column. The \code{FDR} is Benjamini-Hochberg
#'   adjusted \emph{within each group-pair's} module-vs-module comparison grid
#'   (as computed by \code{compare_modules_by_overlap()}), not pooled globally
#'   across all group-pairs; use \code{"FDR"} to damp false positives from the
#'   many module comparisons made within a pair.
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
#' @param layout Character string.
#' Custom layouts; one of "gephi", "square", "square2", "petal",
#' "petal2", "heart_centered","diamond", "star", "star_concentric","rectangle,
#' "rightiso_layers" etc.
#' @param node_add Integer (default = 7).
#' Number of nodes to add in each layer of the layout.
#' @param r Numeric (default = 1).
#' Radius increment for concentric or layered layouts.
#' @param center Logical (default = TRUE).
#' Whether to place a node at the center of the layout.
#' @param idx Optional.
#' Index of nodes to be emphasized or centered in the layout
#' @param shrink Numeric (default = 1).
#' Shrinkage factor applied to the center points.
#' @param inner_shrink Numeric (default = 1).
#' Intra-module compactness factor for \code{layout = "WGCNA"} only.
#' See \code{\link{ggNetView}} for details.  Ignored by other layouts.
#' @param k_nn Numeric (default = 8).
#' Number of nearest neighbors used to build the local adjacency graph.
#' @param push_others_delta Numeric (default = 0).
#' Radial offset applied to the "Others" module to slightly
#' @param layout_module Character  (default = "random")
#' - random : modules are distributed more randomly and independently.
#' - adjacent : modules are positioned close to each other, minimizing inter-module gaps.
#' - order : modules are distributed by order, applicable to `Bipartite, Tripartite, Quadripartite, Multipartite, Pentapartite Layout`
#' @param group_by Character (default = "Modularity").
#' Change group for nodes
#' @param node_fill Character (default = "Modularity").
#' Change fill for nodes
#' @param node_fill_values Named vector of colors for node/module fill (e.g. \code{c("M1" = "red", "M2" = "blue")}).
#' If \code{NULL} (default), uses viridis discrete fill scale
#' (\code{scale_fill_viridis_d});
#' if provided, uses \code{scale_fill_manual(values = node_fill_values)}.
#' @param node_color_values Color setting for node/module border.
#' Supports either a single color string (fixed border color) or a named vector
#' (module-to-color mapping, similar to \code{node_fill_values}).
#' If \code{NULL}, mapped borders use viridis discrete color scale
#' (\code{scale_color_viridis_d}).
#' @param node_size_range Numeric vector of length 2 (default = \code{c(1, 5)}).
#' The range of point size when mapping \code{Degree} to size.
#' First value is minimum size, second is maximum size.
#' @param node_jitter Logical (default = FALSE).
#' Whether to apply node_jitter to points.
#' @param node_jitter_sd  Integer  (default = 0.1).
#' The standard deviation of the node_jitter applied when `node_jitter = TRUE`.
#' @param edge_color Character (default = "grey70").
#' Within-group edge colour: an edge column name (mapping; \code{"corr_direction"}
#' colours positive/negative edges red/blue and adds the counts to the group
#' label) or a single colour (constant).
#' @param edge_color_values Named colour vector or NULL (default = NULL).
#' Manual palette for a categorical \code{edge_color} mapping.
#' @param edge_alpha  Integer  (default = 0.25).
#' Change  line alpha.
#' @param edge_curve Logical (default = FALSE).
#' Whether to draw within-group edges as curves.
#' @param edge_curvature Numeric (default = 0.12).
#' Curvature for within-group edges when \code{edge_curve = TRUE}.
#' @param edge_curve_adaptive Logical (default = TRUE).
#' Whether to adapt within-group edge curvature by edge length when
#' \code{edge_curve = TRUE}.
#' @param edge_curve_adaptive_range Numeric vector of length 2 (default = c(0.7, 1.3)).
#' Multipliers applied to \code{edge_curvature} for shortest and longest
#' within-group edges.
#' @param edge_curve_adaptive_bins Integer (default = 7).
#' Number of bins used to approximate per-edge adaptive curvature for
#' within-group edges.
#' @param module_outline Logical or Character (default = "circle").
#' Add outer boundaries for matched modules.
#' Supported values: \code{"circle"} (use \code{ggforce::geom_mark_circle}),
#' \code{"manual"} (use smoothed polygon boundary like \code{ggNetView}),
#' and \code{"none"} (disable). Logical \code{TRUE}/\code{FALSE} are
#' accepted and mapped to \code{"circle"}/\code{"none"}.
#' @param module_outline_q Numeric (default = 0.88).
#' HDR coverage of the outer boundary when \code{module_outline = "manual"}: the
#' contour is drawn at the density level whose iso-density region contains a
#' fraction \code{module_outline_q} of the module's empirical probability mass.
#' @param module_outline_expand Numeric (default = 1.02).
#' Multiplicative scaling applied to each polygon from its own centroid
#' when \code{module_outline = "manual"}.
#' @param module_outline_bandwidth Numeric (default = 1.0).
#' Multiplier on the robust normal-reference 2D KDE bandwidth used to build
#' the outer boundary when \code{module_outline = "manual"}.
#' @param module_outline_width Numeric (default = 1.25).
#' Line width for module outer boundaries.
#' @param module_outline_linetype Integer or character (default = 2).
#' Linetype for module outer boundaries (e.g. 1 = solid, 2 = dashed).
#' @param module_outline_alpha Numeric (default = 0.5).
#' Alpha for module outer boundaries.
#' @param link_level Character (default = "Module").
#' Cross-group link granularity. One of \code{"None"}, \code{"Module"}, \code{"Node"},
#' \code{"NodeinModule"}, \code{"Module&Node"}, or \code{"Module&Node2"}.
#' \code{"None"} draws no cross-group links.
#' \code{"Module"} links significant module matches (module centroids).
#' \code{"Node"} links shared node names across groups.
#' \code{"NodeinModule"} links both module centroids and nodes within significantly overlapping
#' modules (module links show which modules overlap; node links show shared nodes).
#' \code{"Module&Node"} draws both module-level links and node-level links.
#' \code{"Module&Node2"} same as \code{"Module&Node"} but adds outer boundaries for all modules
#' (not only matched ones).
#' @param link_curve Logical (default = FALSE).
#' Whether to draw cross-group links as curves (\code{geom_curve}) instead of straight segments.
#' @param link_curvature Numeric (default = 0.2).
#' Curvature used when \code{link_curve = TRUE}.
#' @param link_curve_mode Character (default = "outward").
#' Curve direction strategy used when \code{link_curve = TRUE}.
#' \code{"outward"} bends links away from the global center.
#' \code{"inward"} bends links toward the global center.
#' \code{"cross"} follows a cross-axis rule: left links bend left, right links bend right,
#' upper links bend up, and lower links bend down.
#' @param link_curve_adaptive Logical (default = TRUE).
#' Whether to adapt link curvature by link length when \code{link_curve = TRUE}.
#' Longer links get larger curvature and shorter links get smaller curvature.
#' @param link_curve_adaptive_range Numeric vector of length 2 (default = c(0.7, 1.3)).
#' Multipliers applied to \code{link_curvature} for shortest and longest links.
#' The first value is the minimum multiplier; the second is the maximum multiplier.
#' @param link_curve_adaptive_bins Integer (default = 7).
#' Number of bins used to approximate per-link adaptive curvature.
#' @param link_color_node Character or NULL (default = NULL).
#' Colors for node-to-node cross-group links. \code{NULL} = use default palette.
#' A single value, named vector (e.g. \code{c("WT|KO" = "red")}), or unnamed vector (by pair index).
#' @param link_color_module Character or NULL (default = NULL).
#' Colors for module-to-module cross-group links. Same rules as \code{link_color_node}.
#' @param link_width_node Numeric (default = 1).
#' Line width for node-to-node cross-group links. Single value or vector (by pair index/named).
#' @param link_width_module Numeric (default = 1).
#' Line width for module-to-module cross-group links. Single value or vector.
#' @param link_linetype_node Integer or character (default = 2).
#' Linetype for node-to-node cross-group links (e.g. 2 = dashed, 1 = solid). Single value or vector.
#' @param link_linetype_module Integer or character (default = 1).
#' Linetype for module-to-module cross-group links. Single value or vector.
#' @param link_alpha_node Numeric (default = 0.25).
#' Alpha (transparency) for node-to-node cross-group links. Single value or vector.
#' @param link_alpha_module Numeric (default = 0.5).
#' Alpha (transparency) for module-to-module cross-group links. Single value or vector.
#' @param drop_others Logical (default = FALSE).
#' If TRUE, remove nodes in the \code{"Others"} module from each group's
#' \code{graph_obj} before layout and plotting.
#' This is a display-only switch: the module-overlap comparison always runs on the
#' complete network (including \code{"Others"} nodes), so the set of cross-group
#' module links is identical for \code{drop_others = TRUE} and \code{FALSE};
#' only the plotted nodes differ.
#' Note that \code{"Others"} is a display bucket for every module ranked below
#' \code{top_modules}, not a community, and never takes part in module links.
#' @param calculate_topology Logical (default = FALSE).
#' Whether to compute topology for each group using
#' \code{get_network_topology_parallel()} and
#' \code{get_sample_subgraph_topology_parallel()}.
#' @param comparisons Logical (default = TRUE).
#' Whether to perform cross-group comparisons and draw links.
#' If \code{FALSE}, no module overlap or node comparison is done, and no cross-group links are drawn.
#' @param comparisons_groups List or NULL (default = NULL).
#' When \code{comparisons = TRUE}, constrains which group pairs are compared.
#' Each element must be a length-2 character vector, e.g. \code{list(c("WT", "OE"), c("WT", "KO"))}.
#' Group names must exist in \code{group_info$Group} or \code{names(graph_obj_list)}.
#' If \code{NULL}, all pairwise group comparisons are performed.
#' @param order Character vector or NULL (default = NULL).
#' Order of groups for layout positions. Groups are placed evenly on a circle:
#' 1st = top, 2nd = next clockwise, etc. (e.g. 3 groups = triangle, 4 = square).
#' Must contain all unique groups from \code{group_info$Group} or \code{names(graph_obj_list)} exactly once.
#' If \code{NULL}, uses \code{unique(group_info$Group)} or \code{names(graph_obj_list)} order (first occurrence).
#' @param group_layout Character (default = "circle").
#' Arrangement of groups in the multi-group plot.
#' \code{"circle"}: groups placed evenly on a circle (default).
#' \code{"row"}: groups in a grid, filled row-by-row; uses \code{nrow} and \code{ncol}.
#' \code{"column"}: groups in a grid, filled column-by-column; uses \code{nrow} and \code{ncol}.
#' \code{"square"}: 4 groups at corners of a square (top-left, top-right, bottom-right, bottom-left).
#' \code{"diamond"}: 4 groups at top, right, bottom, left (like a rotated square).
#' \code{"triangle"}: 3 groups at vertices of an upright triangle (point up).
#' \code{"triangle_down"}: 3 groups at vertices of an inverted triangle (point down).
#' \code{"snake"}: groups in a grid with snake-like ordering; first row left-to-right,
#' second row right-to-left, third row left-to-right, etc. Uses \code{nrow} and \code{ncol}.
#' \code{"snake_vertical"}: groups along a smooth sine curve; uses \code{sine_period}.
#' \code{"snake_vertical_sin"}: same as \code{"snake_vertical"}.
#' \code{"snake_vertical_cos"}: groups along a smooth cosine curve; uses \code{sine_period}.
#' \code{"snake_vertical_neg_sin"}: smooth \code{-sin} curve.
#' \code{"snake_vertical_neg_cos"}: smooth \code{-cos} curve.
#' \code{"sin"}: vertex-only; group 1 at center; groups 2+ alternate peaks and troughs.
#' \code{"cos"}: vertex-only; group 1 at peak; groups 2+ alternate peak and trough.
#' \code{"-sin"}: vertex-only; group 1 at center; groups 2+ alternate down, up.
#' \code{"-cos"}: vertex-only; group 1 at trough; groups 2+ alternate trough and peak.
#' \code{"center_pairs"}: group 1 at center (y=0); remaining groups form pairs on peak (y=1) and trough (y=-1),
#' with one \code{anchor_dist} gap between center and first pair. Pairs alternate right/left of center.
#' If a pair has only one group, it goes on top (peak). Empty pair slots keep their positions.
#' @param sine_period Numeric (default = 4).
#' Groups per wavelength for \code{snake_vertical*}; ignored for \code{sin}, \code{cos}, \code{-sin}, \code{-cos}.
#' @param scale_groups Logical (default = TRUE).
#' Whether to normalize each group to a comparable coordinate scale before
#' placing groups on anchors for cross-group visual comparison.
#' @param orientation Character string.
#' Custom orientation; one of "up","down","left","right".
#' @param angle Integer  (default = 0).
#' Change  orientation angle.
#' @param anchor_dist Numeric (default = 6).
#' Distance between groups when placing each group network on the outer circle
#' in \code{ggNetView_multi_link}.
#' @param layout_anchor_dist Numeric (default = NULL).
#' Anchor distance passed to the single-group layout function (module spacing
#' within each group). If \code{NULL}, it falls back to \code{anchor_dist}
#' for backward compatibility.
#' @param nrow Integer (default = NULL).
#' Number of rows. Used by: (1) \code{group_layout = "row"}, \code{"column"}, or \code{"snake"} for group grid;
#' (2) layout functions like \code{"consensus_module_equal_gephi"} for module grid.
#' @param ncol Integer (default = NULL).
#' Number of columns. Used by: (1) \code{group_layout = "row"}, \code{"column"}, or \code{"snake"} for group grid;
#' (2) layout functions like \code{"consensus_module_equal_gephi"} for module grid.
#' @param group_label_offset Numeric (default = 0.2).
#' Vertical offset of group labels above each group's layout (added to max y).
#' @param group_label_size Numeric (default = 4).
#' Font size for group labels (Group, Node, Edge, etc.).
#' @param network_outline Logical (default = FALSE).
#' Whether to add a circle boundary around each group (mimics \code{ggforce::geom_mark_circle}).
#' @param network_outline_expand Numeric (default = 2).
#' Expansion in mm for the group circle to account for point size; passed to \code{geom_mark_circle(expand = ...)}.
#' @param network_outline_color Character (default = "grey50").
#' Color of the group outer circle border. A single value applies to all groups;
#' a named vector maps group names to colors (e.g. \code{c("WT" = "blue", "KO" = "red")});
#' an unnamed vector is used by index (recycled if needed).
#' @param network_outline_fill Character or NULL (default = NULL).
#' Fill color of the group outer circle. \code{NULL} = no fill (transparent).
#' A single value, named vector, or unnamed vector works like \code{network_outline_color}.
#' @param network_outline_fill_alpha Numeric (default = 0.2).
#' Alpha (transparency) of the group outer circle fill; 0 = fully transparent, 1 = opaque.
#' @param network_outline_linetype Integer or character (default = 1).
#' Linetype of the group outer circle (e.g. 1 = solid, 2 = dashed).
#' @param network_outline_width Numeric (default = 0.5).
#' Line width of the group outer circle.
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility.
#' @param ring_n,layout.module,group.by,fill.by,fill,color,pointsize,jitter,jitter_sd,mapping_line,linealpha,linecolor,inner_curve,inner_curvature,inner_curve_adaptive,inner_curve_adaptive_range,inner_curve_adaptive_bins,add_outer,q_outer,expand_outer,bandwidth_scale,outerwidth,outerlinetype,outeralpha,link_linewidth_node,link_linewidth_module,link_linealpha_node,link_linealpha_module,dropOthers,label_offset,label_size,add_group_outer,add_group_outer_expand,add_group_outer_color,add_group_outer_fill,add_group_outer_fill_alpha,add_group_outer_linetype,add_group_outer_linewidth
#'   \lifecycle{deprecated} Pre-0.2.0 argument names, kept for backward
#'   compatibility. They emit a deprecation warning and are forwarded to the
#'   new argument (same renaming scheme as \code{\link{ggNetView}}; in addition
#'   \code{inner_curve*} -> \code{edge_curve*}, \code{link_linewidth_*} ->
#'   \code{link_width_*}, \code{link_linealpha_*} -> \code{link_alpha_*},
#'   \code{label_offset}/\code{label_size} -> \code{group_label_offset}/
#'   \code{group_label_size}).
#'
#' @returns A list containing plot, module-overlap info, link info, group graphs,
#' and optional topology results.
#' @export
#'
#' @examples
#' \dontrun{
#' # Provide either `mat` + `group_info` or a pre-built `graph_obj_list`.
#' p <- ggNetView_multi_link(
#'   mat        = mat,
#'   group_info = group_info,
#'   method     = "cor"
#' )
#' }
ggNetView_multi_link <- function(mat = NULL,
                                 group_info = NULL,
                                 graph_obj_list = NULL,
                                 transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy", "rrarefy_relative"),
                                 r.threshold = 0.7,
                                 p.threshold = 0.05,
                                 method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
                                 cor.method = c("pearson", "kendall", "spearman"),
                                 proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                 sig_by = c("pvalue", "FDR"),
                                 module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
                                 SpiecEasi.method = c("mb", "glasso"),
                                 sparcc_R = 20,
                                 node_annotation = NULL,
                                 top_modules = 15,
                                 layout = NULL,
                                 node_add = 7,
                                 r = 1,
                                 center = TRUE,
                                 idx = NULL,
                                 shrink = 1,
                                 inner_shrink = 1,
                                 k_nn = 12,
                                 push_others_delta = 0,
                                 layout_module = c("random", "adjacent", "order"),
                                 group_by = "Modularity",
                                 node_fill = "Modularity",
                                 node_fill_values = NULL,
                                 node_color_values = NULL,
                                 node_size_range = c(1, 5),
                                 node_jitter = FALSE,
                                 node_jitter_sd = 0.01,
                                 edge_color = "grey70",
                                 edge_color_values = NULL,
                                 edge_alpha = 0.25,
                                 edge_curve = FALSE,
                                 edge_curvature = 0.12,
                                 edge_curve_adaptive = TRUE,
                                 edge_curve_adaptive_range = c(0.7, 1.3),
                                 edge_curve_adaptive_bins = 7,
                                 module_outline = "circle",
                                 module_outline_q = 0.88,
                                 module_outline_expand = 1.02,
                                 module_outline_bandwidth = 1.0,
                                 module_outline_width = 1.25,
                                 module_outline_linetype = 2,
                                 module_outline_alpha = 0.5,
                                 link_level = "Module",
                                 link_curve = FALSE,
                                 link_curvature = 0.2,
                                 link_curve_mode = "outward",
                                 link_curve_adaptive = TRUE,
                                 link_curve_adaptive_range = c(0.7, 1.3),
                                 link_curve_adaptive_bins = 7,
                                 link_color_node = NULL,
                                 link_color_module = NULL,
                                 link_width_node = 1,
                                 link_width_module = 1,
                                 link_linetype_node = 2,
                                 link_linetype_module = 1,
                                 link_alpha_node = 0.25,
                                 link_alpha_module = 0.5,
                                 drop_others = FALSE,
                                 calculate_topology = FALSE,
                                 comparisons = TRUE,
                                 comparisons_groups = NULL,
                                 order = NULL,
                                 group_layout = "circle",
                                 scale_groups = TRUE,
                                 orientation = "up",
                                 angle = 0,
                                 anchor_dist = 6,
                                 layout_anchor_dist = NULL,
                                 nrow = NULL,
                                 ncol = NULL,
                                 sine_period = 4,
                                 group_label_offset = 0.2,
                                 group_label_size = 4,
                                 network_outline = FALSE,
                                 network_outline_expand = 2,
                                 network_outline_color = "grey50",
                                 network_outline_fill = NULL,
                                 network_outline_fill_alpha = 0.2,
                                 network_outline_linetype = 1,
                                 network_outline_width = 0.5,
                                 seed = 1115,
                                 # ---- deprecated (< 0.2.0) names ----
                                 ring_n = deprecated(),
                                 layout.module = deprecated(),
                                 group.by = deprecated(),
                                 fill.by = deprecated(),
                                 fill = deprecated(),
                                 color = deprecated(),
                                 pointsize = deprecated(),
                                 jitter = deprecated(),
                                 jitter_sd = deprecated(),
                                 mapping_line = deprecated(),
                                 linealpha = deprecated(),
                                 linecolor = deprecated(),
                                 inner_curve = deprecated(),
                                 inner_curvature = deprecated(),
                                 inner_curve_adaptive = deprecated(),
                                 inner_curve_adaptive_range = deprecated(),
                                 inner_curve_adaptive_bins = deprecated(),
                                 add_outer = deprecated(),
                                 q_outer = deprecated(),
                                 expand_outer = deprecated(),
                                 bandwidth_scale = deprecated(),
                                 outerwidth = deprecated(),
                                 outerlinetype = deprecated(),
                                 outeralpha = deprecated(),
                                 link_linewidth_node = deprecated(),
                                 link_linewidth_module = deprecated(),
                                 link_linealpha_node = deprecated(),
                                 link_linealpha_module = deprecated(),
                                 dropOthers = deprecated(),
                                 label_offset = deprecated(),
                                 label_size = deprecated(),
                                 add_group_outer = deprecated(),
                                 add_group_outer_expand = deprecated(),
                                 add_group_outer_color = deprecated(),
                                 add_group_outer_fill = deprecated(),
                                 add_group_outer_fill_alpha = deprecated(),
                                 add_group_outer_linetype = deprecated(),
                                 add_group_outer_linewidth = deprecated()
){
  # ---- lifecycle: translate deprecated (< 0.2.0) argument names ----------
  .fn_env <- environment()
  .old_args <- .ggnv_collect_deprecated(.fn_env, map = .ggnv_multi_link_deprecated_arg_map)
  if (length(.old_args) > 0L) {
    .explicit_new <- setdiff(names(match.call(expand.dots = FALSE))[-1L],
                             names(.ggnv_multi_link_deprecated_arg_map))
    .renamed <- .ggnv_rename_args(.old_args, fn = "ggNetView_multi_link",
                                  env = .fn_env, user_env = parent.frame(),
                                  map = .ggnv_multi_link_deprecated_arg_map)
    for (.nm in names(.renamed)) {
      if (.nm %in% .explicit_new) next
      assign(.nm, .renamed[[.nm]], envir = .fn_env)
    }
  }

  method <- match.arg(method)
  sig_by <- match.arg(sig_by)


  color_v <- c('#1f78b4', '#e7298a', '#41ab5d', '#807dba',
               '#fb9a99', '#4eb3d3', '#bc80bd', '#e31a1c',
               '#a6cee3', '#b2df8a', '#33a02c', '#fdbf6f',
               '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99',
               '#8dd3c7', '#ffffb3', '#bebada', '#fb8072',
               '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
               '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f'
               )

  if (!is.character(edge_color) || length(edge_color) != 1 || is.na(edge_color) ||
      trimws(edge_color) == "") {
    stop("`edge_color` must be a single colour or a single edge column name.")
  }
  edge_color_is_sign <- identical(edge_color, "corr_direction")
  if (!is.logical(calculate_topology) || length(calculate_topology) != 1 || is.na(calculate_topology)) {
    stop("`calculate_topology` must be TRUE or FALSE.")
  }

  use_graph_obj_list <- !is.null(graph_obj_list)
  if (use_graph_obj_list) {
    if (!is.list(graph_obj_list) || length(graph_obj_list) == 0) {
      stop("`graph_obj_list` must be a non-empty named list of tbl_graph objects.")
    }
    if (is.null(names(graph_obj_list)) || any(names(graph_obj_list) == "")) {
      stop("`graph_obj_list` must be a named list (e.g. list(WT = g1, KO = g2)).")
    }
    for (i in seq_along(graph_obj_list)) {
      if (!inherits(graph_obj_list[[i]], "tbl_graph")) {
        stop(sprintf("`graph_obj_list[[%s]]` must be a tbl_graph object.", names(graph_obj_list)[i]))
      }
    }
    valid_groups <- names(graph_obj_list)
    if (isTRUE(calculate_topology)) {
      warning("`calculate_topology` requires mat and group_info; topology will be skipped when using graph_obj_list.")
      calculate_topology <- FALSE
    }
  } else {
    if (is.null(mat) || is.null(group_info)) {
      stop("Provide either (mat, group_info) or graph_obj_list.")
    }
    if (!is.data.frame(group_info) || !all(c("Sample", "Group") %in% colnames(group_info))) {
      stop("`group_info` must be a data.frame containing columns `Sample` and `Group`.")
    }

    has_underscore <- any(grepl("_", group_info$Sample, fixed = TRUE)) ||
      any(grepl("_", group_info$Group, fixed = TRUE))
    if (has_underscore) {
      warning("`group_info` contains underscores in Sample or Group; replacing with hyphens to avoid parsing issues.")
      group_info <- group_info %>%
        dplyr::mutate(
          Sample = gsub("_", "-", as.character(Sample), fixed = TRUE),
          Group = gsub("_", "-", as.character(Group), fixed = TRUE)
        )
      mat <- as.data.frame(mat)
      if (any(grepl("_", colnames(mat), fixed = TRUE))) {
        colnames(mat) <- gsub("_", "-", colnames(mat), fixed = TRUE)
      }
      if (!is.null(order)) {
        order <- gsub("_", "-", as.character(order), fixed = TRUE)
      }
      if (!is.null(comparisons_groups) && is.list(comparisons_groups)) {
        comparisons_groups <- lapply(comparisons_groups, function(x) gsub("_", "-", as.character(x), fixed = TRUE))
      }
    }
    valid_groups <- unique(group_info$Group)
  }

  if (!is.logical(comparisons) || length(comparisons) != 1 || is.na(comparisons)) {
    stop("`comparisons` must be TRUE or FALSE.")
  }
  groups_to_process <- if (is.null(order)) {
    valid_groups
  } else {
    if (!is.character(order) || length(order) == 0) {
      stop("`order` must be NULL or a character vector of group names.")
    }
    if (!all(order %in% valid_groups)) {
      stop("All elements in `order` must be group names from group_info$Group or names(graph_obj_list).")
    }
    if (!all(valid_groups %in% order)) {
      stop("`order` must contain all unique groups from group_info$Group or names(graph_obj_list).")
    }
    if (length(order) != length(valid_groups)) {
      stop("`order` must contain each group exactly once (no duplicates, no missing).")
    }
    order
  }
  if (!is.null(node_fill_values) && !is.character(node_fill_values)) {
    stop("`node_fill_values` must be NULL or a character vector (preferably named).")
  }
  if (!is.null(node_color_values) && !is.character(node_color_values)) {
    stop("`node_color_values` must be NULL, a single color string, or a named character vector.")
  }
  if (is.logical(module_outline)) {
    if (length(module_outline) != 1 || is.na(module_outline)) {
      stop("`module_outline` must be a single logical or character string.")
    }
    module_outline <- if (isTRUE(module_outline)) "circle" else "none"
  } else if (is.character(module_outline)) {
    if (length(module_outline) != 1 || is.na(module_outline) || trimws(module_outline) == "") {
      stop("`module_outline` must be a single logical or character string.")
    }
    module_outline <- tolower(trimws(module_outline))
    if (!module_outline %in% c("circle", "manual", "none")) {
      stop("`module_outline` must be one of: 'circle', 'manual', 'none'.")
    }
  } else {
    stop("`module_outline` must be a single logical or character string.")
  }
  if (!is.logical(network_outline) || length(network_outline) != 1 || is.na(network_outline)) {
    stop("`network_outline` must be TRUE or FALSE.")
  }
  if (!is.numeric(network_outline_expand) || length(network_outline_expand) != 1 || is.na(network_outline_expand)) {
    stop("`network_outline_expand` must be a single numeric value.")
  }
  if (!is.null(network_outline_color) && !is.character(network_outline_color)) {
    stop("`network_outline_color` must be NULL or a character vector (color string(s)).")
  }
  if (!is.null(network_outline_fill) && !is.character(network_outline_fill)) {
    stop("`network_outline_fill` must be NULL or a character vector (color string(s)).")
  }
  if (!is.numeric(network_outline_fill_alpha) || length(network_outline_fill_alpha) != 1 || is.na(network_outline_fill_alpha)) {
    stop("`network_outline_fill_alpha` must be a single numeric value between 0 and 1.")
  }
  if (network_outline_fill_alpha < 0 || network_outline_fill_alpha > 1) {
    stop("`network_outline_fill_alpha` must be between 0 and 1.")
  }
  if (!is.numeric(network_outline_width) || length(network_outline_width) != 1 || is.na(network_outline_width) || network_outline_width < 0) {
    stop("`network_outline_width` must be a single non-negative numeric value.")
  }
  if (!is.character(link_level) || length(link_level) != 1 || is.na(link_level)) {
    stop("`link_level` must be one of: 'None', 'Module', 'Node', 'NodeinModule', 'Module&Node', 'Module&Node2'.")
  }
  link_level <- tolower(trimws(link_level))
  if (!link_level %in% c("none", "module", "node", "nodeinmodule", "module&node", "module&node2")) {
    stop("`link_level` must be one of: 'None', 'Module', 'Node', 'NodeinModule', 'Module&Node', 'Module&Node2'.")
  }
  if (!is.logical(link_curve) || length(link_curve) != 1 || is.na(link_curve)) {
    stop("`link_curve` must be TRUE or FALSE.")
  }
  if (!is.numeric(link_curvature) || length(link_curvature) != 1 || is.na(link_curvature)) {
    stop("`link_curvature` must be a single numeric value.")
  }
  if (!is.character(link_curve_mode) || length(link_curve_mode) != 1 || is.na(link_curve_mode)) {
    stop("`link_curve_mode` must be one of: 'outward', 'inward', 'cross'.")
  }
  link_curve_mode <- tolower(trimws(link_curve_mode))
  if (!link_curve_mode %in% c("outward", "inward", "cross")) {
    stop("`link_curve_mode` must be one of: 'outward', 'inward', 'cross'.")
  }
  if (!is.logical(link_curve_adaptive) || length(link_curve_adaptive) != 1 || is.na(link_curve_adaptive)) {
    stop("`link_curve_adaptive` must be TRUE or FALSE.")
  }
  if (!is.numeric(link_curve_adaptive_range) || length(link_curve_adaptive_range) != 2 ||
      any(is.na(link_curve_adaptive_range))) {
    stop("`link_curve_adaptive_range` must be a numeric vector of length 2.")
  }
  if (link_curve_adaptive_range[1] <= 0 || link_curve_adaptive_range[2] <= 0 ||
      link_curve_adaptive_range[1] > link_curve_adaptive_range[2]) {
    stop("`link_curve_adaptive_range` must be positive and increasing, e.g. c(0.7, 1.3).")
  }
  if (!is.numeric(link_curve_adaptive_bins) || length(link_curve_adaptive_bins) != 1 || is.na(link_curve_adaptive_bins)) {
    stop("`link_curve_adaptive_bins` must be a single integer >= 1.")
  }
  link_curve_adaptive_bins <- as.integer(link_curve_adaptive_bins)
  if (link_curve_adaptive_bins < 1) {
    stop("`link_curve_adaptive_bins` must be >= 1.")
  }
  if (!is.null(link_color_node) && !is.character(link_color_node)) {
    stop("`link_color_node` must be NULL or a character vector (color string(s)).")
  }
  if (!is.null(link_color_module) && !is.character(link_color_module)) {
    stop("`link_color_module` must be NULL or a character vector (color string(s)).")
  }
  if (!is.numeric(link_width_node) || any(is.na(link_width_node)) || any(link_width_node < 0)) {
    stop("`link_width_node` must be a numeric vector of non-negative values.")
  }
  if (!is.numeric(link_width_module) || any(is.na(link_width_module)) || any(link_width_module < 0)) {
    stop("`link_width_module` must be a numeric vector of non-negative values.")
  }
  if (!is.logical(edge_curve) || length(edge_curve) != 1 || is.na(edge_curve)) {
    stop("`edge_curve` must be TRUE or FALSE.")
  }
  if (!is.numeric(edge_curvature) || length(edge_curvature) != 1 || is.na(edge_curvature)) {
    stop("`edge_curvature` must be a single numeric value.")
  }
  if (!is.logical(edge_curve_adaptive) || length(edge_curve_adaptive) != 1 || is.na(edge_curve_adaptive)) {
    stop("`edge_curve_adaptive` must be TRUE or FALSE.")
  }
  if (!is.numeric(edge_curve_adaptive_range) || length(edge_curve_adaptive_range) != 2 ||
      any(is.na(edge_curve_adaptive_range))) {
    stop("`edge_curve_adaptive_range` must be a numeric vector of length 2.")
  }
  if (edge_curve_adaptive_range[1] <= 0 || edge_curve_adaptive_range[2] <= 0 ||
      edge_curve_adaptive_range[1] > edge_curve_adaptive_range[2]) {
    stop("`edge_curve_adaptive_range` must be positive and increasing, e.g. c(0.7, 1.3).")
  }
  if (!is.numeric(edge_curve_adaptive_bins) || length(edge_curve_adaptive_bins) != 1 || is.na(edge_curve_adaptive_bins)) {
    stop("`edge_curve_adaptive_bins` must be a single integer >= 1.")
  }
  edge_curve_adaptive_bins <- as.integer(edge_curve_adaptive_bins)
  if (edge_curve_adaptive_bins < 1) {
    stop("`edge_curve_adaptive_bins` must be >= 1.")
  }
  if (!is.numeric(anchor_dist) || length(anchor_dist) != 1 || is.na(anchor_dist)) {
    stop("`anchor_dist` must be a single numeric value.")
  }
  if (!is.null(layout_anchor_dist) &&
      (!is.numeric(layout_anchor_dist) || length(layout_anchor_dist) != 1 || is.na(layout_anchor_dist))) {
    stop("`layout_anchor_dist` must be NULL or a single numeric value.")
  }
  if (!is.numeric(group_label_offset) || length(group_label_offset) != 1 || is.na(group_label_offset)) {
    stop("`group_label_offset` must be a single numeric value.")
  }
  if (!is.numeric(group_label_size) || length(group_label_size) != 1 || is.na(group_label_size) || group_label_size <= 0) {
    stop("`group_label_size` must be a single positive numeric value.")
  }
  if (!is.numeric(sine_period) || length(sine_period) != 1 || is.na(sine_period) || sine_period <= 0) {
    stop("`sine_period` must be a single positive numeric value.")
  }
  layout_module <- match.arg(layout_module)
  group_layout <- match.arg(group_layout, choices = c("circle", "row", "column", "square", "diamond", "triangle", "triangle_down", "snake", "snake_vertical", "snake_vertical_sin", "snake_vertical_cos", "snake_vertical_neg_sin", "snake_vertical_neg_cos", "sin", "cos", "-sin", "-cos", "center_pairs"))
  layout_anchor_dist_use <- if (is.null(layout_anchor_dist)) anchor_dist else layout_anchor_dist
  if (is.null(layout) || length(layout) != 1 || is.na(layout) || trimws(as.character(layout)) == "") {
    stop("`layout` must be a non-empty character string (e.g. 'gephi', 'square', 'petal').")
  }
  layout <- as.character(layout)

  graph_list <- list()

  graph_info <- list()

  # Node -> module lookup captured from the FULL graph (before `drop_others`
  # removes the "Others" nodes). The module-overlap test must always run on the
  # complete network so that `drop_others` stays a display-only switch and does
  # not silently re-run the hypergeometric test on a shrunken universe.
  module_ref <- list()

  graph_stat <- list()
  topology_network <- list()
  topology_sample <- list()

  for (g in groups_to_process) {
    if (use_graph_obj_list) {
      graph <- graph_obj_list[[g]]
    } else {
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
    }

    # Snapshot name/module from the complete graph before anything is dropped;
    # `Module_information` below is computed from this, never from the possibly
    # Others-free plotting graph.
    full_node_tbl <- graph %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
    module_ref[[g]] <- full_node_tbl

    # drop_others acts on the source graph_obj BEFORE layout:
    # it removes "Others" nodes first, then downstream layout/plot are rebuilt.
    if (isTRUE(drop_others)) {
      node_tbl <- graph %>%
        tidygraph::activate(nodes) %>%
        tidygraph::as_tibble()

      module_candidates <- c("Modularity", "modularity3", "modularity2")
      module_col <- module_candidates[module_candidates %in% colnames(node_tbl)]
      module_col <- if (length(module_col) > 0) module_col[[1]] else NULL

      if (!is.null(module_col)) {
        if ("name" %in% colnames(node_tbl)) {
          keep_names <- node_tbl %>%
            dplyr::filter(as.character(.data[[module_col]]) != "Others") %>%
            dplyr::pull(name)

          graph <- graph %>%
            tidygraph::activate(nodes) %>%
            tidygraph::filter(name %in% keep_names)
        } else {
          graph <- graph %>%
            tidygraph::activate(nodes) %>%
            tidygraph::filter(as.character(.data[[module_col]]) != "Others")
        }
      } else {
        warning("`drop_others = TRUE` but no module column found in `graph_obj` nodes.")
      }
    }


    # find layout function
    func_name <- paste0("create_layout_", layout)

    # find layout functions from ggNetView package
    lay_func <- utils::getFromNamespace(func_name, "ggNetView")

    # get ly1 (auto-pass only arguments supported by the chosen layout function)
    # `inner_shrink` is currently consumed only by create_layout_WGCNA; the
    # auto-filter below silently drops it for any other layout.
    lay_args <- list(
      graph_obj = graph,
      node_add = node_add,
      r = r,
      inner_shrink = inner_shrink,
      scale = scale_groups,
      anchor_dist = layout_anchor_dist_use,
      orientation = orientation,
      angle = angle,
      nrow = nrow,
      ncol = ncol
    )
    lay_args <- lay_args[names(lay_args) %in% names(formals(lay_func))]
    ly1 <- do.call(lay_func, lay_args)

    # get ly1_1

    if (layout_module == "random") {
      ly1_1 <- module_layout(graph,
                             layout = ly1,
                             center = center,
                             idx = idx,
                             shrink = shrink,
                             node_jitter,
                             node_jitter_sd# ,
                             # seed = seed
      )
    }

    if (layout_module == "adjacent") {
      k_nn_try <- k_nn
      k_nn_cap <- max(1, nrow(ly1) - 1)
      ly1_1 <- NULL
      while (is.null(ly1_1)) {
        ly_try <- tryCatch(
          module_layout3(graph,
                         layout = ly1,
                         center = center,
                         k_nn = k_nn_try,
                         push_others_delta = push_others_delta,
                         shrink = shrink,
                         jitter = node_jitter,
                         jitter_sd = node_jitter_sd
                         # seed = seed
          ),
          error = function(e) e
        )

        if (!inherits(ly_try, "error")) {
          ly1_1 <- ly_try
          break
        }

        err_msg <- conditionMessage(ly_try)
        is_slot_error <- grepl("consecutive slots", err_msg, ignore.case = TRUE)
        if (!is_slot_error || k_nn_try >= k_nn_cap) {
          stop(ly_try)
        }

        next_k <- min(k_nn_cap, max(k_nn_try + 20, ceiling(k_nn_try * 1.25)))
        if (next_k <= k_nn_try) {
          stop(ly_try)
        }
        warning(sprintf(
          "`layout_module = 'adjacent'` failed at k_nn = %d; retrying with k_nn = %d.",
          k_nn_try, next_k
        ))
        k_nn_try <- next_k
      }
    }

    if (layout_module == "order" && func_name != "create_layout_multirings") {
      ly1_1 <- module_layout4(graph,
                              layout = ly1,
                              center = center,
                              k_nn = k_nn,
                              push_others_delta = push_others_delta,
                              shrink = shrink,
                              node_jitter,
                              node_jitter_sd
                              # seed = seed
      )
    }

    if (layout_module == "order" & func_name == "create_layout_multirings") {
      ly1_1 <- module_layout5(graph,
                              layout = ly1,
                              center = center,
                              k_nn = k_nn,
                              push_others_delta = push_others_delta,
                              shrink = shrink,
                              node_jitter,
                              node_jitter_sd
                              # seed = seed
      )
    }

    graph_list[[g]] <- graph

    graph_info[[g]] <- ly1_1$ggplot_data

    graph_stat[[g]] <- stat_graph(graph, edge_color_is_sign)

    if (isTRUE(calculate_topology)) {
      topology_network[[g]] <- tryCatch(
        get_network_topology_parallel(
          graph_obj = graph,
          mat = mat_sub,
          transfrom.method = transfrom.method,
          r.threshold = r.threshold,
          p.threshold = p.threshold,
          method = method,
          cor.method = cor.method,
          proc = proc,
          SpiecEasi.method = SpiecEasi.method,
          sparcc_R = sparcc_R,
          seed = seed
        ),
        error = function(e) {
          warning(paste0("Topology (network) failed for group `", g, "`: ", e$message), call. = FALSE)
          NULL
        }
      )

      topology_sample[[g]] <- tryCatch(
        get_sample_subgraph_topology_parallel(
          graph_obj = graph,
          mat = mat_sub,
          transfrom.method = transfrom.method,
          r.threshold = r.threshold,
          p.threshold = p.threshold,
          method = method,
          cor.method = cor.method,
          proc = proc,
          SpiecEasi.method = SpiecEasi.method,
          sparcc_R = sparcc_R,
          seed = seed
        ),
        error = function(e) {
          warning(paste0("Topology (sample subgraph) failed for group `", g, "`: ", e$message), call. = FALSE)
          NULL
        }
      )
    }
  }

  graph_list_length <- length(graph_list)

  if (!isTRUE(comparisons)) {
    compare_matrix <- matrix(character(0), nrow = 2, ncol = 0)
  } else if (graph_list_length >= 2) {
    compare_matrix <- utils::combn(names(graph_list), 2)
    if (!is.null(comparisons_groups)) {
      if (!is.list(comparisons_groups)) {
        stop("`comparisons_groups` must be a list, e.g. list(c('WT', 'OE'), c('WT', 'KO')).")
      }
      valid_groups <- names(graph_list)
      comp_cols <- list()
      seen_pairs <- character(0)
      for (comp in comparisons_groups) {
        if (!is.character(comp) || length(comp) != 2) {
          stop("Each element of `comparisons_groups` must be a character vector of length 2 (e.g. c('WT', 'OE')).")
        }
        if (!all(comp %in% valid_groups)) {
          warning(sprintf("Skipping comparison c('%s', '%s'): group(s) not in data.", comp[1], comp[2]))
          next
        }
        pair_key <- paste(sort(comp), collapse = "_")
        if (pair_key %in% seen_pairs) {
          warning(sprintf("Skipping duplicate comparison c('%s', '%s').", comp[1], comp[2]))
          next
        }
        seen_pairs <- c(seen_pairs, pair_key)
        comp_cols[[length(comp_cols) + 1]] <- matrix(comp, nrow = 2)
      }
      if (length(comp_cols) > 0) {
        compare_matrix <- do.call(cbind, comp_cols)
      } else {
        compare_matrix <- matrix(character(0), nrow = 2, ncol = 0)
      }
    }
  } else {
    compare_matrix <- matrix(character(0), nrow = 2, ncol = 0)
  }

  if (link_level %in% c("module", "nodeinmodule", "module&node", "module&node2") && ncol(compare_matrix) > 0) {
    compare_out_list <- list()
    for (i in seq_len(ncol(compare_matrix))) {
      # always compare on the full network (see `module_ref` above)
      tmp <- compare_modules_by_overlap(module_ref[[compare_matrix[1, i]]] %>%
                                          dplyr::select(name, Modularity) %>%
                                          dplyr::mutate(Group = compare_matrix[1, i]),
                                        module_ref[[compare_matrix[2, i]]] %>%
                                          dplyr::select(name, Modularity) %>%
                                          dplyr::mutate(Group = compare_matrix[2, i])) %>%
        dplyr::mutate(Group = stringr::str_c(compare_matrix[1, i],
                                             "to",
                                             compare_matrix[2, i],
                                             sep = "_"),
                      # set GroupA/GroupB directly from the comparison pair
                      # instead of splitting `Group` back apart on "_" -- group
                      # names supplied via `graph_obj_list` may themselves
                      # contain underscores, which made tidyr::separate() mis-parse
                      # them into the wrong GroupA/GroupB.
                      GroupA = compare_matrix[1, i],
                      GroupB = compare_matrix[2, i])

      compare_out_list[[stringr::str_c(compare_matrix[1, i],
                                       "to",
                                       compare_matrix[2, i],
                                       sep = "_")]] <- tmp
    }

    # significance filter: `sig_by = "pvalue"` (default, unchanged behaviour) uses
    # the raw p-value; `sig_by = "FDR"` uses the BH-adjusted column from
    # compare_modules_info(), which is the statistically appropriate choice when
    # many module pairs are tested.
    Module_information <- do.call(rbind, compare_out_list) %>%
      dplyr::filter(.data[[sig_by]] < 0.05)
  } else {
    Module_information <- tibble::tibble(
      modA = character(),
      modB = character(),
      overlap = integer(),
      sizeA = integer(),
      sizeB = integer(),
      overlap_coef = numeric(),
      pvalue = numeric(),
      FDR = numeric(),
      Group = character(),
      GroupA = character(),
      GroupB = character()
    )
  }

  Module_information_plot <- if (link_level %in% c("module", "nodeinmodule", "module&node", "module&node2")) {
    # "Others" is a display bucket (everything below `top_modules`), not a real
    # community, so it never takes part in cross-group module links -- regardless
    # of `drop_others`. Dropping a pair needs BOTH sides to be real modules.
    Module_information %>%
      dplyr::filter(modA != "Others", modB != "Others")
  } else {
    Module_information
  }

  # if (isTRUE(scale_groups)) {
  #   anchor_dist = 1
  # }else{
  #   anchor_dist = anchor_dist * 30
  # }


  n_grp <- graph_list_length
  if (group_layout == "circle") {
    angles <- pi/2 - 2 * pi * (0:(n_grp - 1)) / n_grp
    anchors <- lapply(angles, function(a) {
      c(anchor_dist * cos(a), anchor_dist * sin(a))
    })
  } else if (group_layout == "row") {

    nr <- if (!is.null(nrow)) as.integer(nrow) else 1L
    nc <- if (!is.null(ncol)) as.integer(ncol) else n_grp
    if (is.null(nrow) && !is.null(ncol)) nr <- max(1L, as.integer(ceiling(n_grp / nc)))
    if (!is.null(nrow) && is.null(ncol)) nc <- max(1L, as.integer(ceiling(n_grp / nr)))
    nr <- max(1L, nr)
    nc <- max(1L, nc)
    anchors <- lapply(seq_len(n_grp), function(i) {
      ii <- i - 1L
      r <- ii %/% nc
      c <- ii %% nc
      x <- (c - (nc - 1) / 2) * anchor_dist
      y <- -((r - (nr - 1) / 2) * anchor_dist)
      c(x, y)
    })
  } else if (group_layout == "column") {

    nr <- if (!is.null(nrow)) as.integer(nrow) else n_grp
    nc <- if (!is.null(ncol)) as.integer(ncol) else 1L
    if (is.null(nrow) && !is.null(ncol)) nr <- max(1L, as.integer(ceiling(n_grp / nc)))
    if (!is.null(nrow) && is.null(ncol)) nc <- max(1L, as.integer(ceiling(n_grp / nr)))
    nr <- max(1L, nr)
    nc <- max(1L, nc)
    anchors <- lapply(seq_len(n_grp), function(i) {
      ii <- i - 1L
      r <- ii %% nr
      c <- ii %/% nr
      x <- (c - (nc - 1) / 2) * anchor_dist
      y <- -((r - (nr - 1) / 2) * anchor_dist)
      c(x, y)
    })
  } else if (group_layout == "snake") {

    nr <- if (!is.null(nrow)) as.integer(nrow) else 1L
    nc <- if (!is.null(ncol)) as.integer(ncol) else n_grp
    if (is.null(nrow) && !is.null(ncol)) nr <- max(1L, as.integer(ceiling(n_grp / nc)))
    if (!is.null(nrow) && is.null(ncol)) nc <- max(1L, as.integer(ceiling(n_grp / nr)))
    nr <- max(1L, nr)
    nc <- max(1L, nc)
    anchors <- lapply(seq_len(n_grp), function(i) {
      ii <- i - 1L
      r <- ii %/% nc
      c <- if (r %% 2L == 0L) ii %% nc else (nc - 1L) - (ii %% nc)
      x <- (c - (nc - 1) / 2) * anchor_dist
      y <- -((r - (nr - 1) / 2) * anchor_dist)
      c(x, y)
    })
  } else if (group_layout %in% c("snake_vertical", "snake_vertical_sin", "snake_vertical_cos", "snake_vertical_neg_sin", "snake_vertical_neg_cos")) {

    t_vals <- 2 * pi * (0:(n_grp - 1)) / sine_period
    x_centered <- (0:(n_grp - 1) - (n_grp - 1) / 2) * anchor_dist
    layout_curve <- if (group_layout == "snake_vertical_sin") "snake_vertical" else group_layout
    y_fun <- switch(layout_curve,
      snake_vertical         = sin,
      snake_vertical_cos     = cos,
      snake_vertical_neg_sin = function(x) -sin(x),
      snake_vertical_neg_cos = function(x) -cos(x)
    )
    y_vals <- anchor_dist * y_fun(t_vals)
    anchors <- lapply(seq_len(n_grp), function(i) c(x_centered[i], y_vals[i]))
  } else if (group_layout %in% c("sin", "cos", "-sin", "-cos")) {


    x_centered <- (0:(n_grp - 1) - (n_grp - 1) / 2) * anchor_dist
    y_vals <- switch(group_layout,
      sin = c(0, anchor_dist * rep(c(1, -1), length.out = n_grp - 1)),
      cos = anchor_dist * rep(c(1, -1), length.out = n_grp),
      "-sin" = c(0, anchor_dist * rep(c(-1, 1), length.out = n_grp - 1)),
      "-cos" = anchor_dist * rep(c(-1, 1), length.out = n_grp)
    )
    anchors <- lapply(seq_len(n_grp), function(i) c(x_centered[i], y_vals[i]))
  } else if (group_layout == "center_pairs") {


    anchors <- vector("list", n_grp)
    anchors[[1]] <- c(0, 0)
    if (n_grp >= 2) {
      k <- (n_grp - 1) %/% 2
      x_vals <- seq_len(k) * anchor_dist
      pair_idx <- 0
      i <- 2
      while (i <= n_grp) {
        pair_idx <- pair_idx + 1
        x_pos <- x_vals[pair_idx]
        if (i + 1 <= n_grp) {
          anchors[[i]] <- c(x_pos, anchor_dist)
          anchors[[i + 1]] <- c(x_pos, -anchor_dist)
          i <- i + 2
        } else {
          anchors[[i]] <- c(x_pos, anchor_dist)
          i <- i + 1
        }
      }
    }
  } else if (group_layout %in% c("square", "diamond", "triangle", "triangle_down")) {





    base_angle <- switch(group_layout,
      square        = 3 * pi / 4,
      diamond       = pi / 2,
      triangle      = pi / 2,
      triangle_down = -pi / 2
    )
    angles <- base_angle - 2 * pi * (0:(n_grp - 1)) / n_grp
    anchors <- lapply(angles, function(a) {
      c(anchor_dist * cos(a), anchor_dist * sin(a))
    })
  } else {
    stop("`group_layout` must be one of: 'circle', 'row', 'column', 'square', 'diamond', 'triangle', 'triangle_down', 'snake', 'snake_vertical', 'snake_vertical_sin', 'snake_vertical_cos', 'snake_vertical_neg_sin', 'snake_vertical_neg_cos', 'sin', 'cos', '-sin', '-cos', 'center_pairs'.")
  }

  anchors_df <- do.call(rbind, anchors) %>%
    as.data.frame() %>%
    purrr::set_names(c("x", "y"))



  if (isTRUE(scale_groups)) {

    for (i in seq_along(names(graph_info))) {
      graph_info[[i]]$ggplot_node_df <- graph_info[[i]]$ggplot_node_df %>%
        dplyr::mutate(Group = names(graph_info)[i]) %>%
        # dplyr::select(name, x, y, Group) %>%
        dplyr::mutate(
          ymin = min(y),
          ymax = max(y),
          xmin = min(x),
          xmax = max(x),
          xmind = (max(x) + min(x)) / 2,
          ymind = (max(y) + min(y)) / 2,
          scale_v = max(xmax - xmin, ymax - ymin),
          x = (x - xmind)/scale_v,
          y = (y - ymind)/scale_v
        ) %>%
        dplyr::mutate(
          x = x + anchors_df[i,1],
          y = y + anchors_df[i,2]
        )

      graph_info[[i]]$ggplot_edge_df <- graph_info[[i]]$ggplot_edge_df %>%
        dplyr::mutate(xmid = unique(graph_info[[i]]$ggplot_node_df$xmind),
                      ymid = unique(graph_info[[i]]$ggplot_node_df$ymind),
                      scale_v = unique(graph_info[[i]]$ggplot_node_df$scale_v),
                      from_x = (from_x - xmid)/scale_v,
                      from_y = (from_y - ymid)/scale_v,
                      to_x = (to_x - xmid)/scale_v,
                      to_y = (to_y - ymid)/scale_v) %>%
        dplyr::mutate(
          from_x = from_x + anchors_df[i,1],
          from_y = from_y + anchors_df[i,2],
          to_x = to_x + anchors_df[i,1],
          to_y = to_y + anchors_df[i,2]
        )

    }
    # node_jitter TRUE
    if (isTRUE(node_jitter)) {
      for (i in seq_along(names(graph_info))) {
        graph_info[[i]]$ggplot_node_df <- graph_info[[i]]$ggplot_node_df %>%
          dplyr::mutate(
            x = x + stats::rnorm(dplyr::n(), mean = 0, sd = node_jitter_sd),
            y = y + stats::rnorm(dplyr::n(), mean = 0, sd = node_jitter_sd)
          )

        graph_info[[i]]$ggplot_edge_df <- graph_info[[i]]$ggplot_edge_df %>%
          dplyr::select(from, to, weight, correlation, corr_direction, from_id, to_id) %>%
          dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                             dplyr::select(name, x, y),
                           by = c("from_id" = "name")) %>%
          dplyr::rename(from_x = x,
                        from_y = y) %>%
          dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                             dplyr::select(name, x, y),
                           by = c("to_id" = "name")) %>%
          dplyr::rename(to_x = x,
                        to_y = y)
      }
    }
  }

  if (!isTRUE(scale_groups)) {
    for (i in seq_along(names(graph_info))) {
      graph_info[[i]]$ggplot_node_df <- graph_info[[i]]$ggplot_node_df %>%
        dplyr::mutate(
          x = x + anchors_df[i,1],
          y = y + anchors_df[i,2]
        )

      graph_info[[i]]$ggplot_edge_df <- graph_info[[i]]$ggplot_edge_df %>%
        dplyr::select(from, to, weight, correlation, corr_direction, from_id, to_id) %>%
        dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                           dplyr::select(name, x, y),
                         by = c("from_id" = "name")) %>%
        dplyr::rename(from_x = x,
                      from_y = y) %>%
        dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                           dplyr::select(name, x, y),
                         by = c("to_id" = "name")) %>%
        dplyr::rename(to_x = x,
                      to_y = y)

    # node_jitter TRUE
    if (isTRUE(node_jitter)) {
      for (i in seq_along(names(graph_info))) {
        graph_info[[i]]$ggplot_node_df <- graph_info[[i]]$ggplot_node_df %>%
          dplyr::mutate(
            x = x + stats::rnorm(dplyr::n(), mean = 0, sd = node_jitter_sd),
            y = y + stats::rnorm(dplyr::n(), mean = 0, sd = node_jitter_sd)
          )

        graph_info[[i]]$ggplot_edge_df <- graph_info[[i]]$ggplot_edge_df %>%
          dplyr::select(from, to, weight, correlation, corr_direction, from_id, to_id) %>%
          dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                             dplyr::select(name, x, y),
                           by = c("from_id" = "name")) %>%
          dplyr::rename(from_x = x,
                        from_y = y) %>%
          dplyr::left_join(graph_info[[i]]$ggplot_node_df %>%
                             dplyr::select(name, x, y),
                           by = c("to_id" = "name")) %>%
          dplyr::rename(to_x = x,
                        to_y = y)
       }
     }
    }
  }




  # plot data

  # raw data to plot data

  p <- ggplot()

  add_inner_link_geom <- function(edge_data, mapping_obj, color_fixed = NULL) {
    if (!isTRUE(edge_curve)) {
      return(
        ggplot2::geom_segment(
          data = edge_data,
          mapping = mapping_obj,
          alpha = edge_alpha,
          color = color_fixed
        )
      )
    }

    if (!isTRUE(edge_curve_adaptive) || nrow(edge_data) < 2) {
      return(
        ggplot2::geom_curve(
          data = edge_data,
          mapping = mapping_obj,
          curvature = edge_curvature,
          alpha = edge_alpha,
          color = color_fixed
        )
      )
    }

    edge_dist <- sqrt((edge_data$to_x - edge_data$from_x)^2 + (edge_data$to_y - edge_data$from_y)^2)
    if (all(!is.finite(edge_dist)) || diff(range(edge_dist, na.rm = TRUE)) <= .Machine$double.eps) {
      return(
        ggplot2::geom_curve(
          data = edge_data,
          mapping = mapping_obj,
          curvature = edge_curvature,
          alpha = edge_alpha,
          color = color_fixed
        )
      )
    }

    n_bins <- min(edge_curve_adaptive_bins, length(edge_dist))
    breaks <- unique(stats::quantile(edge_dist, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    if (length(breaks) <= 2) {
      return(
        ggplot2::geom_curve(
          data = edge_data,
          mapping = mapping_obj,
          curvature = edge_curvature,
          alpha = edge_alpha,
          color = color_fixed
        )
      )
    }

    edge_data <- edge_data %>%
      dplyr::mutate(.edge_dist = edge_dist,
                    .curve_bin = as.integer(cut(.edge_dist, breaks = breaks, include.lowest = TRUE)))

    bins_present <- sort(unique(stats::na.omit(edge_data$.curve_bin)))
    if (length(bins_present) == 0) {
      return(
        ggplot2::geom_curve(
          data = edge_data,
          mapping = mapping_obj,
          curvature = edge_curvature,
          alpha = edge_alpha,
          color = color_fixed
        )
      )
    }

    bin_mult <- stats::setNames(
      seq(edge_curve_adaptive_range[1], edge_curve_adaptive_range[2], length.out = length(bins_present)),
      bins_present
    )

    layers <- lapply(bins_present, function(b) {
      ggplot2::geom_curve(
        data = edge_data %>% dplyr::filter(.curve_bin == b),
        mapping = mapping_obj,
        curvature = edge_curvature * bin_mult[as.character(b)],
        alpha = edge_alpha,
        color = color_fixed
      )
    })

    return(layers)
  }

  # link_info, cross_channel_lookup, add_link_layer (needed for node links before nodes)
  link_info <- tibble::tibble(
    link_level = character(),
    group_a = character(),
    group_b = character(),
    source = character(),
    target = character(),
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric()
  )

  cross_channel_lookup <- list()
  circle_centers_x <- circle_centers_y <- circle_radius <- NULL
  circle_pad <- 0
  module_circle_tbl <- NULL
  module_circle_pad <- 0
  if (identical(link_curve_mode, "cross")) {
    gnames <- names(graph_info)
    ngrp <- length(gnames)
    cx <- vapply(gnames, function(g) mean(graph_info[[g]]$ggplot_node_df$x), numeric(1))
    cy <- vapply(gnames, function(g) mean(graph_info[[g]]$ggplot_node_df$y), numeric(1))
    names(cx) <- gnames
    names(cy) <- gnames
    circle_centers_x <- cx
    circle_centers_y <- cy
    circle_radius <- vapply(
      gnames,
      function(g) {
        node_df <- graph_info[[g]]$ggplot_node_df
        d <- sqrt((node_df$x - cx[g])^2 + (node_df$y - cy[g])^2)
        module_outline_q <- as.numeric(stats::quantile(d, probs = 0.98, na.rm = TRUE))
        r_max <- max(d, na.rm = TRUE)
        max(max(module_outline_q, r_max * 0.98), 1e-6)
      },
      numeric(1)
    )
    circle_pad <- max(stats::median(circle_radius, na.rm = TRUE) * 0.06, 1e-6)
    module_rows <- list()
    rid <- 1L
    for (g in gnames) {
      node_df <- graph_info[[g]]$ggplot_node_df
      if (!"Modularity" %in% colnames(node_df)) next
      mods <- unique(as.character(node_df$Modularity))
      mods <- mods[!is.na(mods)]
      for (m in mods) {
        md <- node_df[as.character(node_df$Modularity) == m, , drop = FALSE]
        if (nrow(md) == 0) next
        mcx <- mean(md$x)
        mcy <- mean(md$y)
        mr <- max(sqrt((md$x - mcx)^2 + (md$y - mcy)^2))
        if (!is.finite(mr) || mr <= 0) next
        module_rows[[rid]] <- data.frame(
          group = g,
          module = m,
          key = paste0(g, "|", m),
          cx = mcx,
          cy = mcy,
          r = mr,
          stringsAsFactors = FALSE
        )
        rid <- rid + 1L
      }
    }
    if (length(module_rows) > 0) {
      module_circle_tbl <- dplyr::bind_rows(module_rows)
      module_circle_pad <- max(stats::median(module_circle_tbl$r, na.rm = TRUE) * 0.10, 1e-6)
    }
    if (ngrp >= 2L) {
      if (ngrp == 3L) {
        top_ix <- which.max(cy)
        bottom_ix <- setdiff(seq_len(3L), top_ix)
        left_ix <- bottom_ix[which.min(cx[bottom_ix])]
        right_ix <- setdiff(bottom_ix, left_ix)
        top_g <- gnames[top_ix]
        left_g <- gnames[left_ix]
        right_g <- gnames[right_ix]
        outward_tl <- (cx[top_g] + cx[left_g])/2 * (cy[top_g] - cy[left_g]) + (cy[top_g] + cy[left_g])/2 * (cx[left_g] - cx[top_g])
        cross_channel_lookup[[paste(sort(c(top_g, left_g)), collapse = "|")]] <- as.integer(if (outward_tl != 0) -sign(outward_tl) else -1L)
        outward_tr <- (cx[top_g] + cx[right_g])/2 * (cy[top_g] - cy[right_g]) + (cy[top_g] + cy[right_g])/2 * (cx[right_g] - cx[top_g])
        cross_channel_lookup[[paste(sort(c(top_g, right_g)), collapse = "|")]] <- as.integer(if (outward_tr != 0) -sign(outward_tr) else 1L)
        outward_lr <- (cx[left_g] + cx[right_g])/2 * (cy[left_g] - cy[right_g]) + (cy[left_g] + cy[right_g])/2 * (cx[right_g] - cx[left_g])
        cross_channel_lookup[[paste(sort(c(left_g, right_g)), collapse = "|")]] <- as.integer(if (outward_lr != 0) -sign(outward_lr) else 1L)
      } else {
        cx_cent <- mean(cx)
        cy_cent <- mean(cy)
        for (i in seq_len(ngrp - 1L)) {
          for (j in (i + 1L):ngrp) {
            gA <- gnames[i]
            gB <- gnames[j]
            dx <- cx[gB] - cx[gA]
            dy <- cy[gB] - cy[gA]
            mx <- (cx[gA] + cx[gB]) / 2
            my <- (cy[gA] + cy[gB]) / 2
            s <- as.integer(sign((cx_cent - mx) * dy - (cy_cent - my) * dx))
            key_ij <- paste(sort(c(gA, gB)), collapse = "|")
            cross_channel_lookup[[key_ij]] <- if (is.finite(s) && s != 0L) s else 1L
          }
        }
      }
    }
  }


  resolve_link_color <- function(color_vec, default_palette, pair_key, pair_idx) {
    if (is.null(color_vec) || length(color_vec) == 0L) {
      return(default_palette[((pair_idx - 1L) %% length(default_palette)) + 1L])
    }
    if (length(color_vec) == 1L) return(color_vec[1L])
    if (!is.null(names(color_vec)) && pair_key %in% names(color_vec)) return(color_vec[pair_key])
    return(color_vec[((pair_idx - 1L) %% length(color_vec)) + 1L])
  }

  resolve_link_style <- function(val, pair_key, pair_idx) {
    if (is.null(val) || length(val) == 0L) return(NULL)
    if (length(val) == 1L) return(val[1L])
    if (!is.null(names(val)) && pair_key %in% names(val)) return(val[pair_key])
    return(val[((pair_idx - 1L) %% length(val)) + 1L])
  }

  # add_link_layer: complex version with group + module circle avoidance for node links.
  # Uses cross_channel_lookup, circle_centers, circle_radius, module_circle_tbl from above.
  add_link_layer <- (function() {
    add_link_layer_impl <- function(p_obj, df_link, col_link, link_type = c("node", "module"), cross_channel_sign = NULL, link_lt = NULL, link_al = NULL, link_lw = NULL) {
      link_type <- match.arg(link_type)
      link_lt <- if (!is.null(link_lt)) link_lt else (if (link_type == "module") link_linetype_module else link_linetype_node)
      link_al <- if (!is.null(link_al)) link_al else (if (link_type == "module") link_alpha_module else link_alpha_node)
      link_lw <- if (!is.null(link_lw)) link_lw else (if (link_type == "module") link_width_module else link_width_node)
      if (nrow(df_link) == 0) return(p_obj)
      if (isTRUE(link_curve)) {
        df_draw <- df_link
        dx <- df_draw$xend - df_draw$x
        dy <- df_draw$yend - df_draw$y
        mx <- (df_draw$x + df_draw$xend) / 2
        my <- (df_draw$y + df_draw$yend) / 2
        dist_link <- sqrt(dx^2 + dy^2)

        base_abs_curv <- abs(link_curvature)
        curv_abs <- rep(base_abs_curv, nrow(df_draw))
        curv_sign <- rep(1L, nrow(df_draw))
        if (identical(link_curve_mode, "cross")) {
          preferred_sign <- if (!is.null(cross_channel_sign) && length(cross_channel_sign) == 1 && is.finite(cross_channel_sign)) {
            as.integer(ifelse(cross_channel_sign >= 0, 1L, -1L))
          } else {
            NA_integer_
          }
          curv_sign <- rep(ifelse(is.na(preferred_sign), 1L, preferred_sign), nrow(df_draw))
          curv_abs <- rep(max(base_abs_curv, .Machine$double.eps), nrow(df_draw))
          if (!is.null(circle_radius) && length(circle_radius) > 0) {
            t_seq <- seq(0.05, 0.95, length.out = 81)
            endpoint_allow <- 0.26
            ctrl_scale <- 0.42
            candidate_abs <- max(base_abs_curv, 0.02) * c(0.10, 0.14, 0.20, 0.27, 0.36, 0.48, 0.64, 0.84, 1.1, 1.45, 1.9, 2.5, 3.1, 4.0, 5.2, 6.6)
            path_inside_count <- function(i, k_signed) {
              x0 <- df_draw$x[i]
              y0 <- df_draw$y[i]
              x1 <- df_draw$xend[i]
              y1 <- df_draw$yend[i]
              dd <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
              if (!is.finite(dd) || dd <= .Machine$double.eps) return(0L)
              midx <- (x0 + x1) / 2
              midy <- (y0 + y1) / 2
              ux <- -(y1 - y0) / dd
              uy <- (x1 - x0) / dd
              cx_ctrl <- midx + ux * (k_signed * dd * ctrl_scale)
              cy_ctrl <- midy + uy * (k_signed * dd * ctrl_scale)
              bx <- (1 - t_seq)^2 * x0 + 2 * (1 - t_seq) * t_seq * cx_ctrl + t_seq^2 * x1
              by <- (1 - t_seq)^2 * y0 + 2 * (1 - t_seq) * t_seq * cy_ctrl + t_seq^2 * y1
              gsrc <- gtgt <- NA_character_
              if (!is.null(df_draw$group_a) && !is.null(df_draw$group_b)) {
                gsrc <- as.character(df_draw$group_a[i])
                gtgt <- as.character(df_draw$group_b[i])
              }
              inside_any <- rep(FALSE, length(t_seq))
              for (g in names(circle_radius)) {
                rr <- (circle_radius[[g]] + circle_pad)^2
                dist2 <- (bx - circle_centers_x[[g]])^2 + (by - circle_centers_y[[g]])^2
                inside <- dist2 < rr
                if (!is.na(gsrc) && identical(g, gsrc)) inside <- inside & (t_seq > endpoint_allow)
                if (!is.na(gtgt) && identical(g, gtgt)) inside <- inside & (t_seq < (1 - endpoint_allow))
                inside_any <- inside_any | inside
              }
              sum(inside_any)
            }
            for (ii in seq_len(nrow(df_draw))) {
              found <- FALSE
              best_hit <- Inf
              best_sign <- ifelse(is.na(preferred_sign), 1L, preferred_sign)
              best_abs <- max(base_abs_curv, .Machine$double.eps)
              for (cand_abs in candidate_abs) {
                sign_trials <- if (cand_abs <= .Machine$double.eps) c(ifelse(is.na(preferred_sign), 1L, preferred_sign)) else
                  if (is.na(preferred_sign)) c(-1L, 1L) else c(preferred_sign, -preferred_sign)
                for (sgn_i in sign_trials) {
                  hit_n <- path_inside_count(ii, sgn_i * cand_abs)
                  if (hit_n < best_hit) { best_hit <- hit_n; best_sign <- sgn_i; best_abs <- cand_abs }
                  if (hit_n == 0L) {
                    curv_sign[ii] <- sgn_i
                    curv_abs[ii] <- max(cand_abs, .Machine$double.eps)
                    found <- TRUE
                    break
                  }
                }
                if (found) break
              }
              if (!found) { curv_sign[ii] <- best_sign; curv_abs[ii] <- max(best_abs, .Machine$double.eps) }
            }
          }
        } else if (identical(link_curve_mode, "inward")) {
          outward_score <- mx * (-dy) + my * dx
          curv_sign <- ifelse(outward_score >= 0, 1, -1)
        } else {
          outward_score <- mx * (-dy) + my * dx
          curv_sign <- ifelse(outward_score >= 0, -1, 1)
        }

        if (identical(link_curve_mode, "cross") && !is.null(circle_radius) && length(circle_radius) > 0) {
          ctrl_scale <- 0.42
          lerp2 <- function(p, q, t) p + (q - p) * t
          split_quad <- function(P0, P1, P2, t) {
            P01 <- lerp2(P0, P1, t); P12 <- lerp2(P1, P2, t); P0112 <- lerp2(P01, P12, t)
            list(left = list(P0, P01, P0112), right = list(P0112, P12, P2))
          }
          bezier_ctrl <- vector("list", nrow(df_draw))
          for (ii in seq_len(nrow(df_draw))) {
            x0 <- df_draw$x[ii]; y0 <- df_draw$y[ii]; x1 <- df_draw$xend[ii]; y1 <- df_draw$yend[ii]
            dd <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
            if (!is.finite(dd) || dd <= .Machine$double.eps) next
            midx <- (x0 + x1) / 2; midy <- (y0 + y1) / 2
            ux <- -(y1 - y0) / dd; uy <- (x1 - x0) / dd
            k_signed <- curv_sign[ii] * curv_abs[ii]
            P0 <- c(x0, y0); P1 <- c(midx + ux * (k_signed * dd * ctrl_scale), midy + uy * (k_signed * dd * ctrl_scale)); P2 <- c(x1, y1)
            sp1 <- split_quad(P0, P1, P2, 0); R0 <- sp1$right[[1]]; R1 <- sp1$right[[2]]; R2 <- sp1$right[[3]]
            sp2 <- split_quad(R0, R1, R2, 1); S0 <- sp2$left[[1]]; S1 <- sp2$left[[2]]; S2 <- sp2$left[[3]]
            bezier_ctrl[[ii]] <- data.frame(x = c(S0[1], S1[1], S2[1]), y = c(S0[2], S1[2], S2[2]), .curve_id = ii, .pt = 1:3)
          }
          bezier_df <- dplyr::bind_rows(bezier_ctrl)
          if (nrow(bezier_df) == 0) return(p_obj)
          bezier_df <- bezier_df[order(bezier_df$.curve_id, bezier_df$.pt), , drop = FALSE]
          return(p_obj + ggforce::geom_bezier2(data = bezier_df, mapping = ggplot2::aes(x = x, y = y, group = .curve_id),
            linetype = link_lt, linewidth = link_lw, alpha = link_al, color = col_link, lineend = "round"))
        }

        df_draw$.curv_abs <- curv_abs
        df_pos <- df_draw[curv_sign > 0, , drop = FALSE]
        df_neg <- df_draw[curv_sign <= 0, , drop = FALSE]
        dist_pos <- dist_link[curv_sign > 0]
        dist_neg <- dist_link[curv_sign <= 0]
        add_curve_side <- function(p_side, df_side, dist_side, sign_side) {
          if (nrow(df_side) == 0) return(p_side)
          curv_med <- sign_side * stats::median(df_side$.curv_abs, na.rm = TRUE)
          p_side + ggplot2::geom_curve(data = df_side, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            curvature = curv_med, linetype = link_lt, linewidth = link_lw, alpha = link_al, color = col_link)
        }
        p_new <- add_curve_side(p_obj, df_pos, dist_pos, 1)
        p_new <- add_curve_side(p_new, df_neg, dist_neg, -1)
        p_new
      } else {
        p_obj + ggplot2::geom_segment(data = df_link, mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          linetype = link_lt, linewidth = link_lw, alpha = link_al, color = col_link)
      }
    }
    add_link_layer_impl
  })()

  # Node-level links (node, nodeinmodule): draw BEFORE nodes and module visualization
  if (isTRUE(comparisons)) {
    if ((link_level == "node" || link_level == "module&node" || link_level == "module&node2") && ncol(compare_matrix) > 0) {
      tmpi <- 1
      for (i in seq_len(ncol(compare_matrix))) {
        gA <- compare_matrix[1, i]
        gB <- compare_matrix[2, i]

        node_links <- graph_info[[gA]]$ggplot_node_df %>%
          dplyr::select(name, x, y) %>%
          dplyr::inner_join(
            graph_info[[gB]]$ggplot_node_df %>%
              dplyr::select(name, x, y),
            by = "name",
            suffix = c("_A", "_B")
          ) %>%
          dplyr::transmute(
            x = x_A,
            y = y_A,
            xend = x_B,
            yend = y_B,
            group_a = gA,
            group_b = gB
          )

        node_links_info <- graph_info[[gA]]$ggplot_node_df %>%
          dplyr::select(name, x, y) %>%
          dplyr::inner_join(
            graph_info[[gB]]$ggplot_node_df %>%
              dplyr::select(name, x, y),
            by = "name",
            suffix = c("_A", "_B")
          ) %>%
          dplyr::transmute(
            link_level = "node",
            group_a = gA,
            group_b = gB,
            source = as.character(name),
            target = as.character(name),
            x = x_A,
            y = y_A,
            xend = x_B,
            yend = y_B
          )

        link_info <- dplyr::bind_rows(link_info, node_links_info)

        pair_key <- paste(sort(c(gA, gB)), collapse = "|")
        col_i <- resolve_link_color(link_color_node, color_v, pair_key, tmpi)
        lt_i <- resolve_link_style(link_linetype_node, pair_key, tmpi)
        al_i <- resolve_link_style(link_alpha_node, pair_key, tmpi)
        lw_i <- resolve_link_style(link_width_node, pair_key, tmpi)
        cross_sign <- cross_channel_lookup[[pair_key]]
        p <- add_link_layer(p, node_links, col_i, link_type = "node", cross_channel_sign = cross_sign, link_lt = lt_i, link_al = al_i, link_lw = lw_i)
        tmpi <- tmpi + 1
      }
    }

    if (link_level == "nodeinmodule" && nrow(Module_information_plot) > 0) {
      tmpi <- 1
      for (ri in seq_len(nrow(Module_information_plot))) {
        row_i <- Module_information_plot[ri, , drop = FALSE]
        gA <- row_i$GroupA
        gB <- row_i$GroupB
        modA <- as.character(row_i$modA)
        modB <- as.character(row_i$modB)

        nodes_A <- graph_info[[gA]]$ggplot_node_df %>%
          dplyr::filter(as.character(Modularity) == modA) %>%
          dplyr::select(name, x, y)
        nodes_B <- graph_info[[gB]]$ggplot_node_df %>%
          dplyr::filter(as.character(Modularity) == modB) %>%
          dplyr::select(name, x, y)

        nodeinmodule_links <- nodes_A %>%
          dplyr::inner_join(nodes_B, by = "name", suffix = c("_A", "_B")) %>%
          dplyr::transmute(
            x = x_A,
            y = y_A,
            xend = x_B,
            yend = y_B,
            group_a = gA,
            group_b = gB
          )

        nodeinmodule_links_info <- nodes_A %>%
          dplyr::inner_join(nodes_B, by = "name", suffix = c("_A", "_B")) %>%
          dplyr::transmute(
            link_level = "nodeinmodule",
            group_a = gA,
            group_b = gB,
            source = as.character(name),
            target = as.character(name),
            x = x_A,
            y = y_A,
            xend = x_B,
            yend = y_B
          )

        if (nrow(nodeinmodule_links) > 0) {
          link_info <- dplyr::bind_rows(link_info, nodeinmodule_links_info)
          pair_key <- paste(sort(c(gA, gB)), collapse = "|")
          col_i <- resolve_link_color(link_color_node, color_v, pair_key, tmpi)
          lt_i <- resolve_link_style(link_linetype_node, pair_key, tmpi)
          al_i <- resolve_link_style(link_alpha_node, pair_key, tmpi)
          lw_i <- resolve_link_style(link_width_node, pair_key, tmpi)
          cross_sign <- cross_channel_lookup[[pair_key]]
          p <- add_link_layer(p, nodeinmodule_links, col_i, link_type = "node", cross_channel_sign = cross_sign, link_lt = lt_i, link_al = al_i, link_lw = lw_i)
          tmpi <- tmpi + 1
        }
      }
    }
  }


  for (index in seq_along(names(graph_info))) {
    edge_df <- graph_info[[index]]$ggplot_edge_df


    if (isTRUE(network_outline) && nrow(graph_info[[index]]$ggplot_node_df) > 0) {
      gname <- names(graph_info)[index]
      color_vec <- if (is.null(network_outline_color) || length(network_outline_color) == 0L) "grey50" else network_outline_color

      color_grp <- if (length(color_vec) == 1L) {
        color_vec
      } else if (!is.null(names(color_vec)) && gname %in% names(color_vec)) {
        color_vec[gname]
      } else {
        color_vec[((index - 1L) %% length(color_vec)) + 1L]
      }

      fill_grp <- if (is.null(network_outline_fill) || length(network_outline_fill) == 0L) {
        NA
      } else if (length(network_outline_fill) == 1L) {
        network_outline_fill
      } else if (!is.null(names(network_outline_fill)) && gname %in% names(network_outline_fill)) {
        network_outline_fill[gname]
      } else {
        network_outline_fill[((index - 1L) %% length(network_outline_fill)) + 1L]
      }
      group_circle_df <- graph_info[[index]]$ggplot_node_df %>%
        dplyr::mutate(.group_outer = 1L)
      circle_n_grp <- max(40, min(300, as.integer(round(8 * sqrt(nrow(group_circle_df))))))
      alpha_grp <- if (is.na(fill_grp)) 1 else network_outline_fill_alpha
      p <- p +
        ggforce::geom_mark_circle(
          data = group_circle_df,
          mapping = ggplot2::aes(x = x, y = y, group = .group_outer),
          fill = fill_grp,
          alpha = alpha_grp,
          color = color_grp,
          linetype = network_outline_linetype,
          linewidth = network_outline_width,
          n = circle_n_grp,
          expand = grid::unit(network_outline_expand, "mm")
        )
    }

    # plot link
    edge_color_mapped <- edge_color %in% colnames(edge_df)
    if (!edge_color_mapped && edge_color_is_sign) {
      stop("`edge_color = \"corr_direction\"` requires `corr_direction` in edge data.")
    }
    if (!edge_color_mapped) {
      p <- p + add_inner_link_geom(
        edge_data = edge_df,
        mapping_obj = ggplot2::aes(x = from_x, xend = to_x, y = from_y, yend = to_y),
        color_fixed = edge_color
      )
    } else {
      line_values <- edge_df[[edge_color]]
      line_scale <- if (is.numeric(line_values)) {
        ggplot2::scale_color_gradient(low = "#4393c3", high = "#d6604d", name = edge_color)
      } else if (!is.null(edge_color_values)) {
        ggplot2::scale_color_manual(values = edge_color_values, name = edge_color)
      } else if (edge_color_is_sign) {
        ggplot2::scale_color_manual(values = c("Positive" = "#d6604d", "Negative" = "#4393c3"),
                                    name = edge_color)
      } else {
        scale_color_ggnetview(unique(as.character(line_values)), name = edge_color)
      }

      p <- p +
        ggnewscale::new_scale_color() +
        add_inner_link_geom(
          edge_data = edge_df,
          mapping_obj = ggplot2::aes(x = from_x, xend = to_x, y = from_y, yend = to_y, colour = .data[[edge_color]])
        ) +
        line_scale
    }

    module_targets <- character(0)
    if (link_level %in% c("module", "nodeinmodule", "module&node") && nrow(Module_information_plot) > 0) {
      module_targets <- Module_information_plot %>%
        # match by exact GroupA/GroupB equality only; the previous
        # str_detect(Group, ...) treated the group name as a regex and was
        # redundant with this line -- group names containing regex metacharacters
        # (or being substrings of each other) would error or mis-match.
        dplyr::filter(GroupA == names(graph_list)[index] | GroupB == names(graph_list)[index]) %>%
        dplyr::mutate(mod_target = dplyr::case_when(
          GroupA == names(graph_list)[index] ~ modA,
          GroupB == names(graph_list)[index] ~ modB,
          .default = NA
        )) %>%
        dplyr::pull(mod_target) %>%
        unique()

      module_targets <- module_targets[module_targets != "Others"]
    } else {
      # For link_level = "none"/"node" (or when no module match exists),
      # outer boundaries should still be visible for the current group's modules.
      module_targets <- graph_info[[index]]$ggplot_node_df %>%
        dplyr::mutate(Modularity = as.character(Modularity)) %>%
        dplyr::pull(Modularity) %>%
        unique()
      module_targets <- module_targets[module_targets != "Others"]
    }

    outer_node_df <- graph_info[[index]]$ggplot_node_df %>%
      dplyr::filter(name %in% (graph_list[[index]] %>%
                                 tidygraph::activate(nodes) %>%
                                 tidygraph::as_tibble() %>%
                                 dplyr::mutate(Modularity = as.character(Modularity)) %>%
                                 dplyr::filter(Modularity %in% module_targets) %>%
                                 dplyr::pull(name)))

    module_levels_group <- graph_info[[index]]$ggplot_node_df %>%
      dplyr::mutate(Modularity = as.character(Modularity)) %>%
      dplyr::pull(Modularity) %>%
      unique()

    fill_scale_group <- if (is.null(node_fill_values)) {
      ggplot2::scale_fill_viridis_d(drop = FALSE, limits = module_levels_group)
    } else {
      ggplot2::scale_fill_manual(values = node_fill_values, drop = FALSE)
    }

    color_is_fixed <- !is.null(node_color_values) && length(node_color_values) == 1 && (is.null(names(node_color_values)) || names(node_color_values)[1] == "")
    color_scale_group <- if (is.null(node_color_values)) {
      ggplot2::scale_color_viridis_d(drop = FALSE, limits = module_levels_group)
    } else if (isTRUE(color_is_fixed)) {
      NULL
    } else {
      ggplot2::scale_color_manual(values = node_color_values, drop = FALSE)
    }

    if (is.null(node_color_values)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(data = graph_info[[index]]$ggplot_node_df,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   fill = Modularity,
                                                   size = Degree),
                            shape = 21) +
        fill_scale_group +
        ggnewscale::new_scale_fill()
    } else if (isTRUE(color_is_fixed)) {
      p <- p +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(data = graph_info[[index]]$ggplot_node_df,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   fill = Modularity,
                                                   size = Degree),
                            shape = 21,
                            color = node_color_values) +
        fill_scale_group +
        ggnewscale::new_scale_fill()
    } else {
      p <- p +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_point(data = graph_info[[index]]$ggplot_node_df,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   fill = Modularity,
                                                   color = Modularity,
                                                   size = Degree),
                            shape = 21) +
        fill_scale_group +
        color_scale_group +
        ggnewscale::new_scale_fill()
    }

    p <- p + ggplot2::scale_size(range = node_size_range, guide = "none")

    if (module_outline == "circle" && nrow(outer_node_df) > 0) {
      circle_n <- max(40, min(300, as.integer(round(8 * sqrt(nrow(outer_node_df))))))
      p <- p +
        ggnewscale::new_scale_color() +
        {
          if (isTRUE(color_is_fixed)) {
            ggforce::geom_mark_circle(
              data = outer_node_df,
              mapping = ggplot2::aes(x = x, y = y, fill = Modularity),
              color = node_color_values,
              n = circle_n,
              expand = grid::unit(1, "mm")
            )
          } else {
            ggforce::geom_mark_circle(
              data = outer_node_df,
              mapping = ggplot2::aes(x = x, y = y, fill = Modularity, color = Modularity),
              n = circle_n,
              expand = grid::unit(1, "mm")
            )
          }
        } +
        fill_scale_group +
        {
          if (isTRUE(color_is_fixed)) ggplot2::guides(color = "none") else color_scale_group
        }
    }

    if (module_outline == "manual" && nrow(outer_node_df) > 2) {
      maskTable <- generateMask_ggnetview(
        dims = outer_node_df %>% dplyr::select(x, y),
        clusters = outer_node_df %>% dplyr::pull(Modularity),
        q = module_outline_q,
        expand = module_outline_expand,
        bandwidth_scale = module_outline_bandwidth
      ) %>%
        dplyr::mutate(cluster = as.character(cluster))
      p <- p +
        ggnewscale::new_scale_color() +
        {
          if (isTRUE(color_is_fixed)) {
            ggplot2::geom_polygon(
              data = maskTable,
              mapping = ggplot2::aes(x = x, y = y,
                                     group = interaction(cluster, polygon_id),
                                     fill = cluster),
              color = node_color_values,
              linewidth = module_outline_width,
              linetype = module_outline_linetype,
              alpha = module_outline_alpha,
              show.legend = FALSE
            )
          } else {
            ggplot2::geom_polygon(
              data = maskTable,
              mapping = ggplot2::aes(x = x, y = y,
                                     group = interaction(cluster, polygon_id),
                                     fill = cluster, color = cluster),
              linewidth = module_outline_width,
              linetype = module_outline_linetype,
              alpha = module_outline_alpha,
              show.legend = FALSE
            )
          }
        } +
        fill_scale_group +
        {
          if (isTRUE(color_is_fixed)) ggplot2::guides(color = "none") else color_scale_group
        }
    }

    p <- p +
      ggplot2::annotate(
        geom = "text",
        x = mean(graph_info[[index]]$ggplot_node_df$x),
        y = max(graph_info[[index]]$ggplot_node_df$y) + group_label_offset,
        label = if (isTRUE(edge_color_is_sign)) {
          paste0("Group = ", names(graph_info[index]), "\n",
                 "Node = ", graph_stat[[index]]$node, "\n",
                 "Edge = ", graph_stat[[index]]$edge, "\n",
                 "Positive = ", graph_stat[[index]]$position_edge, "\n",
                 "Negative = ", graph_stat[[index]]$negative_edge, "\n")
        } else {
          paste0("Group = ", names(graph_info[index]), "\n",
                 "Node = ", graph_stat[[index]]$node, "\n",
                 "Edge = ", graph_stat[[index]]$edge, "\n")
        },
        size = group_label_size,
        fontface = "bold"
      )

  }

  # When link_curve_mode == "cross" (for module links, cross_channel_lookup already set above): assign one curvature sign per group-pair so
  # same-pair links form a bundle along a channel (left/right/bottom arc), not through center.
  # Cross uses the OPPOSITE of the inward rule per pair so that cross = outer arc vs inward = inner arc (they then look distinct).
  cross_channel_lookup <- list()
  circle_centers_x <- circle_centers_y <- circle_radius <- NULL
  circle_pad <- 0
  module_circle_tbl <- NULL
  module_circle_pad <- 0
  if (identical(link_curve_mode, "cross")) {
    gnames <- names(graph_info)
    ngrp <- length(gnames)
    cx <- vapply(gnames, function(g) mean(graph_info[[g]]$ggplot_node_df$x), numeric(1))
    cy <- vapply(gnames, function(g) mean(graph_info[[g]]$ggplot_node_df$y), numeric(1))
    names(cx) <- gnames
    names(cy) <- gnames
    circle_centers_x <- cx
    circle_centers_y <- cy
    # Forbidden region should cover the visible group orbit, not only the center hole.
    circle_radius <- vapply(
      gnames,
      function(g) {
        node_df <- graph_info[[g]]$ggplot_node_df
        d <- sqrt((node_df$x - cx[g])^2 + (node_df$y - cy[g])^2)
        module_outline_q <- as.numeric(stats::quantile(d, probs = 0.98, na.rm = TRUE))
        r_max <- max(d, na.rm = TRUE)
        max(max(module_outline_q, r_max * 0.98), 1e-6)
      },
      numeric(1)
    )
    circle_pad <- max(stats::median(circle_radius, na.rm = TRUE) * 0.06, 1e-6)
    # Also forbid entering module circles (small circles in each group).
    module_rows <- list()
    rid <- 1L
    for (g in gnames) {
      node_df <- graph_info[[g]]$ggplot_node_df
      if (!"Modularity" %in% colnames(node_df)) next
      mods <- unique(as.character(node_df$Modularity))
      mods <- mods[!is.na(mods)]
      for (m in mods) {
        md <- node_df[as.character(node_df$Modularity) == m, , drop = FALSE]
        if (nrow(md) == 0) next
        mcx <- mean(md$x)
        mcy <- mean(md$y)
        mr <- max(sqrt((md$x - mcx)^2 + (md$y - mcy)^2))
        if (!is.finite(mr) || mr <= 0) next
        module_rows[[rid]] <- data.frame(
          group = g,
          module = m,
          key = paste0(g, "|", m),
          cx = mcx,
          cy = mcy,
          r = mr,
          stringsAsFactors = FALSE
        )
        rid <- rid + 1L
      }
    }
    if (length(module_rows) > 0) {
      module_circle_tbl <- dplyr::bind_rows(module_rows)
      module_circle_pad <- max(stats::median(module_circle_tbl$r, na.rm = TRUE) * 0.10, 1e-6)
    }
    top_g <- left_g <- right_g <- NULL
    if (ngrp >= 2L) {
      if (ngrp == 3L) {
        top_ix <- which.max(cy)
        bottom_ix <- setdiff(seq_len(3L), top_ix)
        left_ix <- bottom_ix[which.min(cx[bottom_ix])]
        right_ix <- setdiff(bottom_ix, left_ix)
        top_g <- gnames[top_ix]
        left_g <- gnames[left_ix]
        right_g <- gnames[right_ix]
        # Use opposite of inward so cross = outer arc (bend away from center), inward = inner arc.
        outward_tl <- (cx[top_g] + cx[left_g])/2 * (cy[top_g] - cy[left_g]) + (cy[top_g] + cy[left_g])/2 * (cx[left_g] - cx[top_g])
        cross_channel_lookup[[paste(sort(c(top_g, left_g)), collapse = "|")]] <- as.integer(if (outward_tl != 0) -sign(outward_tl) else -1L)
        outward_tr <- (cx[top_g] + cx[right_g])/2 * (cy[top_g] - cy[right_g]) + (cy[top_g] + cy[right_g])/2 * (cx[right_g] - cx[top_g])
        cross_channel_lookup[[paste(sort(c(top_g, right_g)), collapse = "|")]] <- as.integer(if (outward_tr != 0) -sign(outward_tr) else 1L)
        outward_lr <- (cx[left_g] + cx[right_g])/2 * (cy[left_g] - cy[right_g]) + (cy[left_g] + cy[right_g])/2 * (cx[right_g] - cx[left_g])
        cross_channel_lookup[[paste(sort(c(left_g, right_g)), collapse = "|")]] <- as.integer(if (outward_lr != 0) -sign(outward_lr) else 1L)
      } else {
        cx_cent <- mean(cx)
        cy_cent <- mean(cy)
        for (i in seq_len(ngrp - 1L)) {
          for (j in (i + 1L):ngrp) {
            gA <- gnames[i]
            gB <- gnames[j]
            dx <- cx[gB] - cx[gA]
            dy <- cy[gB] - cy[gA]
            mx <- (cx[gA] + cx[gB]) / 2
            my <- (cy[gA] + cy[gB]) / 2
            s <- as.integer(sign((cx_cent - mx) * dy - (cy_cent - my) * dx))
            key_ij <- paste(sort(c(gA, gB)), collapse = "|")
            cross_channel_lookup[[key_ij]] <- if (is.finite(s) && s != 0L) s else 1L
          }
        }
      }
    }
  }

  if (isTRUE(comparisons) && (link_level == "module" || link_level == "module&node" || link_level == "module&node2" || link_level == "nodeinmodule") && nrow(Module_information_plot) > 0) {
    tmpi <- 1
    for (link_type in names(table(Module_information_plot$Group))) {
      modA_tmp <- Module_information_plot %>%
        dplyr::filter(Group == link_type) %>%
        dplyr::select(modA, GroupA) %>%
        dplyr::pull(modA)

      GroupA_tmp <- Module_information_plot %>%
        dplyr::filter(Group == link_type) %>%
        dplyr::select(modA, GroupA) %>%
        dplyr::pull(GroupA) %>%
        unique()

      modB_tmp <- Module_information_plot %>%
        dplyr::filter(Group == link_type) %>%
        dplyr::select(modB, GroupB) %>%
        dplyr::pull(modB)

      GroupB_tmp <- Module_information_plot %>%
        dplyr::filter(Group == link_type) %>%
        dplyr::select(modB, GroupB) %>%
        dplyr::pull(GroupB) %>%
        unique()

      Module_location <- Module_information_plot %>%
        dplyr::filter(Group == link_type) %>%
        dplyr::left_join(
          graph_info[[GroupA_tmp]]$ggplot_node_df %>%
            dplyr::mutate(Modularity = as.character(Modularity)) %>%
            dplyr::filter(Modularity %in% modA_tmp) %>%
            dplyr::select(x, y, Modularity) %>%
            dplyr::group_by(Modularity) %>%
            dplyr::summarise(x_center = mean(x), y_center = mean(y)) %>%
            purrr::set_names(c("GroupA_Module", "GroupA_x_center", "GroupA_y_center")),
          by = c("modA" = "GroupA_Module")
        ) %>%
        dplyr::left_join(
          graph_info[[GroupB_tmp]]$ggplot_node_df %>%
            dplyr::mutate(Modularity = as.character(Modularity)) %>%
            dplyr::filter(Modularity %in% modB_tmp) %>%
            dplyr::select(x, y, Modularity) %>%
            dplyr::group_by(Modularity) %>%
            dplyr::summarise(x_center = mean(x), y_center = mean(y)) %>%
            purrr::set_names(c("GroupB_Module", "GroupB_x_center", "GroupB_y_center")),
          by = c("modB" = "GroupB_Module")
        )

      Module_location_plot <- Module_location %>%
        dplyr::transmute(
          x = GroupA_x_center,
          y = GroupA_y_center,
          xend = GroupB_x_center,
          yend = GroupB_y_center,
          group_a = GroupA,
          group_b = GroupB
        )

      link_info <- dplyr::bind_rows(
        link_info,
        Module_location %>%
          dplyr::transmute(
            link_level = "module",
            group_a = GroupA,
            group_b = GroupB,
            source = as.character(modA),
            target = as.character(modB),
            x = GroupA_x_center,
            y = GroupA_y_center,
            xend = GroupB_x_center,
            yend = GroupB_y_center
          )
      )

      pair_key <- paste(sort(c(GroupA_tmp, GroupB_tmp)), collapse = "|")
      col_i <- resolve_link_color(link_color_module, color_v, pair_key, tmpi)
      lt_i <- resolve_link_style(link_linetype_module, pair_key, tmpi)
      al_i <- resolve_link_style(link_alpha_module, pair_key, tmpi)
      lw_i <- resolve_link_style(link_width_module, pair_key, tmpi)
      cross_sign <- cross_channel_lookup[[pair_key]]
      p <- add_link_layer(p, Module_location_plot, col_i, link_type = "module", cross_channel_sign = cross_sign, link_lt = lt_i, link_al = al_i, link_lw = lw_i)
      tmpi <- tmpi + 1
    }
  }

  p <- p +
    coord_fixed(clip = FALSE) +
    theme_ggnetview() +
    theme(legend.position = "none")


  return(
    list(p = p,
         info = as_tibble(Module_information),
         link_info = as_tibble(link_info),
         graph = graph_list,
         topology = if (isTRUE(calculate_topology)) {
           list(network = topology_network,
                sample = topology_sample)
         } else {
           NULL
         }
         )
    )


}
