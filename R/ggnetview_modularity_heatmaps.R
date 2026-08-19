#' Compute module eigengene (PC1) from OTU abundance matrix
#'
#' For each module, performs PCA on the module members' abundance matrix
#' (samples x OTUs) and returns the first principal component as the
#' module eigengene. This represents the "overall activity/state" of
#' each module across samples.
#'
#' @param otu_mat Numeric matrix. Rows = OTUs/ASVs, columns = samples.
#'   Must have rownames (OTU IDs) and colnames (sample IDs).
#' @param graph_obj A tbl_graph from build_graph_from_mat or build_graph_from_df.
#'   Must have node attribute \code{name} and a module column.
#' @param module_col Character. Name of the module column in graph_obj nodes.
#'   One of \code{"Modularity"}, \code{"modularity3"}, \code{"modularity2"}.
#' @param exclude_others Logical. If TRUE, exclude "Others" module from output.
#' @param scale_pca Logical. If TRUE, scale variables before PCA (recommended).
#'
#' @return A data frame with samples as rows and module eigengenes as columns.
#'   Column names are module names.
#' @keywords internal
get_module_eigengene <- function(otu_mat,
                                 graph_obj,
                                 module_col = "Modularity",
                                 exclude_others = TRUE,
                                 scale_pca = TRUE) {
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  if (!"name" %in% colnames(node_df)) {
    stop("`graph_obj` nodes must have a `name` column.", call. = FALSE)
  }

  mod_candidates <- c("Modularity", "modularity3", "modularity2")
  if (!module_col %in% colnames(node_df)) {
    hit <- mod_candidates[mod_candidates %in% colnames(node_df)]
    module_col <- if (length(hit) > 0) hit[1] else stop("No module column found.", call. = FALSE)
  }

  if (is.null(rownames(otu_mat)) || is.null(colnames(otu_mat))) {
    stop("`otu_mat` must have rownames (OTU IDs) and colnames (sample IDs).", call. = FALSE)
  }

  otu_mat <- as.matrix(otu_mat)
  node_names <- as.character(node_df$name)
  modules <- as.character(node_df[[module_col]])

  if (exclude_others) {
    keep <- modules != "Others"
    node_names <- node_names[keep]
    modules <- modules[keep]
  }

  mod_levels <- unique(modules)
  mod_levels <- mod_levels[mod_levels != "Others"]
  if (length(mod_levels) == 0) {
    stop("No non-Others modules found in graph_obj.", call. = FALSE)
  }

  eigengene_list <- list()
  for (mod in mod_levels) {
    otu_in_mod <- node_names[modules == mod]
    otu_in_mat <- intersect(otu_in_mod, rownames(otu_mat))
    if (length(otu_in_mat) < 2L) next

    mat_mod <- t(otu_mat[otu_in_mat, , drop = FALSE])
    if (any(!is.finite(mat_mod))) {
      mat_mod[!is.finite(mat_mod)] <- 0
    }
    pca_out <- stats::prcomp(mat_mod, center = TRUE, scale. = scale_pca)
    pc1 <- pca_out$x[, 1L, drop = TRUE]
    eigengene_list[[mod]] <- pc1
  }

  if (length(eigengene_list) == 0) {
    stop("Could not compute eigengene for any module (need >= 2 OTUs per module).", call. = FALSE)
  }

  out <- as.data.frame(eigengene_list)
  colnames(out) <- names(eigengene_list)
  rownames(out) <- colnames(otu_mat)
  out
}


#' Compute module relative abundance index
#'
#' For each module, computes the sum or mean of member OTU abundances
#' per sample. Simpler than eigengene but more "abundance-focused".
#'
#' @param otu_mat Numeric matrix. Rows = OTUs, columns = samples.
#' @param graph_obj A tbl_graph with \code{name} and module column.
#' @param module_col Character. Module column name.
#' @param type Character. \code{"sum"} (total abundance) or \code{"mean"}
#'   (average abundance per OTU).
#' @param exclude_others Logical. Exclude "Others" module.
#'
#' @return A data frame with samples as rows and modules as columns.
#' @keywords internal
get_module_abundance <- function(otu_mat,
                                graph_obj,
                                module_col = "Modularity",
                                type = c("sum", "mean"),
                                exclude_others = TRUE) {
  type <- match.arg(type)
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_candidates <- c("Modularity", "modularity3", "modularity2")
  if (!module_col %in% colnames(node_df)) {
    hit <- mod_candidates[mod_candidates %in% colnames(node_df)]
    module_col <- if (length(hit) > 0) hit[1] else stop("No module column found.", call. = FALSE)
  }

  otu_mat <- as.matrix(otu_mat)
  node_names <- as.character(node_df$name)
  modules <- as.character(node_df[[module_col]])

  if (exclude_others) {
    keep <- modules != "Others"
    node_names <- node_names[keep]
    modules <- modules[keep]
  }

  mod_levels <- unique(modules)
  mod_levels <- mod_levels[mod_levels != "Others"]

  ab_list <- list()
  for (mod in mod_levels) {
    otu_in_mod <- node_names[modules == mod]
    otu_in_mat <- intersect(otu_in_mod, rownames(otu_mat))
    if (length(otu_in_mat) == 0L) next

    mat_mod <- otu_mat[otu_in_mat, , drop = FALSE]
    if (type == "sum") {
      ab_list[[mod]] <- colSums(mat_mod)
    } else {
      ab_list[[mod]] <- colMeans(mat_mod)
    }
  }

  out <- as.data.frame(ab_list)
  colnames(out) <- names(ab_list)
  rownames(out) <- colnames(otu_mat)
  out
}


#' Network + module-environment heatmaps with shared Mantel API
#'
#' @description
#' Render a \code{ggNetView} network in the centre, surrounded by up to four
#' environmental-correlation heatmap quadrants, with link segments
#' connecting each module's anchor to the corresponding env-variable points
#' on the diagonals. Each module is represented by a single per-sample
#' summary (eigengene or abundance) for downstream statistics. Supports
#' both Pearson/Spearman/Kendall correlation and \code{vegan::mantel}
#' tests, and exposes the same Mantel API as
#' \code{\link{gglink_heatmaps}}.
#'
#' @details
#' \strong{Pipeline.}
#' \enumerate{
#'   \item Read module membership from \code{graph_obj} (the node attribute
#'     selected by the package's module column, e.g. \code{"Modularity"}).
#'   \item Build a per-sample module summary matrix \code{spec_df} of shape
#'     (samples x modules) using either \code{module_eigengene} (PC1 of the
#'     module's OTU sub-matrix; recommended) or \code{module_abundance}
#'     (sum / mean of OTU abundances inside the module).
#'   \item For each \code{(env_block, modules)} pair, compute either a
#'     correlation or a Mantel test (see \strong{Mantel API} below).
#'   \item Render the network via \code{\link{ggNetView}} (parameters
#'     forwarded through \code{...}), draw four heatmap quadrants for env-env
#'     correlations, and overlay link segments from module anchors to env
#'     diagonals.
#' }
#'
#' \strong{Mantel API (shared with \code{\link{gglink_heatmaps}}).} Two
#' algorithms exposed via \code{mantel_kind}; both go through the helpers
#' in \code{\link{mantel_block_vs_col}} / \code{\link{mantel_pairwise}} so
#' the two top-level functions stay numerically identical:
#' \itemize{
#'   \item \code{"block_vs_col"} (default, ecological standard): for each
#'     module, the OTUs that belong to it are pulled out of \code{otu_mat}
#'     (transposed to samples x OTUs) and turned into ONE community
#'     distance matrix with \code{spec_dist_method}; each env column is
#'     turned into its own distance matrix with \code{env_dist_method};
#'     one Mantel test per (module, env_col).
#'   \item \code{"col_vs_col"} (legacy): the module's representative vector
#'     (eigengene or abundance) is treated as a single variable; its
#'     single-column distance matrix is tested against each env column's
#'     single-column distance matrix. Mathematically close to a rank
#'     correlation. Kept for backwards compatibility / sensitivity
#'     comparisons.
#' }
#' Prior versions of this function used the equivalent of
#' \code{"col_vs_col"} implicitly. The default has been switched to
#' \code{"block_vs_col"} (with a one-time \code{message()} on the first
#' Mantel call) to match the standard ecological interpretation; pass
#' \code{mantel_kind = "col_vs_col"} to reproduce the old numbers.
#'
#' \strong{Output schema is stable across modes.} The returned stats data
#' frame always has columns \code{ID, Type, Correlation, Pvalue, p_signif,
#' spec_block, env_block, method}. \code{ID} is the module name in all
#' cases (\code{"M1"}, \code{"M2"}, ...), \code{Type} is the env column
#' name. This makes \code{drop_nonsig} and \code{comparisons_groups} work
#' identically across \code{relation_method} / \code{mantel_kind}
#' combinations.
#'
#' @section Data inputs:
#' The graph carrying module assignments, the env data table that defines
#' the heatmap quadrants, the OTU abundance matrix used to summarise each
#' module, and the named list that partitions env into blocks.
#'
#' @param graph_obj A \code{tbl_graph} (e.g. from \code{build_graph_from_mat}
#'   or \code{build_graph_from_df}). Must have a node \code{name} attribute
#'   and a module column (one of \code{"Modularity"}, \code{"modularity3"},
#'   \code{"modularity2"}; auto-detected).
#' @param env Data frame or matrix of environmental variables. Rows are
#'   samples (rownames matched against \code{otu_mat} columns), columns
#'   are env factors.
#' @param otu_mat Numeric matrix. Rows = OTUs/ASVs (rownames matched
#'   against \code{graph_obj} node names), columns = samples (colnames
#'   matched against \code{env} rownames). Used to compute module
#'   eigengenes / abundances and, in block-vs-col Mantel mode, to assemble
#'   per-module community distance matrices.
#' @param env_select Named list (required). Column indices or names of
#'   \code{env} that form each env block, one block per heatmap quadrant.
#'   \code{length(env_select)} must equal \code{length(orientation)}. The
#'   list names (\code{names(env_select)}) are used by
#'   \code{comparisons_groups} and in the returned stats. Example:
#'   \code{list(Env01 = 1:5, Env02 = 6:10, Env03 = 11:15, Env04 = 16:20)}.
#'
#' @section Module representation:
#' How each module's per-sample value is computed before being correlated
#' with env: either as the eigengene (PC1 of the OTU sub-matrix) or as a
#' summary (sum / mean) of within-module abundances.
#'
#' @param module_index Character (default \code{"eigengene"}). How each
#'   module is summarised into one per-sample value used downstream.
#'   \code{"eigengene"} = PC1 of the module's OTU sub-matrix
#'   (recommended); \code{"abundance"} = sum or mean of OTU abundances
#'   within the module (controlled by \code{abundance_type}).
#' @param abundance_type Character (default \code{"sum"}). Only used when
#'   \code{module_index = "abundance"}. Either \code{"sum"} or
#'   \code{"mean"}.
#'
#' @section Statistics -- correlation:
#' Parameters that govern the env-env tile correlations and (when
#' \code{relation_method = "correlation"}) the module-env link
#' correlations: the correlation method and missing-value handling.
#'
#' @param relation_method Character (default \code{"correlation"}). One of
#'   \code{"correlation"} or \code{"mantel"}.
#' @param cor.method Character (default \code{"pearson"}). Correlation
#'   method used by \code{psych::corr.test} for env-env (heatmap tiles)
#'   and, when \code{relation_method = "correlation"}, for module-env
#'   links. One of \code{"pearson"}, \code{"kendall"}, \code{"spearman"}.
#' @param cor.use Character (default \code{"everything"}). Missing-value
#'   handling for \code{psych::corr.test}. One of \code{"everything"},
#'   \code{"all"}, \code{"complete"}, \code{"pairwise"}, \code{"na"}.
#'
#' @section Statistics -- Mantel:
#' Parameters used only when \code{relation_method = "mantel"}: the
#' Mantel variant, the dissimilarity / distance metrics, the Mantel
#' correlation method, and the permutation count.
#'
#' @param mantel_kind Character (default \code{"block_vs_col"}). Which
#'   Mantel algorithm to use; see \strong{Details}. The same parameter is
#'   exposed in \code{\link{gglink_heatmaps}}.
#' @param mantel.method2 Character (default \code{"pearson"}). Correlation
#'   coefficient passed to \code{vegan::mantel} as its \code{method}
#'   argument. One of \code{"pearson"}, \code{"kendall"},
#'   \code{"spearman"}.
#' @param spec_dist_method Character (default \code{"bray"}). Dissimilarity
#'   method (\code{vegan::vegdist}) used to convert a module's OTU
#'   sub-matrix into ONE community distance matrix when
#'   \code{mantel_kind = "block_vs_col"}.
#' @param env_dist_method Character (default \code{"euclidean"}). Distance
#'   method (\code{vegan::vegdist}) used to convert each env column into
#'   its own distance matrix when \code{relation_method = "mantel"}.
#' @param permutations Integer (default \code{999L}). Number of
#'   permutations passed to \code{vegan::mantel}.
#' @param mantel.seed Integer (default \code{1115}). Seed forwarded to the
#'   Mantel helpers so the permutation p-values are reproducible across runs.
#'
#' @section What gets analysed / drawn:
#' Filters and selectors that decide what ends up on the plot: dropping
#' non-significant links, and which heatmap quadrants are rendered.
#'
#' @param drop_nonsig Logical (default \code{FALSE}). If \code{TRUE},
#'   non-significant links (p > 0.05) are removed from the plots; the
#'   returned stats data frame is unaffected.
#' @param orientation Character vector (default
#'   \code{c("top_right","bottom_right","top_left","bottom_left")}). Which
#'   heatmap quadrants to draw, in the same order as \code{env_select}.
#'
#' @section Geometry:
#' Spatial parameters that position the env heatmaps relative to the
#' central network and select the network's own layout: heatmap offset,
#' network radius, overall heatmap scale, and the layout/module-ordering
#' choices forwarded to \code{ggNetView()}.
#'
#' @param distance Numeric (default \code{3}). Offset between the central
#'   network's outer boundary and the env heatmaps. Positive pushes
#'   heatmaps outward; \code{0} places them flush; negative values pull
#'   them inward and may overlap the network.
#' @param r Numeric (default \code{6}). Effective radius for scaling the
#'   central network.
#' @param HeatmapScale Numeric (default \code{1}). Global scale for the
#'   overall heatmap layout. \code{>1} enlarges, \code{<1} shrinks.
#' @param layout Character (default \code{"gephi"}). Layout passed to the
#'   underlying \code{\link{ggNetView}} call (e.g. \code{"gephi"},
#'   \code{"square"}, \code{"WGCNA"}).
#' @param layout_module Character (default \code{"random"}). Module
#'   ordering strategy passed through to \code{\link{ggNetView}}. One of
#'   \code{"random"}, \code{"adjacent"}, \code{"order"}.
#' @param layout.module \lifecycle{deprecated} Use \code{layout_module}.
#'
#' @section Heatmap aesthetics:
#' Visual styling of the env-env heatmap tiles: per-quadrant colour
#' palettes, label and significance-mark sizes, label rotation, the
#' central anchor point on each heatmap, and tile border styling.
#'
#' @param HeatmapColorBar \code{NULL} or list (default \code{NULL}).
#'   Per-quadrant colour palettes. Three accepted forms:
#'   \itemize{
#'     \item \code{NULL}: built-in defaults.
#'     \item Length-2 named list \code{list(low = ..., high = ...)}:
#'       applied to all quadrants.
#'     \item List of length \code{length(orientation)}: each element is
#'       either \code{c(low, high)} or \code{list(low = ..., high = ...)}.
#'       Example:
#'       \code{list(c("#2166ac","#b2182b"), c("#1b7837","#762a83"),
#'                  c("#4393c3","#d6604d"), c("#92c5de","#f4a582"))}.
#'   }
#' @param HeatmapLabelSize Numeric (default \code{5}). Text size for the
#'   heatmap row/column labels.
#' @param HeatmapSigSize Numeric (default \code{5}). Text size for the
#'   significance marks (\code{*}, \code{**}, \code{***}) inside heatmap
#'   tiles.
#' @param HeatmapLabelOrient Numeric (default \code{0}). Rotation angle (in
#'   degrees) for heatmap row/column labels. Try 45 or 90 to avoid label
#'   overlap.
#' @param HeatmapPointSize Numeric (default \code{5}). Point size for the
#'   central module anchor where the heatmap link lands.
#' @param HeatmapPointFill Character (default \code{"#de77ae"}). Fill
#'   colour for the central module anchor point.
#' @param HeatmapTileColor Character or \code{NA} (default \code{NA}).
#'   Border colour for heatmap tiles.
#' @param HeatmapTileSize Numeric (default \code{0}). Border line width
#'   for heatmap tiles.
#'
#' @section Link line aesthetics:
#' Visual styling of the module-env link segments: line-width range
#' (mapped from p-value), colour gradient (mapped from correlation /
#' Mantel r), and overall transparency.
#'
#' @param SigLineWidth Numeric vector of length 2 (default
#'   \code{c(0.5, 2)}). Min / max line width for module-env links; mapped
#'   from \code{-log10(p-value)} so smaller p -> thicker line.
#' @param SigLineColor Character vector of length 2 (default
#'   \code{c("#fdbb84", "#d7301f")}). Colour gradient for module-env
#'   links, mapped from low / high correlation (or Mantel r).
#' @param SigLineAlpha Numeric in \code{[0, 1]} (default \code{0.5}).
#'   Transparency for module-env link segments.
#'
#' @section Forwarded to ggNetView:
#' Extra arguments captured via \code{...} and forwarded verbatim to the
#' underlying \code{\link{ggNetView}} call. Use them to customise the
#' central network's appearance (labels, fills, jitter, outer rings,
#' point sizes) without leaving this wrapper.
#'
#' @param ... Additional arguments forwarded to the underlying
#'   \code{\link{ggNetView}} network call. Commonly used:
#'   \code{shrink}, \code{inner_shrink} (intra-module compactness, only
#'   for \code{layout = "WGCNA"}), \code{node_jitter}, \code{module_outline},
#'   \code{network_outline}, \code{module_label} (logical or character --
#'   module labels in ggNetView style), \code{module_label_size},
#'   \code{module_label_segment_width}, \code{module_label_segment_alpha},
#'   \code{node_fill_values}, \code{node_size_range}. Deprecated pre-0.2.0
#'   names (\code{add_outer}, \code{label}, \code{fill}, ...) are still
#'   accepted with a lifecycle warning.
#'
#' @return A list of length 3:
#' \describe{
#'   \item{\code{[[1]]}}{ggplot object with straight link segments
#'     (\code{geom_segment}).}
#'   \item{\code{[[2]]}}{ggplot object with curved link segments
#'     (\code{geom_curve}).}
#'   \item{\code{[[3]]}}{Data frame of module-env stats (unfiltered, not
#'     affected by \code{drop_nonsig}). Columns: \code{ID} (module name),
#'     \code{Type} (env column name), \code{Correlation}, \code{Pvalue},
#'     \code{p_signif}, \code{spec_block}, \code{env_block}, \code{method}
#'     (\code{"correlation"} or \code{"mantel"}). Schema is identical
#'     across all \code{relation_method} / \code{mantel_kind} combinations.}
#' }
#'
#' @seealso
#' \code{\link{gglink_heatmaps}} for the spec-select counterpart that
#' shares the same Mantel API;
#' \code{\link{mantel_block_vs_col}}, \code{\link{mantel_pairwise}} for
#' the underlying Mantel implementations;
#' \code{\link{ggNetView}} for the central network rendering and the
#' parameters forwarded via \code{...}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Default: correlation, eigengene as module summary, four quadrants.
#' p <- ggnetview_modularity_heatmaps(
#'   graph_obj  = g,
#'   env        = env,
#'   otu_mat    = otu_mat,
#'   env_select = list(Env01 = 1:5, Env02 = 6:10,
#'                     Env03 = 11:15, Env04 = 16:20)
#' )
#' p[[1]]            # straight links
#' head(p[[3]])      # stats data frame
#'
#' # Ecologically standard Mantel: one test per (module, env_col).
#' p2 <- ggnetview_modularity_heatmaps(
#'   graph_obj        = g,
#'   env              = env,
#'   otu_mat          = otu_mat,
#'   env_select       = list(Env01 = 1:5, Env02 = 6:10,
#'                           Env03 = 11:15, Env04 = 16:20),
#'   relation_method  = "mantel",
#'   mantel_kind      = "block_vs_col",
#'   spec_dist_method = "bray",
#'   env_dist_method  = "euclidean",
#'   permutations     = 999
#' )
#'
#' # Reproduce legacy column-vs-column Mantel results.
#' p3 <- ggnetview_modularity_heatmaps(
#'   graph_obj       = g,
#'   env             = env,
#'   otu_mat         = otu_mat,
#'   env_select      = list(Env01 = 1:5, Env02 = 6:10,
#'                          Env03 = 11:15, Env04 = 16:20),
#'   relation_method = "mantel",
#'   mantel_kind     = "col_vs_col"
#' )
#' }
ggnetview_modularity_heatmaps <- function(
    graph_obj,
    env,
    otu_mat,
    env_select = NULL,
    module_index = c("eigengene", "abundance"),
    abundance_type = c("sum", "mean"),
    relation_method = c("correlation", "mantel"),
    cor.method = c("pearson", "kendall", "spearman"),
    cor.use = c("everything", "all", "complete", "pairwise", "na"),
    mantel.method2 = c("pearson", "kendall", "spearman"),
    mantel_kind = c("block_vs_col", "col_vs_col"),
    spec_dist_method = "bray",
    env_dist_method = "euclidean",
    permutations = 999L,
    mantel.seed = 1115,
    drop_nonsig = FALSE,
    layout = "gephi",
    layout_module = c("random", "adjacent", "order"),
    orientation = c("top_right", "bottom_right", "top_left", "bottom_left"),
    distance = 3,
    r = 6,
    HeatmapScale = 1,
    SigLineAlpha = 0.5,
    HeatmapLabelSize = 5,
    HeatmapSigSize = 5,
    HeatmapColorBar = NULL,
    HeatmapLabelOrient = 0,
    SigLineWidth = c(0.5, 2),
    SigLineColor = c("#fdbb84", "#d7301f"),
    HeatmapPointSize = 5,
    HeatmapPointFill = "#de77ae",
    HeatmapTileColor = NA,
    HeatmapTileSize = 0,
    ...,
    layout.module = deprecated()) {

  module_index <- match.arg(module_index)
  abundance_type <- match.arg(abundance_type)
  relation_method <- match.arg(relation_method)
  cor.method <- match.arg(cor.method)
  cor.use <- match.arg(cor.use)
  mantel.method2 <- match.arg(mantel.method2)
  mantel_kind <- match.arg(mantel_kind)
  if (!is.character(spec_dist_method) || length(spec_dist_method) != 1L || !nzchar(spec_dist_method)) {
    stop("`spec_dist_method` must be a single non-empty character string (passed to `vegan::vegdist`).", call. = FALSE)
  }
  if (!is.character(env_dist_method) || length(env_dist_method) != 1L || !nzchar(env_dist_method)) {
    stop("`env_dist_method` must be a single non-empty character string (passed to `vegan::vegdist`).", call. = FALSE)
  }
  permutations <- as.integer(permutations)
  if (length(permutations) != 1L || is.na(permutations) || permutations < 1L) {
    stop("`permutations` must be a single positive integer.", call. = FALSE)
  }
  if (relation_method == "mantel") {
    message("Using `mantel_kind = \"", mantel_kind, "\"`. ",
            "Note: prior versions ran the equivalent of `\"col_vs_col\"` when ",
            "`relation_method = \"mantel\"`. The new default is the ecologically ",
            "standard `\"block_vs_col\"` (community matrix vs each env column). ",
            "Pass `mantel_kind = \"col_vs_col\"` to reproduce old results.")
  }
  orientation <- match.arg(orientation, several.ok = TRUE)

  distance <- as.numeric(distance)
  if (length(distance) != 1 || is.na(distance) || distance < 0) {
    stop("`distance` must be a single non-negative numeric value.", call. = FALSE)
  }
  r <- as.numeric(r)
  if (length(r) != 1 || is.na(r) || r <= 0) {
    stop("`r` must be a single positive numeric value.", call. = FALSE)
  }
  HeatmapScale <- as.numeric(HeatmapScale)
  SigLineAlpha <- as.numeric(SigLineAlpha)
  if (length(HeatmapScale) != 1 || !is.finite(HeatmapScale) || HeatmapScale <= 0) {
    stop("`HeatmapScale` must be a single positive numeric value.", call. = FALSE)
  }
  if (length(SigLineAlpha) != 1 || !is.finite(SigLineAlpha) || SigLineAlpha < 0 || SigLineAlpha > 1) {
    stop("`SigLineAlpha` must be a single numeric value between 0 and 1.", call. = FALSE)
  }

  if (is.null(env_select)) {
    stop("`env_select` must be provided (named list of env column indices/names).", call. = FALSE)
  }

  env <- as.data.frame(env)
  if (is.data.frame(otu_mat)) otu_mat <- as.matrix(otu_mat)

  if (length(env_select) != length(orientation)) {
    stop("`env_select` must have the same length as `orientation` (one env block per heatmap quadrant).", call. = FALSE)
  }
  env_list <- purrr::map(env_select, ~ env[, .x, drop = FALSE])
  env_block_names <- names(env_list)
  if (is.null(env_block_names)) env_block_names <- paste0("Env", seq_along(env_list))
  names(env_list) <- orientation

  k_vec <- purrr::map_int(env_list, ncol)
  length_dist <- max(k_vec)
  k_gap <- length_dist - k_vec
  names(k_gap) <- orientation

  sample_ids <- rownames(env)
  if (is.null(sample_ids)) sample_ids <- as.character(seq_len(nrow(env)))

  otu_samples <- colnames(otu_mat)
  common_samples <- intersect(sample_ids, otu_samples)
  if (length(common_samples) < 3L) {
    stop("`env` and `otu_mat` must share at least 3 sample IDs.", call. = FALSE)
  }

  env <- env[common_samples, , drop = FALSE]
  otu_mat <- otu_mat[, common_samples, drop = FALSE]

  if (module_index == "eigengene") {
    spec_df <- get_module_eigengene(otu_mat, graph_obj, exclude_others = TRUE)
  } else {
    spec_df <- get_module_abundance(otu_mat, graph_obj, type = abundance_type, exclude_others = TRUE)
  }

  spec_df <- spec_df[common_samples, , drop = FALSE]
  spec_list <- list(Modules = spec_df)
  spec_block_names <- "Modules"

  env_list <- purrr::map(env_select, ~ env[common_samples, .x, drop = FALSE])
  names(env_list) <- env_block_names

  pairs_to_compute <- expand.grid(
    env_block = env_block_names,
    spec_block = spec_block_names,
    stringsAsFactors = FALSE
  )
  pairs_to_compute <- as.matrix(pairs_to_compute[, c("env_block", "spec_block")])

  if (relation_method == "correlation") {
    cor_parts <- list()
    for (col in seq_len(nrow(pairs_to_compute))) {
      env_blk <- pairs_to_compute[col, 1]
      spec_blk <- pairs_to_compute[col, 2]
      j <- which(spec_block_names == spec_blk)
      p <- which(env_block_names == env_blk)
      if (length(j) != 1 || length(p) != 1) next
      cor_out <- psych::corr.test(spec_list[[j]], env_list[[p]], use = cor.use, method = cor.method)
      cor_r <- cor_out$r %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID, names_to = "Type", values_to = "Correlation")
      cor_p <- cor_out$p %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "ID") %>%
        tidyr::pivot_longer(cols = -ID, names_to = "Type", values_to = "Pvalue") %>%
        dplyr::mutate(p_signif = dplyr::case_when(
          Pvalue > 0.05 ~ "",
          Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
          Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
          Pvalue < 0.001 ~ "***",
          TRUE ~ ""
        ))
      cor_parts[[length(cor_parts) + 1L]] <- cbind(
        cor_r,
        cor_p %>% dplyr::select(Pvalue, p_signif)
      ) %>% dplyr::mutate(spec_block = spec_blk, env_block = env_blk)
    }
    cor_spec_env_list_out <- if (length(cor_parts) > 0) {
      do.call(rbind, cor_parts) %>% dplyr::mutate(method = "correlation")
    } else {
      tibble::tibble(ID = character(), Type = character(), Correlation = numeric(),
                    Pvalue = numeric(), p_signif = character(),
                    spec_block = character(), env_block = character(), method = "correlation")
    }
  } else {
    # Mantel test. Two algorithms exposed via `mantel_kind`, both shared with
    # gglink_heatmaps() through helpers in mantel_utils.R:
    #   - "block_vs_col": for each module, pull out that module's OTUs from
    #     otu_mat, transpose to (samples x OTUs), and let mantel_block_vs_col()
    #     build ONE community-distance matrix per module that is tested
    #     against each env column. ID = module name (e.g. "M1").
    #   - "col_vs_col": treat each module's representative vector
    #     (eigengene/abundance, already in spec_df) as a single variable and
    #     run column-vs-column mantel via mantel_pairwise(). ID = module name
    #     (column name in spec_df).
    cor_parts <- list()

    if (mantel_kind == "block_vs_col") {
      module_members <- get_module_members(graph_obj, exclude_others = TRUE)
      module_names <- names(module_members)
      # Restrict to modules that appear in spec_df (e.g. those with >= 2 OTUs
      # actually computed for eigengene); preserves alignment with the plot.
      module_names <- intersect(module_names, colnames(spec_df))
    }

    for (col in seq_len(nrow(pairs_to_compute))) {
      env_blk <- pairs_to_compute[col, 1]
      spec_blk <- pairs_to_compute[col, 2]
      j <- which(spec_block_names == spec_blk)
      p <- which(env_block_names == env_blk)
      if (length(j) != 1 || length(p) != 1) next

      env_block_df <- env_list[[p]]

      if (mantel_kind == "block_vs_col") {
        # One Mantel test per (module, env_col) pair
        for (mod in module_names) {
          otu_in_mod <- intersect(module_members[[mod]], rownames(otu_mat))
          if (length(otu_in_mod) < 2L) next
          # samples x OTUs community submatrix
          comm_sub <- t(otu_mat[otu_in_mod, , drop = FALSE])
          mout <- mantel_block_vs_col(
            spec_df          = comm_sub,
            env_df           = env_block_df,
            block_name       = mod,
            method           = mantel.method2,
            spec_dist_method = spec_dist_method,
            env_dist_method  = env_dist_method,
            permutations     = permutations,
            na_omit          = TRUE,
            seed             = mantel.seed
          )
          if (base::nrow(mout) == 0L) next
          mout <- mout %>%
            dplyr::mutate(
              spec_block = spec_blk,
              env_block  = env_blk,
              p_signif = dplyr::case_when(
                is.na(Pvalue) ~ "",
                Pvalue > 0.05 ~ "",
                Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
                Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
                Pvalue < 0.001 ~ "***",
                TRUE ~ ""
              )
            )
          cor_parts[[length(cor_parts) + 1L]] <- mout
        }
      } else {
        # col_vs_col: legacy; ID = module column name in spec_df
        mout <- mantel_pairwise(
          spec_df      = spec_list[[j]],
          env_df       = env_block_df,
          method       = mantel.method2,
          permutations = permutations,
          na_omit      = TRUE,
          seed         = mantel.seed
        )
        if (base::nrow(mout) == 0L) next
        mout <- mout %>%
          dplyr::mutate(
            spec_block = spec_blk,
            env_block  = env_blk,
            p_signif = dplyr::case_when(
              is.na(Pvalue) ~ "",
              Pvalue > 0.05 ~ "",
              Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
              Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
              Pvalue < 0.001 ~ "***",
              TRUE ~ ""
            )
          )
        cor_parts[[length(cor_parts) + 1L]] <- mout
      }
    }
    cor_spec_env_list_out <- if (length(cor_parts) > 0) {
      do.call(rbind, cor_parts) %>% dplyr::mutate(method = "mantel")
    } else {
      tibble::tibble(ID = character(), Type = character(), Correlation = numeric(),
                    Pvalue = numeric(), p_signif = character(),
                    spec_block = character(), env_block = character(), method = "mantel")
    }
  }

  env_cor_self_list <- list()
  for (i in seq_along(orientation)) {
    ori <- orientation[i]
    cor_out_self <- psych::corr.test(env_list[[i]], use = cor.use, method = cor.method)
    cor_self_r <- cor_out_self$r %>% as.data.frame()
    cor_self_p <- cor_out_self$p %>% as.data.frame()

    if (ori %in% c("top_right", "bottom_right")) {
      cor_self_r[upper.tri(cor_self_r)] <- NA
      cor_self_p[upper.tri(cor_self_p)] <- NA
    } else {
      cor_self_r[lower.tri(cor_self_r)] <- NA
      cor_self_p[lower.tri(cor_self_p)] <- NA
    }

    cor_self_r <- cor_self_r %>%
      tibble::rownames_to_column(var = "ID") %>%
      tidyr::pivot_longer(cols = -ID, names_to = "Type", values_to = "Correlation") %>%
      dplyr::mutate(
        ID = factor(ID, levels = unique(ID), ordered = TRUE),
        Type = factor(Type, levels = if (ori %in% c("top_right", "bottom_left")) rev(unique(Type)) else unique(Type), ordered = TRUE),
        ID2 = as.numeric(ID),
        Type2 = as.numeric(Type)
      ) %>%
      stats::na.omit()

    cor_self_p <- cor_self_p %>%
      tibble::rownames_to_column(var = "ID") %>%
      tidyr::pivot_longer(cols = -ID, names_to = "Type", values_to = "Pvalue") %>%
      dplyr::mutate(
        ID = factor(ID, levels = unique(ID), ordered = TRUE),
        Type = factor(Type, levels = if (ori %in% c("top_right", "bottom_left")) rev(unique(Type)) else unique(Type), ordered = TRUE),
        ID2 = as.numeric(ID),
        Type2 = as.numeric(Type)
      ) %>%
      stats::na.omit() %>%
      dplyr::mutate(p_signif = dplyr::case_when(
        Pvalue > 0.05 ~ "",
        Pvalue > 0.01 & Pvalue <= 0.05 ~ "*",
        Pvalue <= 0.01 & Pvalue >= 0.001 ~ "**",
        Pvalue < 0.001 ~ "***",
        TRUE ~ ""
      ))

    env_cor_self_list[[i]] <- cbind(
      cor_self_r %>% dplyr::select(ID, Type, ID2, Type2, Correlation),
      cor_self_p %>% dplyr::select(Pvalue, p_signif)
    )
  }
  names(env_cor_self_list) <- orientation

  .diag_xy <- function(id_idx, type_idx, ori, k_gap, length_dist, side_anchor, heatmap_step) {
    x_anchor <- if (ori %in% c("top_right", "bottom_right")) side_anchor[["right"]] else side_anchor[["left"]]
    y_anchor <- if (ori %in% c("top_right", "top_left")) side_anchor[["top"]] else side_anchor[["bottom"]]

    x_out <- if (ori %in% c("top_right", "bottom_right")) {
      x_anchor + heatmap_step * (k_gap[[ori]] + id_idx - 1)
    } else {
      -x_anchor - heatmap_step * (length_dist - id_idx)
    }

    y_out <- if (ori %in% c("top_right", "top_left")) {
      y_anchor + heatmap_step * (k_gap[[ori]] + type_idx - 2)
    } else {
      -y_anchor - heatmap_step * (length_dist - type_idx - 1)
    }

    list(x = x_out, y = y_out)
  }

  .make_targets <- function(df, ori, k_gap, length_dist, side_anchor, heatmap_step) {
    df_diag <- df %>%
      dplyr::mutate(ID = as.character(ID), Type = as.character(Type)) %>%
      dplyr::filter(ID == Type)

    xy_diag <- .diag_xy(
      id_idx = df_diag$ID2,
      type_idx = df_diag$Type2,
      ori = ori,
      k_gap = k_gap,
      length_dist = length_dist,
      side_anchor = side_anchor,
      heatmap_step = heatmap_step
    )

    df_diag %>%
      dplyr::transmute(ID, x_to = xy_diag$x, y_to = xy_diag$y)
  }

  if (lifecycle::is_present(layout.module)) {
    lifecycle::deprecate_warn("0.2.0", "ggnetview_modularity_heatmaps(layout.module)",
                              "ggnetview_modularity_heatmaps(layout_module)")
    layout_module <- layout.module
  }
  layout_module <- match.arg(layout_module)
  # translate deprecated (< 0.2.0) ggNetView() argument names in `...`
  ggnetview_args <- .ggnv_rename_args(list(...), fn = "ggnetview_modularity_heatmaps",
                                      env = environment(), user_env = parent.frame())

  gv_call_args <- list(
    graph_obj = graph_obj,
    layout = layout,
    layout_module = layout_module,
    r = 1,
    scale_radius = r,
    return_layout = TRUE,
    node_add = if (!is.null(ggnetview_args$node_add)) ggnetview_args$node_add else 7,
    scale = if (!is.null(ggnetview_args$scale)) ggnetview_args$scale else TRUE,
    anchor_dist = if (!is.null(ggnetview_args$anchor_dist)) ggnetview_args$anchor_dist else 6,
    orientation = if (!is.null(ggnetview_args$orientation)) ggnetview_args$orientation else "up",
    angle = if (!is.null(ggnetview_args$angle)) ggnetview_args$angle else 0,
    center = if (!is.null(ggnetview_args$center)) ggnetview_args$center else TRUE,
    shrink = if (!is.null(ggnetview_args$shrink)) ggnetview_args$shrink else 1,
    inner_shrink = if (!is.null(ggnetview_args$inner_shrink)) ggnetview_args$inner_shrink else 1,
    node_jitter = if (!is.null(ggnetview_args$node_jitter)) ggnetview_args$node_jitter else FALSE,
    node_jitter_sd = if (!is.null(ggnetview_args$node_jitter_sd)) ggnetview_args$node_jitter_sd else 0.1,
    k_nn = if (!is.null(ggnetview_args$k_nn)) ggnetview_args$k_nn else 12,
    push_others_delta = if (!is.null(ggnetview_args$push_others_delta)) ggnetview_args$push_others_delta else 0,
    module_outline = if (!is.null(ggnetview_args$module_outline)) ggnetview_args$module_outline else FALSE,
    node_fill_values = ggnetview_args$node_fill_values,
    seed = if (!is.null(ggnetview_args$seed)) ggnetview_args$seed else 1115
  )
  # only forward node_size_range when the caller supplied one (ggNetView
  # validates it and must otherwise keep its own default)
  if (!is.null(ggnetview_args$node_size_range)) {
    gv_call_args$node_size_range <- ggnetview_args$node_size_range
  }
  if (layout %in% c("consensus_module_equal_gephi", "consensus_module_gephi")) {
    gv_call_args$nrow <- ggnetview_args$nrow
    gv_call_args$ncol <- ggnetview_args$ncol
  }
  gv_formals <- names(formals(ggNetView))
  gv_call_args <- gv_call_args[names(gv_call_args) %in% gv_formals]

  gv_out <- do.call(ggNetView, gv_call_args)
  if (!is.list(gv_out) || !"layout_data" %in% names(gv_out)) {
    stop("ggNetView did not return layout_data. Ensure layout produces modular structure.", call. = FALSE)
  }

  layout_data <- gv_out$layout_data
  graph_ly_scaled <- layout_data$graph_ly_final
  ggplot_data <- layout_data$ggplot_data
  module_centroids <- layout_data$module_centroids %>%
    dplyr::filter(ID %in% colnames(spec_df))

  if (nrow(module_centroids) == 0) {
    stop("No module centroids could be computed. No overlap between modules and spec columns.", call. = FALSE)
  }

  net_bounds <- list(
    xmin = min(graph_ly_scaled$x, na.rm = TRUE),
    xmax = max(graph_ly_scaled$x, na.rm = TRUE),
    ymin = min(graph_ly_scaled$y, na.rm = TRUE),
    ymax = max(graph_ly_scaled$y, na.rm = TRUE)
  )

  # Heatmap step and anchor are fixed (independent of r). When r increases,
  # only the network scales; heatmaps stay at the same position and size,
  # since distance has not changed.
  r_ref <- 6
  heatmap_step <- max(1, 0.5 * r_ref / max(length_dist, 1)) * HeatmapScale
  diag_default <- min(max(abs(graph_ly_scaled$x), na.rm = TRUE),
                      max(abs(graph_ly_scaled$y), na.rm = TRUE)) / sqrt(2)
  .quad_reach <- function(df, ori, default_val) {
    val <- switch(
      ori,
      top_right = {
        idx <- df$x >= 0 & df$y >= 0
        if (!any(idx)) NA_real_ else max(pmin(df$x[idx], df$y[idx]), na.rm = TRUE)
      },
      bottom_right = {
        idx <- df$x >= 0 & df$y <= 0
        if (!any(idx)) NA_real_ else max(pmin(df$x[idx], -df$y[idx]), na.rm = TRUE)
      },
      top_left = {
        idx <- df$x <= 0 & df$y >= 0
        if (!any(idx)) NA_real_ else max(pmin(-df$x[idx], df$y[idx]), na.rm = TRUE)
      },
      bottom_left = {
        idx <- df$x <= 0 & df$y <= 0
        if (!any(idx)) NA_real_ else max(pmin(-df$x[idx], -df$y[idx]), na.rm = TRUE)
      }
    )
    if (!is.finite(val)) default_val else val
  }
  # Anchor stays at fixed position (r_ref scale) regardless of r; distance unchanged
  quad_reach_raw <- stats::setNames(
    vapply(orientation, function(ori) .quad_reach(graph_ly_scaled, ori, diag_default), numeric(1)),
    orientation
  )
  quad_anchor <- (quad_reach_raw / r) * r_ref + distance
  fallback_anchor <- (diag_default / r) * r_ref + distance
  side_anchor <- c(
    right = max(c(quad_anchor[intersect(c("top_right", "bottom_right"), orientation)], fallback_anchor), na.rm = TRUE),
    left = max(c(quad_anchor[intersect(c("top_left", "bottom_left"), orientation)], fallback_anchor), na.rm = TRUE),
    top = max(c(quad_anchor[intersect(c("top_right", "top_left"), orientation)], fallback_anchor), na.rm = TRUE),
    bottom = max(c(quad_anchor[intersect(c("bottom_right", "bottom_left"), orientation)], fallback_anchor), na.rm = TRUE)
  )

  xy_targets <- purrr::imap_dfr(
    env_cor_self_list[orientation],
    .make_targets,
    k_gap = k_gap,
    length_dist = length_dist,
    side_anchor = side_anchor,
    heatmap_step = heatmap_step
  )

  cor_spec_env_location <- cor_spec_env_list_out %>%
    dplyr::mutate(ID = as.character(ID), Type = as.character(Type)) %>%
    dplyr::left_join(module_centroids, by = "ID") %>%
    dplyr::left_join(xy_targets, by = c("Type" = "ID")) %>%
    dplyr::mutate(line_type = dplyr::if_else(.data$Pvalue <= 0.05, "solid", "dashed"))

  link_df <- cor_spec_env_location %>%
    dplyr::mutate(
      # -log10(p): smaller p -> larger value -> thicker line
      # p=0.05->1.3, p=0.01->2, p=0.001->3
      sig_strength = -log10(.data$Pvalue)
    )
  if (isTRUE(drop_nonsig)) {
    link_df <- link_df %>% dplyr::filter(.data$Pvalue <= 0.05)
  }

  .offset_env <- function(df, ori, k_gap, length_dist, side_anchor, heatmap_step,
                          HeatmapLabelOrient = 0, y_top_all = NULL,
                          y_bottom_all = NULL) {
    df <- df %>% dplyr::mutate(ID = as.character(ID), Type = as.character(Type))
    x_anchor <- if (ori %in% c("top_right", "bottom_right")) side_anchor[["right"]] else side_anchor[["left"]]
    y_anchor <- if (ori %in% c("top_right", "top_left")) side_anchor[["top"]] else side_anchor[["bottom"]]

    x_tile <- if (ori %in% c("top_right", "bottom_right")) {
      x_anchor + heatmap_step * (k_gap[[ori]] + df$ID2 - 1)
    } else {
      -x_anchor - heatmap_step * (length_dist - df$ID2)
    }
    y_tile <- if (ori %in% c("top_right", "top_left")) {
      y_anchor + heatmap_step * (k_gap[[ori]] + df$Type2 - 1)
    } else {
      -y_anchor - heatmap_step * (length_dist - df$Type2)
    }
    tile <- df %>% dplyr::mutate(x_tile = x_tile, y_tile = y_tile, orientation = ori)
    diag_df <- df %>% dplyr::filter(ID == Type)
    diag_xy <- .diag_xy(
      id_idx = diag_df$ID2,
      type_idx = diag_df$Type2,
      ori = ori,
      k_gap = k_gap,
      length_dist = length_dist,
      side_anchor = side_anchor,
      heatmap_step = heatmap_step
    )
    x_diag <- diag_xy$x
    y_diag <- diag_xy$y
    diag <- diag_df %>% dplyr::transmute(ID, x_diag, y_diag, orientation = ori)
    if (HeatmapLabelOrient == 0 || is.null(y_top_all) || is.null(y_bottom_all)) {
      y_id_lab <- if (ori %in% c("top_right", "top_left")) {
        y_anchor + heatmap_step * length_dist
      } else {
        -y_anchor - heatmap_step * length_dist
      }
    } else {
      y_id_lab <- if (ori %in% c("top_right", "top_left")) {
        y_top_all + heatmap_step
      } else {
        y_bottom_all - heatmap_step
      }
    }
    x_type_lab <- if (ori %in% c("top_right", "bottom_right")) {
      x_anchor + heatmap_step * length_dist
    } else {
      -x_anchor - heatmap_step * length_dist
    }
    hjust_type <- if (ori %in% c("top_right", "bottom_right")) "left" else "right"
    id_lab <- df %>%
      dplyr::distinct(ID, .keep_all = TRUE) %>%
      dplyr::transmute(
        ID,
        x_id = if (ori %in% c("top_right", "bottom_right")) {
          x_anchor + heatmap_step * (k_gap[[ori]] + ID2 - 1)
        } else {
          -x_anchor - heatmap_step * (length_dist - ID2)
        },
        y_id = y_id_lab,
        orientation = ori
      )
    type_lab <- df %>%
      dplyr::distinct(Type, .keep_all = TRUE) %>%
      dplyr::transmute(
        Type,
        x_type = x_type_lab,
        y_type = if (ori %in% c("top_right", "top_left")) {
          y_anchor + heatmap_step * (k_gap[[ori]] + Type2 - 1)
        } else {
          -y_anchor - heatmap_step * (length_dist - Type2)
        },
        hjust_type = hjust_type,
        orientation = ori
      )
    list(tile = tile, diag = diag, id_lab = id_lab, type_lab = type_lab)
  }

  .compute_y_range <- function(df, ori, k_gap, length_dist, side_anchor, heatmap_step) {
    df <- df %>% dplyr::mutate(ID = as.character(ID), Type = as.character(Type))
    y_anchor <- if (ori %in% c("top_right", "top_left")) side_anchor[["top"]] else side_anchor[["bottom"]]
    y_tile <- if (ori %in% c("top_right", "top_left")) {
      y_anchor + heatmap_step * (k_gap[[ori]] + df$Type2 - 1)
    } else {
      -y_anchor - heatmap_step * (length_dist - df$Type2)
    }
    data.frame(orientation = ori, ymin = min(y_tile, na.rm = TRUE), ymax = max(y_tile, na.rm = TRUE))
  }

  y_ranges <- purrr::imap_dfr(
    env_cor_self_list[orientation],
    ~ .compute_y_range(.x, .y, k_gap, length_dist, side_anchor, heatmap_step)
  )
  y_top_all <- max(y_ranges$ymax[y_ranges$orientation %in% c("top_right", "top_left")], na.rm = TRUE)
  y_bottom_all <- min(y_ranges$ymin[y_ranges$orientation %in% c("bottom_right", "bottom_left")], na.rm = TRUE)

  packs <- purrr::imap(
    env_cor_self_list[orientation],
    ~ .offset_env(.x, .y, k_gap, length_dist, side_anchor, heatmap_step,
                  HeatmapLabelOrient = HeatmapLabelOrient,
                  y_top_all = y_top_all,
                  y_bottom_all = y_bottom_all)
  )

  diag_all <- purrr::map_dfr(packs, "diag")

  pal_default_low <- c("#4d9221", "#8073ac", "#4393c3", "#66bd63")
  pal_default_high <- c("#c51b7d", "#e08214", "#d6604d", "#f46d43")
  n_quad <- length(packs)
  if (!is.null(HeatmapColorBar) && is.list(HeatmapColorBar) &&
      length(HeatmapColorBar) == 2 && all(c("low", "high") %in% names(HeatmapColorBar))) {
    pal_low <- rep(HeatmapColorBar$low, length.out = n_quad)
    pal_high <- rep(HeatmapColorBar$high, length.out = n_quad)
  } else {
    pal_low <- rep(pal_default_low, length.out = n_quad)
    pal_high <- rep(pal_default_high, length.out = n_quad)
  }

  .add_quadrant_layers <- function(p, pack, idx, low_pal, high_pal) {
    tile <- pack$tile
    diag <- pack$diag
    id_lab <- pack$id_lab
    type_lab <- pack$type_lab
    ori_lab <- unique(id_lab$orientation)[1]
    hjust_id <- if (ori_lab %in% c("top_right", "top_left")) 0 else 1
    vjust_id <- if (ori_lab %in% c("top_right", "top_left")) -0.25 else 0.25
    p +
      ggplot2::geom_tile(data = tile, ggplot2::aes(x = x_tile, y = y_tile, fill = Correlation),
                         colour = HeatmapTileColor, linewidth = HeatmapTileSize) +
      ggplot2::geom_text(data = tile, ggplot2::aes(x = x_tile, y = y_tile, label = p_signif), size = HeatmapSigSize) +
      ggplot2::geom_text(data = id_lab, ggplot2::aes(x = x_id, y = y_id, label = ID),
                        size = HeatmapLabelSize, vjust = vjust_id, hjust = hjust_id, angle = HeatmapLabelOrient) +
      ggplot2::geom_text(data = type_lab, ggplot2::aes(x = x_type, y = y_type, label = Type),
                        hjust = type_lab$hjust_type[1], size = HeatmapLabelSize) +
      ggplot2::geom_point(data = diag, ggplot2::aes(x = x_diag, y = y_diag),
                         shape = 21, fill = HeatmapPointFill, size = HeatmapPointSize) +
      ggplot2::scale_fill_gradient2(
        low = low_pal[idx], mid = "#ffffff", high = high_pal[idx],
        midpoint = 0, name = paste0("Env ", idx),
        guide = ggplot2::guide_colorbar(order = idx, direction = "horizontal")
      )
  }

  p_heatmaps <- ggplot2::ggplot()
  for (i in seq_along(packs)) {
    if (i > 1) p_heatmaps <- p_heatmaps + ggnewscale::new_scale_fill()
    p_heatmaps <- .add_quadrant_layers(p_heatmaps, packs[[i]], idx = i, pal_low, pal_high)
  }

  node_plot_df <- ggplot_data[[1]]
  edge_plot_df <- ggplot_data[[2]]
  fill_scale_net <- if (!is.null(ggnetview_args$node_fill_values)) {
    ggplot2::scale_fill_manual(values = ggnetview_args$node_fill_values, guide = "none")
  } else {
    scale_fill_ggnetview(unique(node_plot_df$Modularity), guide = "none")
  }

  add_outer <- isTRUE(ggnetview_args$module_outline)
  add_group_outer <- isTRUE(ggnetview_args$network_outline)
  label_arg <- if ("module_label" %in% names(ggnetview_args)) ggnetview_args$module_label else TRUE
  show_module_label <- isTRUE(label_arg) || (is.character(label_arg) && length(label_arg) == 1 && nchar(trimws(label_arg)) > 0)
  module_label_prefix <- if (is.character(label_arg) && length(label_arg) == 1) trimws(label_arg) else "Modularity"
  if (identical(module_label_prefix, "")) module_label_prefix <- "Modularity"
  labelsize <- if (!is.null(ggnetview_args$module_label_size)) as.numeric(ggnetview_args$module_label_size) else 10
  labelsegmentsize <- if (!is.null(ggnetview_args$module_label_segment_width)) as.numeric(ggnetview_args$module_label_segment_width) else 1
  labelsegmentalpha <- if (!is.null(ggnetview_args$module_label_segment_alpha)) as.numeric(ggnetview_args$module_label_segment_alpha) else 1
  mod_col_lab <- if ("modularity3" %in% colnames(graph_ly_scaled)) "modularity3" else "Modularity"
  lab_df <- NULL
  if (show_module_label) {
    xr_lab <- range(graph_ly_scaled$x, na.rm = TRUE)
    yr_lab <- range(graph_ly_scaled$y, na.rm = TRUE)
    x_mid_lab <- stats::median(graph_ly_scaled$x)
    dx_lab <- diff(xr_lab) * 0.12
    lab_df <- graph_ly_scaled %>%
      dplyr::distinct(.data[[mod_col_lab]], .keep_all = TRUE) %>%
      dplyr::filter(.data[[mod_col_lab]] != "Others") %>%
      dplyr::mutate(side = ifelse(x < x_mid_lab, "left", "right")) %>%
      dplyr::group_by(side) %>%
      dplyr::arrange(y, .by_group = TRUE) %>%
      dplyr::mutate(
        y_rank = dplyr::row_number(),
        y_target = scales::rescale(y_rank, to = yr_lab),
        x_anchor = dplyr::if_else(side == "left", xr_lab[1] - dx_lab, xr_lab[2] + dx_lab),
        nudge_x = x_anchor - x,
        nudge_y = y_target - y,
        hjust = dplyr::if_else(side == "left", 1, 0),
        .label_text = paste0(module_label_prefix, " ", .data[[mod_col_lab]])
      ) %>%
      dplyr::ungroup()
  }
  module_label_fun <- function(x) paste0(module_label_prefix, " ", x)
  q_outer <- if (!is.null(ggnetview_args$module_outline_q)) as.numeric(ggnetview_args$module_outline_q) else 0.88
  expand_outer <- if (!is.null(ggnetview_args$module_outline_expand)) as.numeric(ggnetview_args$module_outline_expand) else 1.02
  bandwidth_scale <- if (!is.null(ggnetview_args$module_outline_bandwidth)) as.numeric(ggnetview_args$module_outline_bandwidth) else 1.0
  outerwidth <- if (!is.null(ggnetview_args$module_outline_width)) as.numeric(ggnetview_args$module_outline_width) else 1.25
  outerlinetype <- if (!is.null(ggnetview_args$module_outline_linetype)) ggnetview_args$module_outline_linetype else 2
  outeralpha <- if (!is.null(ggnetview_args$module_outline_alpha)) as.numeric(ggnetview_args$module_outline_alpha) else 0.5

  p0 <- p_heatmaps +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color()
  p0 <- p0 +
    ggplot2::geom_segment(
      data = edge_plot_df,
      ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
      alpha = 0.25,
      colour = "grey70"
    )

  if (add_group_outer && nrow(node_plot_df) > 0) {
    group_circle_df <- node_plot_df %>% dplyr::mutate(.group_outer = 1L)
    circle_n_grp <- max(40, min(300, as.integer(round(8 * sqrt(nrow(group_circle_df))))))
    fill_grp <- if (!is.null(ggnetview_args$network_outline_fill) && length(ggnetview_args$network_outline_fill) > 0)
      ggnetview_args$network_outline_fill[1L] else NA
    alpha_grp <- if (is.na(fill_grp)) 1 else
      if (!is.null(ggnetview_args$network_outline_fill_alpha)) ggnetview_args$network_outline_fill_alpha else 0.2
    p0 <- p0 +
      ggforce::geom_mark_circle(
        data = group_circle_df,
        mapping = ggplot2::aes(x = x, y = y, group = .group_outer),
        fill = fill_grp,
        alpha = alpha_grp,
        color = if (!is.null(ggnetview_args$network_outline_color)) ggnetview_args$network_outline_color else "grey50",
        linetype = if (!is.null(ggnetview_args$network_outline_linetype)) ggnetview_args$network_outline_linetype else 1,
        linewidth = if (!is.null(ggnetview_args$network_outline_width)) ggnetview_args$network_outline_width else 0.5,
        n = circle_n_grp,
        expand = grid::unit(if (!is.null(ggnetview_args$network_outline_expand)) ggnetview_args$network_outline_expand else 2, "mm")
      )
  }

  if (add_outer) {
    maskTable <- generateMask_ggnetview(
      dims = graph_ly_scaled[, c("x", "y")],
      clusters = graph_ly_scaled$Modularity,
      q = q_outer,
      expand = expand_outer,
      bandwidth_scale = bandwidth_scale
    )
    if (nrow(maskTable) > 0L) {
      maskTable <- maskTable %>%
        dplyr::mutate(cluster = factor(cluster, levels = levels(graph_ly_scaled$Modularity), ordered = TRUE))
      mask_classes <- levels(maskTable$cluster)
      fill_scale_mask <- if (!is.null(ggnetview_args$node_fill_values)) {
        ggplot2::scale_fill_manual(values = ggnetview_args$node_fill_values, guide = "none")
      } else {
        scale_fill_ggnetview(mask_classes, guide = "none")
      }
      color_scale_mask <- if (!is.null(ggnetview_args$node_fill_values)) {
        ggplot2::scale_color_manual(values = ggnetview_args$node_fill_values, guide = "none")
      } else {
        scale_color_ggnetview(mask_classes, guide = "none")
      }
      p0 <- p0 +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        ggplot2::geom_polygon(
          data = maskTable %>% dplyr::filter(cluster != "Others"),
          mapping = ggplot2::aes(x = x, y = y,
                                 group = interaction(cluster, polygon_id),
                                 fill = cluster, color = cluster),
          linewidth = outerwidth,
          linetype = outerlinetype,
          alpha = outeralpha,
          show.legend = FALSE
        ) +
        fill_scale_mask +
        color_scale_mask +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color()
    }
  }

  p0 <- p0 +
    ggplot2::geom_point(
      data = node_plot_df,
      ggplot2::aes(x = x, y = y, fill = Modularity),
      shape = 21,
      size = 3,
      alpha = 0.9
    ) +
    fill_scale_net

  p1 <- p0 +
    ggnewscale::new_scale_color() +
    ggplot2::geom_segment(
      data = link_df,
      ggplot2::aes(x = x, y = y, xend = x_to, yend = y_to,
                   color = Correlation, linetype = line_type, linewidth = sig_strength),
      alpha = SigLineAlpha
    ) +
    ggplot2::scale_color_gradient(
      low = SigLineColor[1], high = SigLineColor[2],
      guide = ggplot2::guide_colorbar(direction = "horizontal")
    ) +
    ggplot2::scale_linewidth_continuous(
      range = SigLineWidth, name = "-log10(Pvalue)",
      guide = ggplot2::guide_legend(nrow = 1)  # horizontal row at top
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::geom_point(data = diag_all, ggplot2::aes(x = x_diag, y = y_diag),
                        shape = 21, fill = HeatmapPointFill, size = HeatmapPointSize) +
    {       if (show_module_label && !is.null(lab_df) && nrow(lab_df) > 0) {
        lab_classes <- sort(unique(lab_df$Modularity))
        color_scale_lab <- if (!is.null(ggnetview_args$node_fill_values)) ggplot2::scale_color_manual(values = ggnetview_args$node_fill_values, labels = module_label_fun) else scale_color_ggnetview(lab_classes, labels = module_label_fun)
        list(
          ggnewscale::new_scale_color(),
          ggrepel::geom_text_repel(
            data = lab_df,
            mapping = ggplot2::aes(x = x, y = y, label = .label_text, color = Modularity),
            size = labelsize,
            nudge_x = lab_df$nudge_x,
            nudge_y = lab_df$nudge_y,
            hjust = lab_df$hjust,
            min.segment.length = 0,
            segment.size = labelsegmentsize,
            segment.alpha = labelsegmentalpha,
            max.overlaps = Inf,
            box.padding = 0.15,
            point.padding = 0.15,
            force = 0.05,
            show.legend = FALSE
          ),
          color_scale_lab
        )
      } else list() } +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      aspect.ratio = 1,
      legend.position = "top",
      legend.box = "horizontal",
      legend.box.just = "left"
    )

  p2 <- p0 +
    ggnewscale::new_scale_color() +
    ggplot2::geom_curve(
      data = link_df,
      ggplot2::aes(x = x, y = y, xend = x_to, yend = y_to,
                   color = Correlation, linetype = line_type, linewidth = sig_strength),
      alpha = SigLineAlpha,
      curvature = 0.25
    ) +
    ggplot2::scale_color_gradient(
      low = SigLineColor[1], high = SigLineColor[2],
      guide = ggplot2::guide_colorbar(direction = "horizontal")
    ) +
    ggplot2::scale_linewidth_continuous(
      range = SigLineWidth, name = "-log10(Pvalue)",
      guide = ggplot2::guide_legend(nrow = 1)  # horizontal row at top
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::geom_point(data = diag_all, ggplot2::aes(x = x_diag, y = y_diag),
                        shape = 21, fill = HeatmapPointFill, size = HeatmapPointSize) +
    {       if (show_module_label && !is.null(lab_df) && nrow(lab_df) > 0) {
        lab_classes <- sort(unique(lab_df$Modularity))
        color_scale_lab <- if (!is.null(ggnetview_args$node_fill_values)) ggplot2::scale_color_manual(values = ggnetview_args$node_fill_values, labels = module_label_fun) else scale_color_ggnetview(lab_classes, labels = module_label_fun)
        list(
          ggnewscale::new_scale_color(),
          ggrepel::geom_text_repel(
            data = lab_df,
            mapping = ggplot2::aes(x = x, y = y, label = .label_text, color = Modularity),
            size = labelsize,
            nudge_x = lab_df$nudge_x,
            nudge_y = lab_df$nudge_y,
            hjust = lab_df$hjust,
            min.segment.length = 0,
            segment.size = labelsegmentsize,
            segment.alpha = labelsegmentalpha,
            max.overlaps = Inf,
            box.padding = 0.15,
            point.padding = 0.15,
            force = 0.05,
            show.legend = FALSE
          ),
          color_scale_lab
        )
      } else list() } +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      aspect.ratio = 1,
      legend.position = "top",
      legend.box = "horizontal",
      legend.box.just = "left"
    )

  list(p1, p2, cor_spec_env_list_out)
}
