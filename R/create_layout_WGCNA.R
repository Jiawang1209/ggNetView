#' WGCNA-style bubble-pack layout (local FR + front-chain packing)
#'
#' Produces a layout that reproduces the classic WGCNA "ME-bubbles"
#' visualization (e.g. Nature Medicine 2024, Fig. 1): each module is drawn
#' as a visually distinct, organic-looking cluster of nodes, and the
#' clusters themselves are packed tightly into a single roughly-circular
#' arrangement with module area proportional to module size.
#'
#' The reason a plain \code{layout = "fr"} does NOT reproduce that style:
#' a global Fruchterman-Reingold run lets between-module edges pull
#' modules toward one another, so the separate bubbles collapse into a
#' single blob.  This layout sidesteps that by laying out each module
#' independently and then placing the pre-laid-out clusters geometrically.
#'
#' Algorithm:
#' \enumerate{
#'   \item Compute \code{mod_levels} (biggest module first, "Others" last)
#'     using exactly the same rule as \code{module_layout4()}, so the
#'     coordinates returned here line up row-by-row with the sorted node
#'     order produced in the second tier.
#'   \item Assign each module a disc radius \code{r_mod = sqrt(n_mod) * r}
#'     (area \eqn{\propto} node count).
#'   \item Pack all module discs via progressive front-chain circle packing
#'     (Wang et al. 2006 variant): circle 1 at origin, circle 2 tangent on
#'     +x axis, circle i>=3 tangent to two consecutive circles on the
#'     current front chain, choosing the non-overlapping candidate closest
#'     to the origin.  This yields a compact, roughly-circular outer shape.
#'   \item For each module, extract the induced sub-graph (only within-module
#'     edges), run a local Fruchterman-Reingold layout, then centre the
#'     coordinates and rescale so every node fits inside the allocated disc
#'     radius (with a small margin).  Modules with no internal edges fall
#'     back to uniform-random disc fill.
#'   \item Apply the standard \code{orientation} / \code{angle} rotation
#'     that every other \code{create_layout_*} function uses.
#' }
#'
#' Pair this layout with \code{layout_module = "order"} (which dispatches
#' to \code{module_layout4()}) so the module grouping established here is
#' preserved during node-to-slot assignment.  Using \code{"random"} or
#' \code{"adjacent"} will re-shuffle the coordinates via kNN and destroy
#' the bubble-pack structure.
#'
#' @noRd
create_layout_WGCNA <- function(
    graph_obj,
    node_add = NULL,
    r = 1,
    inner_shrink = 1,
    scale = TRUE,
    anchor_dist = 10,
    orientation = c("up", "down", "left", "right"),
    angle = 0
){

  orientation <- match.arg(orientation)
  base_angle <- switch(orientation,
                       up = 0, right = -pi/2, down = pi, left = pi/2)
  theta_shift <- base_angle + .normalize_angle(angle)

  # --- Validate inner_shrink ---
  # Default 1 reproduces the historical behaviour exactly (fit_radius = 0.95 * radius).
  # Values < 1 contract intra-module FR/uniform fill toward each module's centre
  # *without* changing the packed disc geometry, so inter-module gaps appear
  # while disc centre-to-centre distances stay invariant.
  if (!is.numeric(inner_shrink) || length(inner_shrink) != 1L ||
      !is.finite(inner_shrink) || inner_shrink <= 0) {
    stop("`inner_shrink` must be a single positive finite numeric value.")
  }

  # --- Extract node data ---
  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  if (!"modularity3" %in% colnames(node_df)) {
    stop("WGCNA layout requires a `modularity3` column on the graph. ",
         "Use build_graph_from_wgcna() or build_graph_from_mat() to create the graph.")
  }

  # --- mod_levels: biggest module first, "Others" last ---
  # This must match module_layout4()'s computation exactly, otherwise the
  # coordinates returned here will not align row-by-row with the sorted
  # node order bound in the second tier.
  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    dplyr::mutate(modularity4 = factor(
      modularity3,
      levels = c(setdiff(modularity3, "Others"), "Others"),
      ordered = TRUE
    )) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels

  # --- Per-module counts in mod_levels order ---
  mod_counts_raw <- node_df %>%
    dplyr::count(modularity3, name = "n") %>%
    dplyr::mutate(modularity3 = as.character(modularity3))

  mod_counts <- mod_counts_raw[match(mod_levels, mod_counts_raw$modularity3), , drop = FALSE]
  rownames(mod_counts) <- NULL

  n_mod <- nrow(mod_counts)
  if (n_mod < 1L) {
    stop("WGCNA layout requires at least 1 module (from column `modularity3`).")
  }

  # --- Module radii: area proportional to node count ---
  # `r` is a user-tunable scale factor. Larger r -> looser packing.
  mod_counts$radius <- sqrt(mod_counts$n) * r

  # --- Pack module discs using progressive front-chain packing ---
  packing <- .pack_wgcna_modules(mod_counts$radius)
  mod_counts$cx <- packing$cx
  mod_counts$cy <- packing$cy

  # --- Build a full igraph once, then induce subgraphs per module ---
  # Using a shared igraph avoids re-coercing graph_obj inside the loop.
  ig_full <- igraph::as.igraph(graph_obj)

  # `modularity3` on vertices is used to select subgraph membership.
  vattr_mod <- igraph::vertex_attr(ig_full, "modularity3")
  if (is.null(vattr_mod)) {
    stop("Internal error: igraph vertices lost the `modularity3` attribute.")
  }

  # --- Lay out each module independently: local FR, then fit into disc ---
  # set.seed() has already been called by ggNetView(), so both the FR
  # random initialization and the uniform-fill fallback are deterministic.
  ly_list <- lapply(seq_len(n_mod), function(i) {
    mod_name <- mod_levels[i]
    keep     <- which(as.character(vattr_mod) == mod_name)

    .fill_disc_fr(
      ig_full = ig_full,
      keep    = keep,
      cx      = mod_counts$cx[i],
      cy      = mod_counts$cy[i],
      radius  = mod_counts$radius[i],
      inner_shrink = inner_shrink
    )
  })

  ly <- dplyr::bind_rows(ly_list)

  # --- Orientation rotation (same idiom as every other create_layout_*) ---
  if (theta_shift != 0) {
    Rm <- matrix(c(cos(theta_shift), -sin(theta_shift),
                   sin(theta_shift),  cos(theta_shift)), nrow = 2)
    xy <- as.matrix(ly[, c("x", "y")])
    ly[, c("x", "y")] <- t(Rm %*% t(xy))
  }

  rownames(ly) <- NULL
  ly
}


# ---------------------------------------------------------------------------
# Internal helpers (not exported, not documented as Rd pages)
# ---------------------------------------------------------------------------

#' Progressive front-chain circle packing
#'
#' Greedy placement that keeps the arrangement compact. Circle 1 is placed
#' at the origin; circle 2 is placed tangent to circle 1 along the +x axis;
#' each subsequent circle is placed tangent to two consecutive circles on
#' the current front chain, choosing the candidate that keeps the
#' \strong{minimum bounding radius measured from the area-weighted
#' centroid of the packing} as small as possible, among all
#' non-overlapping tangent positions.  After the loop finishes, the
#' whole packing is translated so its area-weighted centroid sits at
#' the origin.
#'
#' Why \strong{area-weighted centroid}, not origin?  Scoring from the
#' origin looks fine for uniform-sized circles but gets visibly biased
#' once radii are heterogeneous -- circle 0 is anchored at (0, 0) and
#' circle 1 is anchored on the +x axis, so the true centre of mass
#' drifts sideways from step 2 onward.  Minimising \code{cx^2 + cy^2}
#' then pulls subsequent placements toward a point that is no longer
#' the visual centre, producing a comma / crescent silhouette.
#' Minimising the distance \emph{from the current centroid} directly
#' optimises the property we care about -- a circular outer silhouette
#' -- and making that centroid area-weighted prevents a cluster of many
#' tiny discs from out-voting a single dominant module.  The final
#' translation then makes the plotted graph visually centred.
#'
#' @param radii numeric vector of module radii, ordered largest first.
#' @return data.frame with columns \code{cx}, \code{cy}, one row per input radius.
#' @noRd
.pack_wgcna_modules <- function(radii) {
  n <- length(radii)
  if (n == 0L) return(data.frame(cx = numeric(0), cy = numeric(0)))

  cx <- numeric(n)
  cy <- numeric(n)

  cx[1] <- 0; cy[1] <- 0
  if (n == 1L) return(data.frame(cx = cx, cy = cy))

  cx[2] <- radii[1] + radii[2]; cy[2] <- 0
  if (n == 2L) {
    # Still translate to area-weighted centroid for consistency.
    w <- radii^2
    cxbar <- sum(w * cx) / sum(w)
    cybar <- sum(w * cy) / sum(w)
    return(data.frame(cx = cx - cxbar, cy = cy - cybar))
  }

  # Cyclic list of indices currently on the outer frontier
  front <- c(1L, 2L)

  for (i in 3:n) {
    r_i <- radii[i]
    best_cx <- NA_real_
    best_cy <- NA_real_
    best_score <- Inf
    best_k <- NA_integer_

    # Pre-computed sums for incremental centroid of the first i circles.
    # After placing circle i at (cx_p, cy_p) with radius r_i, the
    # area-weighted centroid (weights = r^2) is:
    #   cxbar = (sum_{q<i} r_q^2 * cx[q] + r_i^2 * cx_p) / W_new
    #   W_new  = sum_{q<i} r_q^2 + r_i^2
    # We then score by the maximum (dist-to-centroid + radius) over all
    # i circles -- i.e. the bounding-radius of the whole packing
    # measured from the true centre of mass, not from the origin.
    r2_prev   <- radii[1:(i - 1)]^2
    W_prev    <- sum(r2_prev)
    Sx_prev   <- sum(r2_prev * cx[1:(i - 1)])
    Sy_prev   <- sum(r2_prev * cy[1:(i - 1)])

    m <- length(front)
    for (k in seq_len(m)) {
      j1 <- front[k]
      j2 <- front[(k %% m) + 1L]

      cand <- .tangent_to_two(
        cx[j1], cy[j1], radii[j1],
        cx[j2], cy[j2], radii[j2],
        r_i
      )
      if (is.null(cand)) next

      # Try both tangent solutions
      for (p in 1:2) {
        cx_p <- cand$cx[p]
        cy_p <- cand$cy[p]

        # Reject if overlapping any already-placed circle
        dd <- (cx[1:(i - 1)] - cx_p)^2 + (cy[1:(i - 1)] - cy_p)^2
        min_allowed <- (radii[1:(i - 1)] + r_i)^2 - 1e-6
        if (any(dd < min_allowed)) next

        # Hypothetical area-weighted centroid with circle i included
        W_new <- W_prev + r_i * r_i
        cxbar <- (Sx_prev + r_i * r_i * cx_p) / W_new
        cybar <- (Sy_prev + r_i * r_i * cy_p) / W_new

        # Bounding radius of the full packing measured from that centroid
        d_prev <- sqrt((cx[1:(i - 1)] - cxbar)^2 + (cy[1:(i - 1)] - cybar)^2) +
                  radii[1:(i - 1)]
        d_new  <- sqrt((cx_p - cxbar)^2 + (cy_p - cybar)^2) + r_i
        score  <- max(max(d_prev), d_new)

        if (score < best_score) {
          best_score <- score
          best_cx <- cx_p
          best_cy <- cy_p
          best_k <- k
        }
      }
    }

    if (is.na(best_cx)) {
      # No valid tangent placement found on current front chain;
      # fall back to pushing the new circle onto the +x perimeter.
      existing_outer <- max(
        sqrt(cx[1:(i - 1)]^2 + cy[1:(i - 1)]^2) + radii[1:(i - 1)]
      )
      best_cx <- existing_outer + r_i
      best_cy <- 0
      best_k <- NA_integer_
    }

    cx[i] <- best_cx
    cy[i] <- best_cy

    if (!is.na(best_k)) {
      # Insert i between front[best_k] and front[best_k + 1]
      front <- append(front, i, after = best_k)
    } else {
      front <- c(front, i)
    }
  }

  # Translate the whole packing so its area-weighted centroid is the
  # origin.  Without this, the fixed anchor at (0, 0) + (d, 0) would
  # leave the plot visually off-centre for every non-trivial radii vector.
  w     <- radii^2
  cxbar <- sum(w * cx) / sum(w)
  cybar <- sum(w * cy) / sum(w)

  data.frame(cx = cx - cxbar, cy = cy - cybar)
}


#' Apollonius tangent construction
#'
#' Given two circles (x1,y1,r1) and (x2,y2,r2), return the two centers
#' where a new circle of radius r3 can sit tangent to both. Uses the law
#' of cosines on the triangle with sides (r1+r3), (r2+r3), d(c1,c2).
#'
#' @return list with numeric vectors \code{cx}, \code{cy} of length 2, or
#'   \code{NULL} when no tangent solution exists.
#' @noRd
.tangent_to_two <- function(x1, y1, r1, x2, y2, r2, r3) {
  d12 <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  a <- r1 + r3
  b <- r2 + r3
  c <- d12

  if (c > a + b) return(NULL)       # circles too far apart
  if (c < abs(a - b)) return(NULL)  # one encloses the other
  if (c <= 0) return(NULL)

  cos_alpha <- (a * a + c * c - b * b) / (2 * a * c)
  cos_alpha <- max(-1, min(1, cos_alpha))
  alpha <- acos(cos_alpha)

  base <- atan2(y2 - y1, x2 - x1)

  list(
    cx = c(x1 + a * cos(base + alpha),
           x1 + a * cos(base - alpha)),
    cy = c(y1 + a * sin(base + alpha),
           y1 + a * sin(base - alpha))
  )
}


#' Run Fruchterman-Reingold on a single module and fit the result into a disc
#'
#' Extracts the induced subgraph on the given vertex set (keeping only
#' within-module edges), then places every vertex inside a disc of
#' radius \code{radius} centred at \code{(cx, cy)}.
#'
#' \strong{Handling of fragmented modules.}  When the induced subgraph is
#' not fully connected (common after aggressive sparsification, e.g.
#' \code{top_k = 1} in \code{trans_TOM_in_WGCNA()}), running FR on the
#' union of all components produces the classic artefact of a couple of
#' tight sub-clusters surrounded by a spray of isolated dots pushed
#' toward the boundary.  We avoid that by:
#' \enumerate{
#'   \item Running FR only on the \strong{largest connected component}
#'     (LCC), scaled into a disc of radius \code{radius * 0.95}.
#'   \item Placing every non-LCC vertex (smaller components and isolates)
#'     by uniform-disc fill within the \strong{same} inner disc.  These
#'     vertices have no FR force to position them informatively, so
#'     blending them into the same area as the LCC gives a cohesive
#'     bubble instead of a ring of orphans.
#' }
#' A 5\% margin is left between the outermost point and the packing
#' boundary so nodes don't graze the disc edge.
#'
#' If the induced subgraph has no edges at all, the function falls back
#' to a single uniform disc fill.  Single nodes are placed at the disc
#' centre; empty modules return 0 rows.
#'
#' \strong{Intra-module compactness (\code{inner_shrink}).}  The effective
#' fill radius is \code{radius * 0.95 * inner_shrink}.  At the default
#' \code{inner_shrink = 1} the behaviour is identical to the original
#' implementation -- nodes spread to fill 95\% of the allocated disc.
#' Smaller values (e.g. 0.65) tighten the FR cluster toward each module's
#' centre, exposing the natural hub/periphery structure produced by FR
#' and leaving inter-module whitespace.  The packed disc geometry
#' (\code{cx, cy, radius}) is unchanged, so module centre-to-centre
#' distances are invariant under \code{inner_shrink}; only the visible
#' point cloud per module shrinks.
#'
#' Vertex ordering is preserved: the rows of the returned data.frame are
#' in the same order as \code{keep}, which matches the order produced by
#' \code{module_layout4()} for this module.
#'
#' @noRd
.fill_disc_fr <- function(ig_full, keep, cx, cy, radius, inner_shrink = 1) {
  n <- length(keep)
  if (n <= 0L) return(data.frame(x = numeric(0), y = numeric(0)))
  if (n == 1L) return(data.frame(x = cx, y = cy))

  sub <- igraph::induced_subgraph(ig_full, vids = keep)

  # No internal edges -> FR degenerates; use uniform disc fill.
  # Historical behaviour (preserved at inner_shrink = 1): the empty-edge
  # branch uses the full disc `radius`, *without* the 0.95 inner margin
  # that the with-edge branch applies.  We multiply by `inner_shrink`
  # so this branch contracts proportionally when the user opts in.
  if (igraph::ecount(sub) == 0L) {
    return(.fill_disc_uniform(n = n, cx = cx, cy = cy,
                              radius = radius * inner_shrink))
  }

  # Inner radius shared by FR scaling and the non-LCC uniform fill.
  # The 0.95 factor keeps a small margin between the outermost point and
  # the disc boundary; `inner_shrink` (default 1) further contracts the
  # fill area toward (cx, cy) without altering the disc-packing geometry.
  fit_radius <- radius * 0.95 * inner_shrink

  # Split vertices into "belongs to largest connected component" vs "rest"
  # (smaller components + isolates).  FR will be run only on the LCC.
  comps  <- igraph::components(sub)
  lcc_id <- which.max(comps$csize)
  mem    <- comps$membership
  is_lcc <- (mem == lcc_id)

  lcc_idx  <- which(is_lcc)     # positions within `keep` / within `sub`
  rest_idx <- which(!is_lcc)

  # Pre-allocate per-vertex coordinates (indexed in the order of `keep`).
  x <- numeric(n)
  y <- numeric(n)

  # --- LCC: local FR, centred, scaled into disc ---
  if (length(lcc_idx) == 1L) {
    # Degenerate LCC of size 1; just drop it at the centre.
    x[lcc_idx] <- cx
    y[lcc_idx] <- cy
  } else {
    sub_lcc <- igraph::induced_subgraph(sub, vids = lcc_idx)
    coords  <- igraph::layout_with_fr(sub_lcc, niter = 500L)
    coords[, 1] <- coords[, 1] - mean(coords[, 1])
    coords[, 2] <- coords[, 2] - mean(coords[, 2])
    rmax <- max(sqrt(coords[, 1]^2 + coords[, 2]^2))
    if (rmax > 0) {
      coords <- coords * (fit_radius / rmax)
    }
    x[lcc_idx] <- coords[, 1] + cx
    y[lcc_idx] <- coords[, 2] + cy
  }

  # --- Non-LCC vertices: blend into the SAME inner disc via uniform fill ---
  # They'd otherwise be painted as a halo around the tight FR cluster.
  if (length(rest_idx) > 0L) {
    fill <- .fill_disc_uniform(
      n      = length(rest_idx),
      cx     = cx,
      cy     = cy,
      radius = fit_radius
    )
    x[rest_idx] <- fill$x
    y[rest_idx] <- fill$y
  }

  data.frame(x = x, y = y)
}


#' Uniform random sample of n points inside a disc centred at (cx, cy)
#'
#' Fallback used by \code{.fill_disc_fr()} when a module has no internal
#' edges.  Uses the inverse-CDF transform \code{r = R * sqrt(U)},
#' \code{theta = 2 * pi * V}, which yields a genuinely uniform density
#' (unlike sampling \code{r} uniformly, which concentrates points at the
#' centre).  Relies on the outer \code{set.seed()} call in \code{ggNetView()}
#' for reproducibility.
#'
#' @noRd
.fill_disc_uniform <- function(n, cx, cy, radius) {
  if (n <= 0L) return(data.frame(x = numeric(0), y = numeric(0)))
  if (n == 1L) return(data.frame(x = cx, y = cy))

  u <- stats::runif(n)
  v <- stats::runif(n)
  r <- radius * sqrt(u)
  theta <- 2 * pi * v

  data.frame(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta)
  )
}
