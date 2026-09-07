#' @noRd
get_neighbors <- function(ly,
                          k = 5,
                          idx = NULL,
                          coord = NULL,
                          # seed = seed,
                          tol = 1e-12) {








  stopifnot(all(c("x","y") %in% names(ly)))
  stopifnot(is.numeric(ly$x), is.numeric(ly$y))
  n <- nrow(ly)
  if (n < 1) stop("ly requires at least one point.")

  # set.seed(seed)



  if (!is.null(coord)) {
    fx <- coord[1]
    fy <- coord[2]
    self_ids <- which(abs(ly$x - fx) <= tol & abs(ly$y - fy) <= tol)
    idx_used <- if (length(self_ids) == 1) self_ids else NA_integer_
  } else {


    if (is.null(idx)) idx <- sample.int(n, 1)
    stopifnot(idx >= 1, idx <= n)
    fx <- ly$x[idx]
    fy <- ly$y[idx]
    self_ids <- idx
    idx_used <- idx
  }


  d <- sqrt((ly$x - fx)^2 + (ly$y - fy)^2)


  cand <- which(is.finite(d))
  o <- cand[order(d[cand], cand)]


  self_pick <- intersect(self_ids, o)
  others <- setdiff(o, self_pick)
  need_self <- length(self_pick)
  k_self <- min(need_self, k)
  k_other <- max(0, k - k_self)
  nn_idx <- c(utils::head(self_pick, k_self), utils::head(others, k_other))


  neighbors <- data.frame(
    node = nn_idx,
    x    = ly$x[nn_idx],
    y    = ly$y[nn_idx],
    dist = d[nn_idx]
  )
  focal <- data.frame(
    node = idx_used,
    x = fx,
    y = fy
  )

  if (nrow(neighbors) < k) {
    message(sprintf("Only %d neighbors returned (fewer than the requested k = %d).", nrow(neighbors), k))
  }

  list(focal = focal, neighbors = neighbors)
}



get_neighbors_hub <- function(ly,
                              k=5,
                              idx=NULL,
                              coord=NULL,
                              # seed=seed,
                              tol=1e-12) {
  out <- get_neighbors(ly = ly,
                       k = k,
                       idx = idx,
                       coord = coord,
                       # seed = seed,
                       tol = tol)
  nb  <- out$neighbors

  cx  <- mean(nb$x); cy <- mean(nb$y)
  r   <- sqrt((nb$x - cx)^2 + (nb$y - cy)^2)
  ord <- order(r, nb$dist, nb$node)
  out$neighbors <- nb[ord, , drop = FALSE]
  out
}

shrink_rings_global <- function(df, shrink){
  cx <- mean(df$x); cy <- mean(df$y)
  r  <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
  th <- atan2(df$y - cy, df$x - cx)
  r2 <- r * shrink
  df$x <- cx + r2 * cos(th)
  df$y <- cy + r2 * sin(th)
  df
}

radial_offset <- function(df, delta){
  r <- sqrt(df$x^2 + df$y^2)
  th <- atan2(df$y, df$x)
  df$x <- (r + delta) * cos(th)
  df$y <- (r + delta) * sin(th)
  df
}



module_layout <- function(graph_obj,
                          layout,
                          center = TRUE,
                          idx = NULL,
                          shrink = 1,
                          jitter,
                          jitter_sd
                          # seed = seed
                          ){

  # set.seed(seed)


  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels


  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))


  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)


  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)


  coord <- if (isTRUE(center)) c(0, 0) else NULL

  neighbors_list <- list()



  for (i in seq_len(nrow(node_df_sorted_number))) {

    if (i == 1) {
      out <- get_neighbors_hub(ly = layout,
                           k = node_df_sorted_number$n[i],
                           coord = coord,
                           idx = idx)
      out_ly <- out$neighbors
      coords <- out_ly %>% dplyr::select(x, y)

      coords <- shrink_rings_global(coords, shrink = shrink)

      neighbors_list[[i]] <- coords
      ly_sub <- layout[-out_ly$node, , drop = FALSE]

    } else if (i == nrow(node_df_sorted_number)) {
        # The last iteration dumps whatever is left in `ly_sub`.  Keep
        # only the (x, y) coordinate columns -- otherwise any
        # module-aware first-tier layout that adds an extra column to
        # its return value (e.g. `group` in
        # `create_layout_circular_modules_*_layout()` or any
        # `*partite_*_layout()` family) would produce a 3-column
        # data.frame here, while the earlier iterations always pushed
        # 2-column data.frames into `neighbors_list`.  The subsequent
        # `do.call(rbind, neighbors_list)` would then fail with
        # "numbers of columns of arguments do not match".
        neighbors_list[[i]] <- ly_sub %>% dplyr::select(x, y)
    } else {
      out <- get_neighbors(ly = ly_sub,
                           k = node_df_sorted_number$n[i])
      out_ly <- out$neighbors
      coords <- out_ly %>% dplyr::select(x, y)

      coords <- shrink_rings_global(coords, shrink = shrink)

      neighbors_list[[i]] <- coords
      ly_sub <- ly_sub[-out_ly$node, , drop = FALSE]
    }
  }

  ly_final <- do.call(rbind, neighbors_list)

  if (isTRUE(jitter)) {
    ly_final <- ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }else{
    ly_final <- ly_final
  }

  # combine

  graph_ly_final <- dplyr::bind_cols(
    ly_final,
    graph_obj_sort %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
  )

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)

  # get location result

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort,
              graph_ly_final = graph_ly_final,
              ggplot_data = ggplot_data
              ))
}



module_layout2 <- function(graph_obj,
                           layout,
                           center = TRUE,
                           shrink = 0.9,
                           k_nn = 8,
                           push_others_delta = 0.2,
                           jitter,
                           jitter_sd
                           # seed = seed
                           ) {

  stopifnot(all(c("x","y") %in% names(layout)))

  # set.seed(seed)


  # A graph left with no nodes (e.g. every edge dropped by the correlation /
  # p-value thresholds) cannot be laid out. Say so, instead of failing later
  # with an opaque "Columns `x` and `y` don't exist" from the empty result.
  if (igraph::vcount(tidygraph::as.igraph(graph_obj)) == 0L) {
    stop("`graph_obj` has no nodes left to lay out: every edge was removed by the ",
         "current thresholds. Relax `r.threshold` / `p.threshold`, or drop this ",
         "sub-network from the comparison.", call. = FALSE)
  }

  xy <- as.matrix(layout[, c("x","y")])
  # FNN::get.knn() does NOT error when k >= nrow(xy): it warns and returns
  # out-of-range indices, which downstream turns into an empty region-growing
  # result and the opaque "Columns `x` and `y` don't exist" failure. Clamp k to
  # the number of available neighbours, and handle the degenerate 0/1-slot case.
  n_slot <- nrow(xy)
  if (n_slot < 2L) {
    adj <- rep(list(integer(0)), n_slot)
  } else {
    k_eff <- max(1L, min(as.integer(k_nn), n_slot - 1L))
    nn    <- FNN::get.knn(xy, k = k_eff)$nn.index
    adj   <- lapply(seq_len(n_slot), function(i) unique(stats::na.omit(nn[i, ])))
  }


  nodes_tb <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_order <- nodes_tb %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::pull(modularity4)

  nodes_sorted <- nodes_tb %>%
    dplyr::mutate(modularity3 = factor(modularity3, levels = mod_order)) %>%
    dplyr::arrange(modularity3, dplyr::desc(Degree))

  need_per_mod <- nodes_sorted %>%
    dplyr::count(modularity3)


  nL    <- nrow(layout)
  free  <- rep(TRUE, nL)
  claim <- rep(NA_character_, nL)
  placed_coords <- setNames(vector("list", nrow(need_per_mod)),
                            as.character(need_per_mod$modularity3))


  get_frontier <- function(){
    used_idx <- which(!free)
    if (length(used_idx) == 0) return(integer())
    uniq <- unique(unlist(adj[used_idx], use.names = FALSE))
    uniq[free[uniq]]
  }


  pick_seed_first <- function(){
    if (center) {
      r <- sqrt(layout$x^2 + layout$y^2)
      cand <- order(r)[seq_len(max(10, ceiling(0.01 * length(r))))]
      sample(cand, 1)
    } else {
      sample(which(free), 1)
    }
  }


  pick_seed_next <- function(){
    fr <- get_frontier()
    if (length(fr) == 0) return(sample(which(free), 1))
    cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
    dd <- (layout$x[fr] - cx)^2 + (layout$y[fr] - cy)^2
    prob <- (max(dd) + 1e-9 - dd); prob <- prob / sum(prob)
    sample(fr, 1, prob = prob)
  }


  grow_region <- function(seed_idx, quota){
    q <- seed_idx
    region <- integer()
    while (length(region) < quota && length(q) > 0) {
      i <- q[1]; q <- q[-1]
      if (!free[i]) next
      region <- c(region, i); free[i] <<- FALSE
      nb <- adj[[i]]
      nb <- nb[free[nb]]
      if (length(nb) == 0) next


      if (any(!is.na(claim))) {
        cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
      } else {
        cx <- 0; cy <- 0
      }
      sc  <- -((layout$x[nb]-cx)^2 + (layout$y[nb]-cy)^2)
      # ord <- order(sc, runif(length(nb)))
      ord <- order(sc, nb)
      nb  <- nb[ord]
      q   <- c(q, nb)
    }
    region
  }


  fill_deficit <- function(region, quota, target_idx){
    need <- quota - length(region)
    if (need <= 0) return(region)
    free_idx <- which(free)
    if (length(free_idx) == 0) return(region)

    if (length(region) > 0) {
      tx <- mean(layout$x[region]); ty <- mean(layout$y[region])
    } else {
      tx <- layout$x[target_idx];    ty <- layout$y[target_idx]
    }
    dd <- (layout$x[free_idx] - tx)^2 + (layout$y[free_idx] - ty)^2
    pick <- head(free_idx[order(dd)], min(need, length(free_idx)))
    free[pick] <<- FALSE
    c(region, pick)
  }


  for (mi in seq_len(nrow(need_per_mod))) {
    mod    <- as.character(need_per_mod$modularity3[mi])
    k_need <- need_per_mod$n[mi]
    if (mod == "Others") next

    if (sum(free) < k_need) {
      stop(sprintf("Not enough free layout points: module %s needs %d, but only %d remain.",
                   mod, k_need, sum(free)))
    }

    seed_idx <- if (mi == 1) pick_seed_first() else pick_seed_next()
    region   <- grow_region(seed_idx, k_need)
    if (length(region) < k_need) {
      region <- fill_deficit(region, k_need, target_idx = seed_idx)
    }

    claim[region] <- mod
    coords <- layout[region, c("x","y"), drop = FALSE]
    coords <- shrink_rings_global(coords, shrink = shrink)


    cx <- mean(coords$x); cy <- mean(coords$y)
    rr <- sqrt((coords$x - cx)^2 + (coords$y - cy)^2)
    ord_pts <- order(rr, coords$x, coords$y)

    nodes_this <- nodes_sorted %>%
      dplyr::filter(modularity3 == mod)

    placed_coords[[mod]] <- dplyr::bind_cols(coords[ord_pts, , drop=FALSE], nodes_this)
  }


  if ("Others" %in% names(placed_coords)) {
    n_oth <- need_per_mod$n[need_per_mod$modularity3 == "Others"]
    if (n_oth > 0) {
      idx_free <- which(free)
      if (length(idx_free) < n_oth) {
        stop(sprintf("`Others` group needs %d points, but only %d remain.", n_oth, length(idx_free)))
      }
      r_free <- sqrt(layout$x[idx_free]^2 + layout$y[idx_free]^2)
      take   <- idx_free[order(-r_free)][seq_len(n_oth)]
      free[take] <- FALSE

      coords <- layout[take, c("x","y")]
      coords <- radial_offset(coords, delta = push_others_delta)

      nodes_oth <- nodes_sorted %>%
        dplyr::filter(modularity3 == "Others")

      placed_coords[["Others"]] <- dplyr::bind_cols(coords, nodes_oth)
    }
  }


  graph_ly_final <- dplyr::bind_rows(placed_coords[as.character(need_per_mod$modularity3)])

  if (isTRUE(jitter)) {
    graph_ly_final <- graph_ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }else{
    graph_ly_final <- graph_ly_final
  }

  ly_final <- graph_ly_final[, c("x","y")]


  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels


  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))


  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)


  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)


  list(
    layout = ly_final,
    graph_obj = graph_obj_sort,
    graph_ly_final = graph_ly_final,
    ggplot_data = ggplot_data
  )
}


module_layout3 <- function(graph_obj,
                           layout,                # data.frame(x, y)
                           center = TRUE,
                           shrink = 0.9,
                           k_nn = 12,
                           push_others_delta = 0.2,
                           jitter,
                           jitter_sd
                           # seed = seed
                           ) {

  stopifnot(all(c("x","y") %in% names(layout)))
  # set.seed(seed)


  # shrink_rings_global <- function(df, shrink){
  #   cx <- mean(df$x); cy <- mean(df$y)
  #   r  <- sqrt((df$x - cx)^2 + (df$y - cy)^2)
  #   th <- atan2(df$y - cy, df$x - cx)
  #   r2 <- r * shrink
  #   df$x <- cx + r2 * cos(th)
  #   df$y <- cy + r2 * sin(th)
  #   df
  # }
  # radial_offset <- function(df, delta){
  #   r <- sqrt(df$x^2 + df$y^2)
  #   th <- atan2(df$y, df$x)
  #   df$x <- (r + delta) * cos(th)
  #   df$y <- (r + delta) * sin(th)
  #   df
  # }


  # A graph left with no nodes (e.g. every edge dropped by the correlation /
  # p-value thresholds) cannot be laid out. Say so, instead of failing later
  # with an opaque "Columns `x` and `y` don't exist" from the empty result.
  if (igraph::vcount(tidygraph::as.igraph(graph_obj)) == 0L) {
    stop("`graph_obj` has no nodes left to lay out: every edge was removed by the ",
         "current thresholds. Relax `r.threshold` / `p.threshold`, or drop this ",
         "sub-network from the comparison.", call. = FALSE)
  }

  xy <- as.matrix(layout[, c("x","y")])
  # FNN::get.knn() does NOT error when k >= nrow(xy): it warns and returns
  # out-of-range indices, which downstream turns into an empty region-growing
  # result and the opaque "Columns `x` and `y` don't exist" failure. Clamp k to
  # the number of available neighbours, and handle the degenerate 0/1-slot case.
  n_slot <- nrow(xy)
  if (n_slot < 2L) {
    adj <- rep(list(integer(0)), n_slot)
  } else {
    k_eff <- max(1L, min(as.integer(k_nn), n_slot - 1L))
    nn    <- FNN::get.knn(xy, k = k_eff)$nn.index
    adj   <- lapply(seq_len(n_slot), function(i) unique(stats::na.omit(nn[i, ])))
  }


  nodes_tb <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()

  mod_order <- nodes_tb %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::pull(modularity4)

  nodes_sorted <- nodes_tb %>%
    dplyr::mutate(modularity3 = factor(modularity3, levels = mod_order)) %>%
    dplyr::arrange(modularity3, dplyr::desc(Degree))

  need_per_mod <- nodes_sorted %>%
    dplyr::count(modularity3)


  nL    <- nrow(layout)
  free  <- rep(TRUE, nL)
  claim <- rep(NA_character_, nL)
  placed_coords <- setNames(vector("list", nrow(need_per_mod)),
                            as.character(need_per_mod$modularity3))


  get_frontier <- function(){
    used_idx <- which(!free)
    if (length(used_idx) == 0) return(integer())
    uniq <- unique(unlist(adj[used_idx], use.names = FALSE))
    uniq[free[uniq]]
  }


  get_free_components <- function(){
    visited <- rep(FALSE, nL)
    comps <- list()
    for (i in which(free)) {
      if (visited[i]) next
      q <- i
      comp <- integer()
      while (length(q) > 0) {
        j <- q[1]; q <- q[-1]
        if (visited[j]) next
        visited[j] <- TRUE
        if (!free[j]) next
        comp <- c(comp, j)
        nb <- adj[[j]]
        nb <- nb[free[nb] & !visited[nb]]
        q <- c(q, nb)
      }
      if (length(comp) > 0) comps <- c(comps, list(comp))
    }
    comps
  }


  pick_seed_first <- function(k_need){
    comps <- get_free_components()
    comps <- comps[order(-sapply(comps, length))]
    ok <- comps[sapply(comps, length) >= k_need]
    if (length(ok) == 0) {
      stop(sprintf("No region with >= %d consecutive slots exists. Please increase k_nn (currently %d) or switch to another layout.", k_need, k_nn))
    }
    cand <- ok[[1]]
    if (center) {
      r <- sqrt(layout$x[cand]^2 + layout$y[cand]^2)
      cand[order(r)[1]]
    } else {
      cand[sample(length(cand), 1)]
    }
  }


  pick_seed_next <- function(k_need){
    comps <- get_free_components()
    ok <- comps[sapply(comps, length) >= k_need]
    if (length(ok) == 0) {
      stop(sprintf("No region with >= %d consecutive slots exists. Please increase k_nn (currently %d) or switch to another layout.", k_need, k_nn))
    }
    cx <- mean(layout$x[!free], na.rm = TRUE)
    cy <- mean(layout$y[!free], na.rm = TRUE)
    best <- ok[[1]]
    best_d <- Inf
    for (cc in ok) {
      dc <- (mean(layout$x[cc]) - cx)^2 + (mean(layout$y[cc]) - cy)^2
      if (dc < best_d) { best_d <- dc; best <- cc }
    }
    fr <- intersect(best, get_frontier())
    if (length(fr) > 0) {
      dd <- (layout$x[fr] - cx)^2 + (layout$y[fr] - cy)^2
      fr[order(dd)[1]]
    } else {
      dd <- (layout$x[best] - cx)^2 + (layout$y[best] - cy)^2
      best[order(dd)[1]]
    }
  }


  grow_region <- function(seed_idx, quota){
    q <- seed_idx
    region <- integer()
    while (length(region) < quota && length(q) > 0) {
      i <- q[1]; q <- q[-1]
      if (!free[i]) next
      region <- c(region, i); free[i] <<- FALSE
      nb <- adj[[i]]
      nb <- nb[free[nb]]
      if (length(nb) == 0) next


      if (any(!is.na(claim))) {
        cx <- mean(layout$x[!free], na.rm = TRUE); cy <- mean(layout$y[!free], na.rm = TRUE)
      } else {
        cx <- 0; cy <- 0
      }
      sc  <- -((layout$x[nb]-cx)^2 + (layout$y[nb]-cy)^2)
      # ord <- order(sc, runif(length(nb)))
      ord <- order(sc, nb)
      nb  <- nb[ord]
      q   <- c(q, nb)
    }
    region
  }


  fill_deficit <- function(region, quota, target_idx){
    need <- quota - length(region)
    if (need <= 0) return(region)
    tx <- if (length(region) > 0) mean(layout$x[region]) else layout$x[target_idx]
    ty <- if (length(region) > 0) mean(layout$y[region]) else layout$y[target_idx]

    while (length(region) < quota) {

      frontier <- unique(unlist(adj[region], use.names = FALSE))
      frontier <- frontier[free[frontier]]

      if (length(frontier) == 0) {
        stop(sprintf(
          paste0(
            "The module requires %d consecutive slots, but the adjacent region can only grow to %d. ",
            "Please increase k_nn (currently %d) or switch to a layout with better connectivity."
          ),
          quota, length(region), k_nn
        ))
      }

      dd <- (layout$x[frontier] - tx)^2 + (layout$y[frontier] - ty)^2
      pick <- frontier[order(dd, frontier)[1]]
      region <- c(region, pick)
      free[pick] <<- FALSE
      tx <- mean(layout$x[region])
      ty <- mean(layout$y[region])
    }
    region
  }


  for (mi in seq_len(nrow(need_per_mod))) {
    mod    <- as.character(need_per_mod$modularity3[mi])
    k_need <- need_per_mod$n[mi]
    if (mod == "Others") next

    if (sum(free) < k_need) {
      stop(sprintf("Not enough free layout points: module %s needs %d, but only %d remain.",
                   mod, k_need, sum(free)))
    }

    seed_idx <- if (mi == 1) pick_seed_first(k_need) else pick_seed_next(k_need)
    region   <- grow_region(seed_idx, k_need)
    if (length(region) < k_need) {
      region <- fill_deficit(region, k_need, target_idx = seed_idx)
    }

    claim[region] <- mod
    coords <- layout[region, c("x","y"), drop = FALSE]
    coords <- shrink_rings_global(coords, shrink = shrink)


    cx <- mean(coords$x); cy <- mean(coords$y)
    rr <- sqrt((coords$x - cx)^2 + (coords$y - cy)^2)
    ord_pts <- order(rr, coords$x, coords$y)

    nodes_this <- nodes_sorted %>%
      dplyr::filter(modularity3 == mod)

    placed_coords[[mod]] <- dplyr::bind_cols(coords[ord_pts, , drop=FALSE], nodes_this)
  }


  if ("Others" %in% names(placed_coords)) {
    n_oth <- need_per_mod$n[need_per_mod$modularity3 == "Others"]
    if (n_oth > 0) {
      idx_free <- which(free)
      if (length(idx_free) < n_oth) {
        stop(sprintf("`Others` group needs %d points, but only %d remain.", n_oth, length(idx_free)))
      }

      r_free <- sqrt(layout$x[idx_free]^2 + layout$y[idx_free]^2)
      take   <- idx_free[order(-r_free)][seq_len(n_oth)]
      free[take] <- FALSE
      coords <- layout[take, c("x","y"), drop = FALSE]
      coords <- radial_offset(coords, delta = push_others_delta)



      nodes_oth <- nodes_sorted %>%
        dplyr::filter(modularity3 == "Others") %>%
        dplyr::arrange(dplyr::desc(Degree))


      rr   <- sqrt(coords$x^2 + coords$y^2)
      ordc <- order(rr, coords$x, coords$y)
      coords <- coords[ordc, , drop = FALSE]


      placed_coords[["Others"]] <- dplyr::bind_cols(coords, nodes_oth)
    }
  }


  graph_ly_final <- dplyr::bind_rows(placed_coords[as.character(need_per_mod$modularity3)])

  if (isTRUE(jitter)) {
    graph_ly_final <- graph_ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }else{
    graph_ly_final <- graph_ly_final
  }

  ly_final <- graph_ly_final[, c("x","y")]


  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels


  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))


  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)


  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)


  list(
    layout = ly_final,
    graph_obj = graph_obj_sort,
    graph_ly_final = graph_ly_final,
    ggplot_data = ggplot_data
  )
}


module_layout4 <- function(graph_obj,
                           layout,                # data.frame(x, y)
                           center = TRUE,
                           shrink = 0.9,
                           k_nn = 8,
                           push_others_delta = 0.2,
                           jitter,
                           jitter_sd
                           # seed = seed
){
  # set.seed(seed)


  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  node_df %>%
    dplyr::count(modularity3, name = "size") %>%
    dplyr::arrange(desc(size)) %>%
    dplyr::mutate(modularity4 = factor(modularity3,
                                       levels = c(setdiff(modularity3, "Others"), "Others"),
                                       ordered = TRUE)) %>%
    dplyr::arrange(modularity4) %>%
    dplyr::mutate(modularity4 = as.character(modularity4)) %>%
    dplyr::pull(modularity4) -> mod_levels


  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))


  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)


  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)


  ly_final <- data.frame(x = layout$x,
                         y = layout$y)

  if (isTRUE(jitter)) {
    ly_final <- ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }else{
    ly_final <- ly_final
  }
  # combine

  graph_ly_final <- dplyr::bind_cols(
    ly_final,
    graph_obj_sort %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
  ) %>%
    dplyr::mutate(Modularity = droplevels(Modularity))

  # Per-module shrink toward each module's centroid.
  #
  # module_layout4 pairs with module-aware first-tier layouts
  # (circular_modules_gephi_layout, bipartite_gephi_layout, WGCNA, ...),
  # where the first tier has already placed nodes into geometrically
  # distinct per-module clusters.  The right "shrink" semantic here is
  # therefore to contract each cluster toward its own centroid -- module
  # positions stay put, module interiors become tighter, inter-module
  # whitespace grows.  (Contrast with module_layout / module_layout3,
  # which use shrink_rings_global() because their first-tier layouts are
  # not module-clustered in coordinate space.)
  if (!is.null(shrink) && shrink != 1) {
    graph_ly_final <- graph_ly_final %>%
      dplyr::group_by(modularity3) %>%
      dplyr::mutate(
        x = mean(x) + (x - mean(x)) * shrink,
        y = mean(y) + (y - mean(y)) * shrink
      ) %>%
      dplyr::ungroup()

    # Keep the standalone `ly_final` data.frame in sync with the shrunk
    # coordinates in `graph_ly_final` (same row order guaranteed by the
    # bind_cols above).  Both are returned and consumed downstream.
    ly_final$x <- graph_ly_final$x
    ly_final$y <- graph_ly_final$y
  }

  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)

  # get location result

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort,
              graph_ly_final = graph_ly_final,
              ggplot_data = ggplot_data
  ))

}


module_layout5 <- function(graph_obj,
                           layout,                # data.frame(x, y)
                           center = TRUE,
                           shrink = 0.9,
                           k_nn = 8,
                           push_others_delta = 0.2,
                           jitter,
                           jitter_sd
                           # seed = seed
){
  # set.seed(seed)


  node_df <- graph_obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble()


  # ---- Determine module order (matches the order used by create_layout_multirings) ----
  # Default: order modules by size ascending (smallest first), with "Others" at the end.
  # If the caller has explicitly customised module order via
  # `update_graph_modules*(levels = ...)`, the graph carries
  # `.modularity_user_ordered = TRUE`; in that case we respect the existing
  # `Modularity` factor levels so that the node order here matches the ring
  # order produced by `create_layout_multirings`.
  user_ordered <- isTRUE(igraph::graph_attr(graph_obj, ".modularity_user_ordered"))

  if (user_ordered && is.factor(node_df$Modularity)) {
    mod_levels <- levels(droplevels(node_df$Modularity))
  } else {
    node_df %>%
      dplyr::count(modularity3, name = "size") %>%
      dplyr::arrange(size) %>%
      dplyr::mutate(modularity4 = factor(modularity3,
                                         levels = c(setdiff(modularity3, "Others"), "Others"),
                                         ordered = TRUE)) %>%
      dplyr::arrange(modularity4) %>%
      dplyr::mutate(modularity4 = as.character(modularity4)) %>%
      dplyr::pull(modularity4) -> mod_levels
  }


  node_df_sorted <- node_df %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels)) %>%
    tidygraph::arrange(modularity3, dplyr::desc(Degree))


  node_df_sorted_number <- node_df_sorted %>%
    dplyr::count(modularity3)


  graph_obj_sort <- graph_obj %>%
    tidygraph::mutate(modularity3 = factor(modularity3, levels = mod_levels, ordered = TRUE)) %>%
    tidygraph::arrange(modularity3)


  ly_final <- data.frame(x = layout$x,
                         y = layout$y)

  if (isTRUE(jitter)) {
    ly_final <- ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }else{
    ly_final <- ly_final
  }

  # combine

  graph_ly_final <- dplyr::bind_cols(
    ly_final,
    graph_obj_sort %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
  ) %>%
    dplyr::mutate(Modularity = droplevels(Modularity))


  ggplot_data <- get_location(graph_ly_final, graph_obj_sort)

  # get location result

  return(list(layout = ly_final,
              graph_obj = graph_obj_sort,
              graph_ly_final = graph_ly_final,
              ggplot_data = ggplot_data
  ))

}


#' @noRd
module_layout_passthrough <- function(graph_obj,
                                      layout,                # data.frame(x, y)
                                      jitter = FALSE,
                                      jitter_sd = 0.01) {
  # A no-op module arranger. Self-arranging first-tier layouts (e.g.
  # create_layout_circlepack) already place BOTH modules and nodes, so the
  # module-arrangement stage must not run module_layout*/2/3/4/5 (which would
  # re-arrange the modules and destroy the packing). This wrapper simply keeps
  # the incoming coordinates and packages them into the standard ly1_1 list.
  #
  # Coordinates in `layout` MUST already be in graph_obj's node order (bind_cols
  # below attaches them to node rows by position).
  ly_final <- data.frame(x = layout$x, y = layout$y)

  if (isTRUE(jitter)) {
    ly_final <- ly_final %>%
      dplyr::mutate(x = x + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd),
                    y = y + stats::rnorm(dplyr::n(), mean = 0, sd = jitter_sd))
  }

  graph_ly_final <- dplyr::bind_cols(
    ly_final,
    graph_obj %>%
      tidygraph::activate(nodes) %>%
      tidygraph::as_tibble()
  ) %>%
    dplyr::mutate(Modularity = droplevels(Modularity))

  ggplot_data <- get_location(graph_ly_final, graph_obj)

  return(list(layout = ly_final,
              graph_obj = graph_obj,
              graph_ly_final = graph_ly_final,
              ggplot_data = ggplot_data
  ))
}
