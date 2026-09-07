#' Create network topology
#'
#' Generates a data frame containing network topology information from either a
#' pre-built graph object or directly from an adjacency matrix.
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#'   The network object to be visualized.
#' @param graph_obj_list A list of graph objects. Optional alternative to
#'   \code{graph_obj}. Each element is analyzed separately.
#' @param mat Numeric Matrix (default = NULL)
#' The matrix to build graph_obj
#' @param graph_mat_list A list of matrices corresponding to
#'   \code{graph_obj_list}. Optional. Each element is paired with the graph
#'   object at the same position and used for topology analysis separately.
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
#' @param SpiecEasi.method Character. Inverse-covariance estimation method
#'   passed to SpiecEasi when \code{method = "SpiecEasi"}.
#'   One of \code{"mb"} (Meinshausen-Buehlmann, default) or \code{"glasso"}.
#' @param sparcc_R Integer.
#' Number of bootstrap/permutation replicates for SparCC p-values (when \code{method = "SPARCC"}).
#' Default 20.
#' @param bootstrap Numeric  (default = 100).
#' Number of bootstrap iterations for stability analysis
#' @param seed Integer (default = 1115).
#' Random seed for reproducibility. The correlation/network estimation
#' (e.g. SparCC, SpiecEasi, rarefaction) and robustness resampling steps use
#' random number generation; setting this makes results deterministic and
#' consistent with \code{get_network_topology_parallel()}.
#'
#' @returns A list containing topology output and robustness output for a single
#'   network. When \code{graph_obj_list} is provided, returns a named list of
#'   such results.
#' @export
#'
#' @examples
#' \donttest{
#' data(ppi_example)
#' obj <- build_graph_from_df(
#'   df              = ppi_example$ppi,
#'   node_annotation = ppi_example$annotation
#' )
#' topo <- get_network_topology(graph_obj = obj, bootstrap = 10)
#' head(topo$topology)
#' }
get_network_topology <- function(graph_obj = NULL,
                                 graph_obj_list = NULL,
                                 mat = NULL,
                                 graph_mat_list = NULL,
                                 transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy", "rrarefy_relative"),
                                 r.threshold = 0.7,
                                 p.threshold = 0.05,
                                 method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
                                 cor.method = c("pearson", "kendall", "spearman"),
                                 proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                 SpiecEasi.method = c("mb", "glasso"),
                                 sparcc_R = 20,
                                 bootstrap = 100,
                                 seed = 1115){

  if (!is.null(graph_obj_list) || !is.null(graph_mat_list)) {
    if (!is.null(graph_obj)) {
      stop("Please provide either `graph_obj` or `graph_obj_list`, not both.", call. = FALSE)
    }
    if (is.null(graph_obj_list)) {
      stop("`graph_obj_list` must be provided when using list input mode.", call. = FALSE)
    }
    if (!is.list(graph_obj_list)) {
      stop("`graph_obj_list` must be a list.", call. = FALSE)
    }
    if (length(graph_obj_list) == 0) {
      stop("`graph_obj_list` must contain at least one graph object.", call. = FALSE)
    }
    if (!is.null(graph_mat_list)) {
      if (!is.list(graph_mat_list)) {
        stop("`graph_mat_list` must be a list.", call. = FALSE)
      }
      if (length(graph_mat_list) != length(graph_obj_list)) {
        stop("`graph_obj_list` and `graph_mat_list` must have the same length.", call. = FALSE)
      }
    }

    result_names <- names(graph_obj_list)
    if (is.null(result_names) && !is.null(graph_mat_list)) {
      result_names <- names(graph_mat_list)
    }
    if (is.null(result_names)) {
      result_names <- paste0("network_", seq_along(graph_obj_list))
    }

    out <- lapply(seq_along(graph_obj_list), function(i) {
      get_network_topology(
        graph_obj = graph_obj_list[[i]],
        mat = if (is.null(graph_mat_list)) NULL else graph_mat_list[[i]],
        transfrom.method = transfrom.method,
        r.threshold = r.threshold,
        p.threshold = p.threshold,
        method = method,
        cor.method = cor.method,
        proc = proc,
        SpiecEasi.method = SpiecEasi.method,
        sparcc_R = sparcc_R,
        bootstrap = bootstrap,
        seed = seed + i - 1L
      )
    })
    names(out) <- result_names
    return(out)
  }

  if (is.null(graph_obj)) {
    stop("`graph_obj` must be provided unless `graph_obj_list` is used.", call. = FALSE)
  }

  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
    stop("`seed` must be a single finite numeric.", call. = FALSE)
  }
  .ggnv_local_seed(seed)

  # create igraph object
  ig <- tidygraph::as.igraph(graph_obj)

  # argument check
  method <- match.arg(method)
  transfrom.method <-  match.arg(transfrom.method)
  cor.method <- match.arg(cor.method)
  proc <- match.arg(proc)
  SpiecEasi.method <- match.arg(SpiecEasi.method)
  sparcc_R <- as.integer(sparcc_R)[1L]
  if (is.na(sparcc_R) || sparcc_R < 1L) {
    stop("`sparcc_R` must be a positive integer.", call. = FALSE)
  }

  # Correct over the n(n-1)/2 unique off-diagonal tests, not all n^2 cells:
  # the duplicated triangle and the uninformative diagonal otherwise distort
  # the ranks that Benjamini-Hochberg depends on. See .ggnv_adjust_p_matrix().
  adjust_p_matrix <- function(p_mat, proc_method) {
    .ggnv_adjust_p_matrix(p_mat, proc_method)
  }

  if (is.null(mat)) {
    graph_obj = graph_obj

    network.raw = NULL

    cohesion_out <- data.frame(
      cohension_position = NA,
      cohension_negative = NA
    )

    Weighted.simu <- data.frame(
      remain.mean = NA,
      remain.sd = NA,
      remain.se = NA
    )

    Unweighted.simu <- data.frame(
      remain.mean = NA,
      remain.sd = NA,
      remain.se = NA
    )



  }else{
    # data transfrom
    mat <- apply_transform_method(mat, transfrom.method)

    # calculate correlation

    # WGCNA
    if (method == "WGCNA") {
      sp.ra <- colMeans(t(mat))

      # WGCNA for correlation
      occor <- WGCNA::corAndPvalue(t(mat), method = cor.method)
      occor.p <- adjust_p_matrix(occor$p, proc)

      # R and pvalue
      occor.r <- occor$cor
      diag(occor.r) <- 0
      occor.r[is.na(occor.p) | is.na(occor.r) | occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
      occor.r[is.na(occor.r)]=0

      network.raw <- occor.r[colSums(abs(occor.r)) > 0, colSums(abs(occor.r)) > 0, drop = FALSE]
      sp.ra2 <- sp.ra[colSums(abs(occor.r))>0]

    }

    # SpiecEasi selects edges via StARS model selection (sparse inverse
    # covariance), NOT a correlation threshold. Use the weighted estimate
    # directly -- partial correlation (glasso) / symmetrised beta (mb) -- whose
    # non-zero support IS the selected edge set. `r.threshold` does not apply.
    if (method == "SpiecEasi") {
      sp.ra <- colMeans(t(mat))
      se_output <- if (SpiecEasi.method == "mb") "beta" else "partial_correlation"
      am <- spieceasi_matrix_rcpp(as.matrix(t(mat)), method = SpiecEasi.method, output = se_output,
                                  lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 50))
      rownames(am) <- rownames(mat)
      colnames(am) <- rownames(mat)
      am2 <- am
      am2[is.na(am2)] <- 0
      diag(am2) <- 0
      network.raw <- am2[colSums(abs(am2)) > 0, colSums(abs(am2)) > 0, drop = FALSE]
      sp.ra2 <- sp.ra[colSums(abs(am2)) > 0]
    }

    # SparCC (sparcc_matrix_rcpp: filter by r.threshold and p.threshold)
    if (method == "SPARCC") {
      sp.ra <- colMeans(t(mat))
      occor.r <- sparcc_matrix_rcpp(as.matrix(t(mat)))
      p_mat <- sparcc_pvalue_rcpp(as.matrix(t(mat)), R = sparcc_R, seed = seed)
      diag(occor.r) <- 0
      occor.r[abs(occor.r) < r.threshold | is.na(p_mat) | p_mat > p.threshold] <- 0
      occor.r[is.na(occor.r)] <- 0
      rownames(occor.r) <- rownames(mat)
      colnames(occor.r) <- rownames(mat)
      SparCC_graph2 <- Matrix::Matrix(occor.r, sparse = TRUE)
      # coerce back to a dense base matrix: downstream .cohension_compute()
      # calls as.data.frame(), which cannot coerce a sparse dgCMatrix/dsCMatrix.
      network.raw <- as.matrix(SparCC_graph2[colSums(abs(SparCC_graph2)) > 0, colSums(abs(SparCC_graph2)) > 0])
      sp.ra2 <- sp.ra[colSums(abs(SparCC_graph2)) > 0]
    }

    # cor
    if (method == "cor") {
      sp.ra <- colMeans(t(mat))
      # Cor for correlation
      occor <- psych::corr.test(t(mat), method = cor.method, adjust = "none")
      occor.p <- adjust_p_matrix(occor$p, proc)

      # R and pvalue
      occor.r <- occor$r
      diag(occor.r) <- 0
      occor.r[is.na(occor.p) | is.na(occor.r) | occor.p > p.threshold | abs(occor.r) < r.threshold] = 0
      occor.r[is.na(occor.r)]=0

      network.raw <- occor.r[colSums(abs(occor.r)) > 0, colSums(abs(occor.r)) > 0, drop = FALSE]
      sp.ra2 <- sp.ra[colSums(abs(occor.r))>0]

    }
  }



  .rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=TRUE){
    id.rm<-sample(seq_len(nrow(netRaw)), round(nrow(netRaw)*rm.percent))
    net.Raw=netRaw
    net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;
    if (abundance.weighted){
      net.stength= net.Raw*sp.ra
    } else {
      net.stength= net.Raw
    }

    sp.meanInteration<-colMeans(net.stength)

    id.rm2<- which(sp.meanInteration<=0)
    remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)

    remain.percent
  }


  .rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=TRUE,bootstrap=bootstrap){
    t(sapply(rm.p.list,function(x){
      remains=sapply(1:bootstrap,function(i){
        .rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
      })
      remain.mean=mean(remains)
      remain.sd=sd(remains)
      remain.se=sd(remains)/(bootstrap^0.5)
      result<-c(remain.mean,remain.sd,remain.se)
      names(result)<-c("remain.mean","remain.sd","remain.se")
      result
    }))
  }


  # Delegates to the shared .compute_cohension() so the cohesion definition and
  # its NA guards stay identical to get_network_topology_parallel(). `mat` is the
  # abundance matrix captured from this function's scope.
  .cohension_compute <- function(network.raw){
    .compute_cohension(network.raw = network.raw, mat = mat)
  }


  # function
  .network.efficiency <- function(ig){
    dd <- 1/igraph::distances(ig)
    diag(dd) <- NA
    efficiency <- mean(dd, na.rm=TRUE)
    return(efficiency)
  }

  .info.centrality.vertex <- function(ig, net=NULL, verbose=FALSE){
    if(is.null(net)) {
      net <- .network.efficiency(ig)
      }
    if(is.numeric(net)==F){
      net <- .network.efficiency(ig)
    }
    count <- c()
    for(i in seq_along(igraph::V(ig))){
      count <- c(count, (net - .network.efficiency(igraph::delete_vertices(ig, i)))/net)
      if(verbose){
        print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
      }
    }
    return(count)
  }

  # get node and edge number
  .get_topology <- function(ig){

    # compute network topology
    # node number
    node_number <- igraph::vcount(ig)

    # edge number
    edge_number <- igraph::ecount(ig)

    # degree
    degree_average <- mean(igraph::degree(ig, mode = "all"))
    degree_average_out <- mean(igraph::degree(ig, mode = "out"))
    degree_average_in <- mean(igraph::degree(ig, mode = "in"))
    degree_average_total <- mean(igraph::degree(ig, mode = "total"))
    # same output

    # mean_distance
    distance_mean <- igraph::mean_distance(ig, directed = FALSE)

    # diameter
    dia <- igraph::diameter(ig, directed = FALSE)

    # transitivity (global/local)
    clust_global <- igraph::transitivity(ig, type = "global")
    clust_local  <- mean(igraph::transitivity(ig, type = "local"), na.rm = TRUE)

    # density
    density <- igraph::edge_density(ig)

    # shortest path between vertices
    network_efficiency <- .network.efficiency(ig)

    # Information Centrality
    network_info.centrality <- sum(.info.centrality.vertex(ig))

    # betweenness
    betweenness_vals <- mean(igraph::betweenness(ig), na.rm = TRUE)
    betweenness_edge_vals <- mean(igraph::edge_betweenness(ig), na.rm = TRUE)

    # closeness
    closeness_vals <- mean(igraph::closeness(ig), na.rm = TRUE)
    closeness_vals_out <- mean(igraph::closeness(ig, mode = "out"), na.rm = TRUE)
    closeness_vals_in <- mean(igraph::closeness(ig, mode = "in"), na.rm = TRUE)
    closeness_vals_all <- mean(igraph::closeness(ig, mode = "all"), na.rm = TRUE)
    closeness_vals_total <- mean(igraph::closeness(ig, mode = "total"), na.rm = TRUE)
    # same output


    # eigen_centrality
    eigen_vals <- mean(igraph::eigen_centrality(ig)$vector, na.rm = TRUE)

    # modularity
    fc <- igraph::cluster_fast_greedy(ig)
    modularity <- igraph::modularity(ig, igraph::membership(fc))

    # K-core decomposition of graphs
    kcore_vals  <- igraph::coreness(ig)
    kcore_mean  <- mean(kcore_vals, na.rm = TRUE)
    kcore_max   <- max(kcore_vals, na.rm = TRUE)
    kcore_min   <- min(kcore_vals, na.rm = TRUE)

    # Cohension from igraph
    cohesion <- igraph::cohesion(ig)
    vertex_cohesion <- igraph::vertex_connectivity(ig)
    edge_cohesion   <- igraph::edge_connectivity(ig)

    # Cohension
    ## Positive Cohension
    position_cohension <- cohesion_out$cohension_position

    ## Negative Cohension
    negative_cohension <- cohesion_out$cohension_negative

    # Robustness
    ## Robustness_weight
    robustness_weight = mean(as.data.frame(Weighted.simu)[["remain.mean"]])

    ## Robustness_unweight
    robustness_unweight = mean(as.data.frame(Unweighted.simu)[["remain.mean"]])

    # Vulenrability
    vulenrability = max(.info.centrality.vertex(ig))


    # Stability

    out <- data.frame(
      Node = node_number,
      Edge = edge_number,
      Degree = degree_average,
      Distance = distance_mean,
      Diameter = dia,
      Density = density,
      Transitivity_global = clust_global,
      Transitivity_local = clust_local,
      Betweenness = betweenness_vals,
      Betweenness_edge = betweenness_edge_vals,
      Closeness = closeness_vals,
      Eigen_centrality = eigen_vals,
      Modularity = modularity,
      K_core_mean = kcore_mean,
      K_core_max = kcore_max,
      K_core_min = kcore_min,
      Network_efficiency = network_efficiency,
      Network_info.centrality = network_info.centrality,
      Cohension_Positive = position_cohension,
      Cohension_Negative = negative_cohension,
      Robustness_weight = robustness_weight,
      Robustness_unweight = robustness_unweight,
      Vulenrability = vulenrability,
      Stability = mean(c(robustness_weight, robustness_unweight))
      )

    return(out)


  }

  # compute
  if (is.null(mat)) {
    Weighted.simu <- data.frame(
      remain.mean = NA,
      remain.sd = NA,
      remain.se = NA
    )

    Unweighted.simu <- data.frame(
      remain.mean = NA,
      remain.sd = NA,
      remain.se = NA
    )

  }else{
    Weighted.simu <- .rmsimu(netRaw=network.raw,
                             rm.p.list=seq(0.05,1,by=0.05),
                             sp.ra=sp.ra2,
                             abundance.weighted=TRUE,
                             bootstrap=bootstrap)

    Unweighted.simu <- .rmsimu(netRaw=network.raw,
                               rm.p.list=seq(0.05,1,by=0.05),
                               sp.ra=sp.ra2,
                               abundance.weighted=FALSE,
                               bootstrap=bootstrap)
  }

  robustness <- data.frame(Proportion.removed = rep(seq(0.05,1,by=0.05),2),
                           rbind(Weighted.simu, Unweighted.simu),
                           weighted=rep(c("weighted", "unweighted"), each=length(seq(0.05,1,by=0.05)))
  )


  if (is.null(network.raw)) {

    cohesion_out <- data.frame(
      cohension_position = NA,
      cohension_negative = NA
    )

  }else{
    cohesion_out <- .cohension_compute(network.raw = network.raw)
  }



  # node topology
  network_topology <- .get_topology(ig = ig) %>%
    tibble::rownames_to_column(var = "ID") %>%
    tidyr::pivot_longer(cols = -ID, names_to = "Kind", values_to = "Value") %>%
    dplyr::select(-ID) %>%
    purrr::set_names(c("Topology", "Target_network"))

  network_topology

  # random topology
  #
  # The null ensemble is drawn through future.apply with `future.seed = TRUE`
  # so that the L'Ecuyer-CMRG streams -- and hence the numbers reported -- are
  # identical to get_network_topology_parallel(), whether that function runs
  # sequentially or on any number of workers. A plain loop here would consume
  # the global Mersenne-Twister stream instead and silently produce a different
  # ensemble from the same `seed`.
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

  random_topology <- future.apply::future_lapply(seq_len(bootstrap), function(i) {

    random_graph <- igraph::sample_gnm(n = igraph::vcount(ig),
                                       m = igraph::ecount(ig),
                                       directed = FALSE,
                                       loops = FALSE)

    .get_topology(ig = random_graph) %>%
      dplyr::mutate(Robustness_weight = NA,
                    Robustness_unweight = NA,
                    Cohension_Positive = NA,
                    Cohension_Negative = NA,
                    Stability = NA
                    )

  }, future.seed = TRUE)

  random_topology_df <- do.call(rbind, random_topology)

  # output
  # when bootstrap == 0 the loop above never ran, so random_topology_df is NULL
  # and colMeans() crashes. Emit a single all-NA row keyed to the target network
  # topology so the downstream left_join still produces the expected schema.
  if (is.null(random_topology_df) || nrow(random_topology_df) == 0L) {
    random_topology_df_mean <- data.frame(
      Topology       = network_topology$Topology,
      Random_nerwork = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    random_topology_df_mean <- colMeans(random_topology_df) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Kind") %>%
      purrr::set_names(c("Topology", "Random_nerwork"))
  }



  # output topology
  topology_df <- dplyr::left_join(network_topology,
                                  random_topology_df_mean,
                                  by = "Topology"
                                  )

  out <-  list(topology = topology_df,
               Robustness = robustness
 )

  return(out)

}
