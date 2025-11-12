#' Extract subgrah from graph object
#'
#' @param graph_obj An graph object from build_graph_from_mat or build_graph_from_df.
#' The network object to be visualized.
#' @param select_module a character vectors
#' Select the module name in graph object
#'
#' @returns list
#' @export
#'
#' @examples NULL
get_subgraph <- function(graph_obj, select_module = NULL){

  # get obj
  obj <- graph_obj

  # get module name
  module_name <- obj %>%
    tidygraph::activate(nodes) %>%
    tidygraph::as_tibble() %>%
    dplyr::pull(Modularity) %>%
    levels()

  # get module list
  module_list <- purrr::map(obj %>%
                              tidygraph::activate(nodes) %>%
                              tidygraph::as_tibble() %>%
                              dplyr::group_split(Modularity),
                            ~.x)
  names(module_list) <- module_name

  # get module ID
  id_list <- purrr::map(obj %>%
                          tidygraph::activate(nodes) %>%
                          tidygraph::as_tibble() %>%
                          dplyr::group_split(Modularity),
                        ~.x[[1]])
  names(id_list) <- module_name

  # create sub_graph object
  sub_graph <- list()

  for (i in module_name) {

    sub_graph[[i]] <- igraph::subgraph(tidygraph::as.igraph(obj),
                                       id_list[[i]]) %>%
      tidygraph::as_tbl_graph()

  }

  # stat
  stat_module <- purrr::map(id_list, ~length(.x)) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Module") %>%
    dplyr::rename(Number = V1)

  print(stat_module)


  if (!is.null(select_module)) {
    graph_select <- obj %>%
      tidygraph::filter(as.character(Modularity) %in% select_module)
  }else{
    graph_select <- NULL
  }


  return(list(sub_graph_all = sub_graph,
              stat_module = stat_module,
              sub_graph_select = graph_select))

}

