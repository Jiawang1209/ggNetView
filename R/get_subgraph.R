get_subgraph <- function(graph_obj){

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

  return(sub_graph)

}

