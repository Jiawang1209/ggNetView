get_location <- function(graph_ly_final, graph_obj_sort){

  ggplot_node_df <- graph_ly_final %>%
    dplyr::select(name, x, y, dplyr::everything()) %>%
    dplyr::mutate(id = 1:nrow(.))

  ggplot_edge_df <- graph_obj_sort %>%
    tidygraph::activate(edges) %>%
    tidygraph::as_tibble() %>%
    dplyr::left_join(ggplot_node_df %>% dplyr::select(name, id, x, y),  by = c("from" = "id")) %>%
    dplyr::rename(from_id = name,
                  from_x = x,
                  from_y = y) %>%
    dplyr::left_join(ggplot_node_df %>% dplyr::select(name, id, x, y),  by = c("to" = "id")) %>%
    dplyr::rename(to_id = name,
                  to_x = x,
                  to_y = y)

  return(list(ggplot_node_df = ggplot_node_df,
              ggplot_edge_df = ggplot_edge_df))
}

