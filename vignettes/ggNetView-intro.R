## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 6,
  fig.height = 5,
  dpi       = 72
)


## ----setup--------------------------------------------------------------------
library(ggNetView)
library(igraph)


## ----build--------------------------------------------------------------------
data(ppi_example)

ig <- igraph::graph_from_data_frame(
  d        = ppi_example$ppi,
  vertices = ppi_example$annotation,
  directed = FALSE
)

graph_obj <- build_graph_from_igraph(
  igraph        = ig,
  module.method = "Fast_greedy"
)

graph_obj


## ----inspect------------------------------------------------------------------
nodes <- get_graph_nodes(graph_obj)
head(nodes)

info <- get_info_from_graph(graph_obj)
lapply(info, head, 3)


## ----plot-fr------------------------------------------------------------------
ggNetView(
  graph_obj,
  layout    = "fr",
  seed      = 1,
  node_size_range = c(2, 8),
  node_fill   = "Modularity",
  module_label     = FALSE
)


## ----plot-circle--------------------------------------------------------------
ggNetView(
  graph_obj,
  layout    = "circle",
  seed      = 1,
  node_size_range = c(2, 8),
  node_fill   = "Modularity"
)


## ----sessioninfo--------------------------------------------------------------
sessionInfo()

