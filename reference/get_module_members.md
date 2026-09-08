# Get module-to-OTU membership from a tbl_graph

Internal helper used by the modularity-based heatmap function. Returns a
named list mapping each module label to its set of OTU/node names.

## Usage

``` r
get_module_members(graph_obj, module_col = "Modularity", exclude_others = TRUE)
```

## Arguments

- graph_obj:

  A `tbl_graph` or `igraph` with a node `name` attribute and a module
  column.

- module_col:

  Module column name. If missing, falls back to one of `"Modularity"`,
  `"modularity3"`, `"modularity2"`.

- exclude_others:

  If `TRUE`, drop nodes labelled `"Others"`.

## Value

Named list. Names are module labels, values are character vectors of
OTU/node names.
