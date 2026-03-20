# Create sample-level subgraph topology from a global network

For each sample (each column in `mat`), OTUs with non-zero abundance are
treated as present (`OTU != 0`). The function extracts a subgraph from
`graph_obj` using
[`igraph::subgraph()`](https://r.igraph.org/reference/subgraph.html),
then computes topology for each sample-specific subgraph via
[`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

## Usage

``` r
get_sample_subgraph_topology(
  graph_obj,
  mat = NULL,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  r.threshold = 0.7,
  p.threshold = 0.05,
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
  cor.method = c("pearson", "kendall", "spearman"),
  proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  SpiecEasi.method = c("mb", "glasso"),
  sparcc_R = 20,
  bootstrap = 100
)
```

## Arguments

- graph_obj:

  An graph object from `build_graph_from_mat` or `build_graph_from_df`.

- mat:

  Numeric matrix used to build `graph_obj`. Rows are OTUs and columns
  are samples.

- transfrom.method:

  Character. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- r.threshold:

  Numeric. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- p.threshold:

  Numeric. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- method:

  Character. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- cor.method:

  Character. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- proc:

  Character. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

- sparcc_R:

  Integer. Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md)
  for SparCC p-values. Default 20.

- bootstrap:

  Numeric (default = 100). Passed to
  [`get_network_topology()`](https://jiawang1209.github.io/ggNetView/reference/get_network_topology.md).

## Value

A list with:

- `subgraph_list`: sample-wise `tbl_graph` subgraphs

- `topology`: merged sample-wise topology table

- `Robustness`: merged sample-wise robustness table

- `sample_stat`: sample-wise node/edge counts and status

## Examples

``` r
NULL
#> NULL
```
