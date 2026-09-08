# Compute per-node centralities and attach them to a graph object

Computes a panel of standard per-node centrality measures and adds them
as new vertex columns on the input \`tbl_graph\`. Use this for
node-importance analysis when you want to rank nodes (or map a
centrality to a visual aesthetic in \[ggNetView()\]) by something more
informative than degree alone.

## Usage

``` r
get_node_centrality(
  graph_obj,
  measures = c("Betweenness", "Closeness", "Eigenvector", "PageRank", "Hub_score",
    "Authority_score", "Coreness", "Harmonic"),
  weighted = FALSE,
  overwrite = TRUE
)
```

## Arguments

- graph_obj:

  A \`tbl_graph\` produced by any \`build_graph_from\_\*()\`
  constructor.

- measures:

  Character vector. Which centralities to compute. Pass \`"all"\` to
  compute every supported measure. Defaults to all eight measures listed
  in **Available measures**.

- weighted:

  Logical (default \`FALSE\`). If \`TRUE\`, the edge \`weight\`
  attribute is used as a distance (\`igraph\` convention: higher weight
  = farther). Because correlation networks have \`weight =
  \|correlation\|\` – where higher means \*closer\* – the distance used
  internally is \`1 / weight\` so that strongly correlated pairs count
  as short paths. Set \`FALSE\` (default) for the textbook unweighted
  versions.

- overwrite:

  Logical (default \`TRUE\`). If a measure column already exists on the
  input graph, controls whether to overwrite it (silent overwrite when
  \`TRUE\`; warning + skip when \`FALSE\`).

## Value

A \`tbl_graph\` whose node table is augmented with one column per
requested measure (using the column names listed in **Available
measures**). Other vertex / edge columns are preserved verbatim.

## Details

Note that \[get_network_topology()\] reports the same families of
metrics but as \*\*network-level summaries\*\* (a single mean / sum per
network). \`get_node_centrality()\` is the per-node counterpart – it
keeps every individual value so you can sort, rank, threshold, or colour
by it.

## Available measures

All measures wrap the corresponding \`igraph\` function:

- \`"Betweenness"\`:

  Number of shortest paths through each node (\`igraph::betweenness\`).

- \`"Closeness"\`:

  Inverse mean shortest-path distance from a node to every other node
  (\`igraph::closeness\`, \`mode = "all"\`). Returns \`NaN\` for nodes
  in their own connected component when the component is a singleton.

- \`"Eigenvector"\`:

  Eigenvector centrality (\`igraph::eigen_centrality\`).

- \`"PageRank"\`:

  Google PageRank score (\`igraph::page_rank\`).

- \`"Hub_score"\`:

  HITS hub score (\`igraph::hub_score\`).

- \`"Authority_score"\`:

  HITS authority score (\`igraph::authority_score\`).

- \`"Coreness"\`:

  k-core membership of each vertex (\`igraph::coreness\`).

- \`"Harmonic"\`:

  Harmonic centrality (\`igraph::harmonic_centrality\`); robust to
  disconnected graphs where Closeness becomes ill-defined.

## See also

\[get_network_topology()\] for network-level summaries of the same
metrics; \[get_node_ivi()\] for an integrative importance score that
combines local, semi-local, and global centralities.

## Examples

``` r
# \donttest{
set.seed(1)
mat <- matrix(stats::rnorm(40 * 20), nrow = 40, ncol = 20)
rownames(mat) <- paste0("feature", seq_len(40))
colnames(mat) <- paste0("sample",  seq_len(20))
obj <- build_graph_from_mat(
  mat = mat, method = "cor", cor.method = "pearson",
  proc = "none", r.threshold = 0.3, p.threshold = 0.05
)
#> The max module in network is 6 we use the 6  modules for next analysis

obj_aug <- get_node_centrality(obj)
obj_aug %>%
  tidygraph::activate(nodes) %>%
  tidygraph::as_tibble() %>%
  dplyr::arrange(dplyr::desc(Betweenness)) %>%
  utils::head(5)
#> # A tibble: 5 × 15
#>   name      modularity modularity2 modularity3 Modularity Degree Strength
#>   <chr>     <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl>
#> 1 feature40 5          5           5           5               4     2.06
#> 2 feature7  3          3           3           3               4     2.20
#> 3 feature36 2          2           2           2               6     3.23
#> 4 feature23 3          3           3           3               2     1.12
#> 5 feature39 4          4           4           4               3     1.56
#> # ℹ 8 more variables: Betweenness <dbl>, Closeness <dbl>, Eigenvector <dbl>,
#> #   PageRank <dbl>, Hub_score <dbl>, Authority_score <dbl>, Coreness <dbl>,
#> #   Harmonic <dbl>
# }
```
