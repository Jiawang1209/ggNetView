# Create network topology

Generates a data frame containing network topology information from
either a pre-built graph object or directly from an adjacency matrix.

## Usage

``` r
get_network_topology(
  graph_obj = NULL,
  graph_obj_list = NULL,
  mat = NULL,
  graph_mat_list = NULL,
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

  An graph object from build_graph_from_mat or build_graph_from_df. The
  network object to be visualized.

- graph_obj_list:

  A list of graph objects. Optional alternative to `graph_obj`. Each
  element is analyzed separately.

- mat:

  Numeric Matrix (default = NULL) The matrix to build graph_obj

- graph_mat_list:

  A list of matrices corresponding to `graph_obj_list`. Optional. Each
  element is paired with the graph object at the same position and used
  for topology analysis separately.

- transfrom.method:

  Character. Data transformation methods applied before correlation
  analysis. Options include: "none" (raw data), "scale" (z-score
  standardization), "center" (mean centering only), "log2" (log2
  transfrom), "log10" (log10 transfrom), "ln" (natural transfrom ),
  "rrarefy" (random rarefaction using
  [`vegan::rrarefy`](https://vegandevs.github.io/vegan/reference/rarefy.html)),
  "rrarefy_relative" (rarefy then convert to relative abundance).

- r.threshold:

  Numeric. Correlation coefficient threshold; edges are kept only if
  \|r\| \>= r.threshold.

- p.threshold:

  Numeric. Significance threshold for correlations; edges are kept only
  if p \< p.threshold.

- method:

  Character. Relationship analysis methods. Options include: "WGCNA",
  "SpiecEasi", "SPARCC" and "cor".

- cor.method:

  Character. Correlation analysis method. Options include "pearson",
  "kendall", and "spearman".

- proc:

  Character. Correlation p-value adjustment methods. Options include:
  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and
  "none".

- SpiecEasi.method:

  Character. Inverse-covariance estimation method passed to SpiecEasi
  when `method = "SpiecEasi"`. One of `"mb"` (Meinshausen-Buehlmann,
  default) or `"glasso"`.

- sparcc_R:

  Integer. Number of bootstrap/permutation replicates for SparCC
  p-values (when `method = "SPARCC"`). Default 20.

- bootstrap:

  Numeric (default = 100). Number of bootstrap iterations for stability
  analysis

## Value

A list containing topology output and robustness output for a single
network. When `graph_obj_list` is provided, returns a named list of such
results.

## Examples

``` r
# \donttest{
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
topo <- get_network_topology(graph_obj = obj, bootstrap = 10)
head(topo$topology)
#> # A tibble: 6 × 3
#>   Topology Target_network Random_nerwork
#>   <chr>             <dbl>          <dbl>
#> 1 Node           100            100     
#> 2 Edge            50             50     
#> 3 Degree           1              1     
#> 4 Distance        59.8            3.43  
#> 5 Diameter       109.             9.2   
#> 6 Density          0.0101         0.0101
# }
```
