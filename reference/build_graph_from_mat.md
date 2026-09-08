# Build a correlation-based network from a matrix

Build a correlation-based network from a matrix

## Usage

``` r
build_graph_from_mat(
  mat,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  r.threshold = 0.7,
  p.threshold = 0.05,
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor", "Hmisc"),
  cor.method = c("pearson", "kendall", "spearman"),
  proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  SpiecEasi.method = c("mb", "glasso"),
  sparcc_R = 20,
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- mat:

  Numeric matrix. A numeric matrix with samples in colums and variables
  in rows

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
  "SpiecEasi", "SPARCC", "cor", and "Hmisc".

- cor.method:

  Character. Correlation analysis method. Options include "pearson",
  "kendall", and "spearman".

- proc:

  Character. Correlation p-value adjustment methods. Options include:
  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and
  "none".

- module.method:

  Character. Network community detection (module identification) method.
  Options include "Fast_greedy", "Walktrap", "Edge_betweenness", and
  "Spinglass".

- SpiecEasi.method:

  Character. Method used in `SpiecEasi` network inference; options
  include "mb" and "glasso".

- sparcc_R:

  Integer. Number of bootstrap/permutation replicates for SparCC
  p-values (when `method = "SPARCC"`). Default 20.

- node_annotation:

  Data frame. Optional node annotation table, containing metadata such
  as taxonomy or functional categories.

- top_modules:

  Integer. Number of top-ranked modules to retain for downstream
  visualization or analysis.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

A graph object representing the correlation-based microbial network.
Node/edge attributes include correlation statistics and (optionally)
module labels.

## Examples

``` r
# \donttest{
set.seed(1)
mat <- matrix(stats::rnorm(40 * 20), nrow = 40, ncol = 20)
rownames(mat) <- paste0("feature", seq_len(40))
colnames(mat) <- paste0("sample",  seq_len(20))
obj <- build_graph_from_mat(
  mat           = mat,
  method        = "cor",
  cor.method    = "pearson",
  proc          = "none",
  r.threshold   = 0.3,
  p.threshold   = 0.05,
  module.method = "Fast_greedy"
)
#> The max module in network is 6 we use the 6  modules for next analysis
obj
#> # A tbl_graph: 38 nodes and 46 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 38 × 7 (active)
#>    name      modularity modularity2 modularity3 Modularity Degree Strength
#>    <chr>     <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl>
#>  1 feature36 2          2           2           2               6    3.23 
#>  2 feature10 2          2           2           2               5    3.08 
#>  3 feature15 2          2           2           2               5    2.79 
#>  4 feature20 2          2           2           2               4    2.34 
#>  5 feature22 2          2           2           2               3    1.63 
#>  6 feature26 2          2           2           2               3    1.60 
#>  7 feature12 2          2           2           2               2    1.17 
#>  8 feature37 2          2           2           2               2    1.05 
#>  9 feature2  2          2           2           2               1    0.449
#> 10 feature24 2          2           2           2               1    0.495
#> # ℹ 28 more rows
#> #
#> # Edge Data: 46 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1    13    14  0.454       0.454 Positive      
#> 2    11    13  0.459      -0.459 Negative      
#> 3     3     9  0.449      -0.449 Negative      
#> # ℹ 43 more rows
# }
```
