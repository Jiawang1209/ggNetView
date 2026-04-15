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
#> The max module in network is 7 we use the 7  modules for next analysis
obj
#> # A tbl_graph: 22 nodes and 18 edges
#> #
#> # An undirected simple graph with 5 components
#> #
#> # Node Data: 22 × 7 (active)
#>    name      modularity modularity2 modularity3 Modularity Degree Strength
#>    <chr>     <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl>
#>  1 feature33 1          1           1           1               3    1.41 
#>  2 feature7  1          1           1           1               2    0.951
#>  3 feature31 1          1           1           1               2    1.00 
#>  4 feature8  1          1           1           1               1    0.531
#>  5 feature13 1          1           1           1               1    0.472
#>  6 feature28 1          1           1           1               1    0.458
#>  7 feature36 2          2           2           2               4    2.25 
#>  8 feature10 2          2           2           2               3    1.80 
#>  9 feature12 2          2           2           2               2    1.17 
#> 10 feature15 2          2           2           2               2    1.21 
#> # ℹ 12 more rows
#> #
#> # Edge Data: 18 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1    13    14  0.454       0.454 Positive      
#> 2     8    15  0.458       0.458 Positive      
#> 3    15    16  0.460       0.460 Positive      
#> # ℹ 15 more rows
# }
```
