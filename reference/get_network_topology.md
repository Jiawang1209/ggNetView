# Create network topology

Generates a data frame containing network topology information from
either a pre-built graph object or directly from an adjacency matrix.

## Usage

``` r
get_network_topology(
  graph_obj,
  mat = NULL,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  r.threshold = 0.7,
  p.threshold = 0.05,
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
  cor.method = c("pearson", "kendall", "spearman"),
  proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH",
    "TSBH"),
  SpiecEasi.method = c("mb", "glasso"),
  sparcc_R = 20,
  bootstrap = 100
)
```

## Arguments

- graph_obj:

  An graph object from build_graph_from_mat or build_graph_from_df. The
  network object to be visualized.

- mat:

  Numeric Matrix (default = NULL) The matrix to build graph_obj

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
  "Bonferroni", "Holm", "Hochberg", " SidakSS", "SidakSD","BH", "BY",
  "ABH", and "TSBH".

- sparcc_R:

  Integer. Number of bootstrap/permutation replicates for SparCC
  p-values (when `method = "SPARCC"`). Default 20.

- bootstrap:

  Numeric (default = 100). Number of bootstrap iterations for stability
  analysis

## Value

data frame of network topolog

## Examples

``` r
NULL
#> NULL
```
