# Visualize network with custom layouts in different samples

Visualize network with custom layouts in different samples

## Usage

``` r
ggNetView_multi(
  mat,
  group_info,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  r.threshold = 0.7,
  p.threshold = 0.05,
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
  cor.method = c("pearson", "kendall", "spearman"),
  proc = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  SpiecEasi.method = c("mb", "glasso"),
  sparcc_R = 20,
  node_annotation = NULL,
  top_modules = 15,
  layout = NULL,
  ...,
  layout_nrow = NULL,
  layout_ncol = NULL,
  seed = 1115,
  nrow = NULL,
  ncol = NULL
)
```

## Arguments

- mat:

  Numeric matrix. A numeric matrix with samples in rows and variables in
  columns.

- group_info:

  DataFrame The group information contains: Sample and Group

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

  p.threshold Significance threshold for correlations; edges are kept
  only if p \< p.threshold.

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

- layout:

  Character string naming the layout passed to
  [`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  (e.g. "gephi", "fr", "circle", "square").

- ...:

  Additional arguments passed to
  [`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md)
  (node\_\*, edge\_\*, module_label\_\*, module_outline\_\*,
  network_outline\_\*, layout geometry, ...). Deprecated pre-0.2.0 names
  (e.g. `fill.by`, `pointsize`) are still accepted with a lifecycle
  warning.

- layout_nrow:

  Integer (default = NULL). Number of layout rows passed to `ggNetView`
  when using consensus-module grid layouts.

- layout_ncol:

  Integer (default = NULL). Number of layout columns passed to
  `ggNetView` when using consensus-module grid layouts.

- seed:

  Integer (default = 1115). Random seed for reproducibility.

- nrow:

  Integer (default = NULL). Number of rows in the combined patchwork
  plot.

- ncol:

  Integer (default = NULL). Number of columns in the combined patchwork
  plot.

## Value

A ggplot object representing the network visualization.

## Examples

``` r
if (FALSE) { # \dontrun{
# `mat` is a numeric matrix (features x samples) and
# `group_info` is a data frame with columns Sample and Group.
p <- ggNetView_multi(
  mat        = mat,
  group_info = group_info,
  method     = "cor",
  layout     = "fr"
)
} # }
```
