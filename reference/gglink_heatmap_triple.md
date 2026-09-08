# Visualize multi-orientation environmental-species correlation heatmaps2

Visualize multi-orientation environmental-species correlation heatmaps2

## Usage

``` r
gglink_heatmap_triple(
  Environment,
  Experiment,
  edge,
  node,
  sample_col = "Sample",
  delim = ",",
  hub_n = NULL,
  r = 6,
  cor.method = c("pearson", "kendall", "spearman"),
  cor.use = c("pairwise", "everything", "all", "complete", "na"),
  env_p_adjust = "none",
  link_p_adjust = "none",
  sig_breaks = c(0.05, 0.01, 0.001)
)
```

## Arguments

- Environment:

  character or data.frame File path or data frame of environment data.

- Experiment:

  character or data.frame File path or data frame of experiment data.

- edge:

  character or data.frame File path or data frame of edge data. Must
  contain columns `from` and `to`; an optional numeric `weight` column
  is mapped to edge colour/width (defaults to `1` when absent).

- node:

  character or data.frame File path or data frame of node data. Must
  contain a `node` column listing every node referenced by `edge`; node
  names matching columns of `Experiment` become the hub nodes anchored
  on the central heatmap. An optional `annotation` column drives node
  fill/shape (when absent it is derived automatically: `"Experiment"`
  for hub nodes, `"Environment"` otherwise).

- sample_col:

  Character (default = "Sample") Column name used as sample ID when
  input is a data frame or file.

- delim:

  Character (default = ",") Delimiter for reading input files.

- hub_n:

  Integer (default = NULL) If `NULL` (recommended), hubs are the
  `Experiment` variables present in `node`. If an integer, the `hub_n`
  highest out-degree nodes are used instead (they must then correspond
  one-to-one to the Experiment variables, and `node` rows must list
  circle nodes first).

- r:

  numeric (default = 6) Radius of the outer node circle.

- cor.method:

  Character (default = "pearson") Correlation method passed to
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html),
  used for both the Environment x Environment heatmap and the
  Environment x Experiment links. One of `"pearson"`, `"kendall"`,
  `"spearman"`.

- cor.use:

  Character (default = "pairwise") Missing-value handling passed to
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html);
  same vocabulary as
  [`gglink_heatmaps()`](https://jiawang1209.github.io/ggNetView/reference/gglink_heatmaps.md).
  Note the default differs from that function (`"everything"`) because
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html)
  itself defaults to `"pairwise"`, which is what this plot has always
  used.

- env_p_adjust:

  Character (default = "none") Multiple-testing correction for the
  Environment x Environment correlations (the significance stars on the
  triangular heatmap). Any method accepted by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), or
  `"none"`.

- link_p_adjust:

  Character (default = "none") Multiple-testing correction for the
  Environment x Experiment correlations (the linetype of the link
  segments). Any method accepted by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), or
  `"none"`. Note that
  [`psych::corr.test()`](https://rdrr.io/pkg/psych/man/corr.test.html)
  defaults to `"holm"` here but still reports raw p-values in `$p`, so
  the previous hard-coded call was in effect uncorrected; `"none"` keeps
  that behaviour.

- sig_breaks:

  Numeric vector of length 3 (default = c(0.05, 0.01, 0.001)) Strictly
  decreasing p-value cut points shared by the heatmap stars (`""` /
  `"*"` / `"**"` / `"***"`) and the link-segment linetype legend.

## Value

a ggplot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Environment / Experiment: samples in rows (with a Sample column),
# variables in columns. Edges connect Experiment variables (hubs) to
# any other nodes.
p <- gglink_heatmap_triple(
  Environment = env_df,   # Sample + environmental variables
  Experiment  = exp_df,   # Sample + experiment variables (become hubs)
  edge        = data.frame(from = c("ExpA", "ExpB"),
                           to   = c("pH", "TN"),
                           weight = c(0.8, 0.5)),
  node        = data.frame(node = c("pH", "TN", "ExpA", "ExpB"))
)
} # }
```
