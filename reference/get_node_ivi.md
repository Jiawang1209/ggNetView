# Compute the Integrated Value of Influence (IVI) per node

Computes the Integrated Value of Influence (IVI) of Salavaty et al.
(2020) and attaches the per-node score as a new vertex column on the
input \`tbl_graph\`. IVI integrates six centrality measures across three
levels – local (degree, local H-index), semi-local (ClusterRank), and
global (betweenness, closeness, collective influence) – and is more
robust than any single centrality for identifying influential nodes in
complex biological networks (microbiome co-occurrence, gene
co-expression, PPI, ...).

## Usage

``` r
get_node_ivi(
  graph_obj,
  weights = NULL,
  mode = c("all", "in", "out"),
  directed = FALSE,
  d = 3,
  scale = c("range", "z-scale", "none"),
  ncores = 1L,
  overwrite = TRUE
)
```

## Arguments

- graph_obj:

  A \`tbl_graph\` produced by any \`build_graph_from\_\*()\`
  constructor.

- weights:

  Numeric vector or \`NULL\` (default). Edge weights to pass through to
  \`influential::ivi()\`. \`NULL\` uses the unweighted IVI definition.
  Pass the string \`"weight"\` to use the graph's existing \`weight\`
  edge attribute as IVI weights (interpreted as distances by
  \`influential::ivi\`).

- mode:

  Character (default \`"all"\`). Edge mode passed to
  \`influential::ivi()\`. One of \`"all"\`, \`"in"\`, \`"out"\`.
  \`"all"\` is the right choice for undirected graphs.

- directed:

  Logical (default \`FALSE\`). Whether the graph should be treated as
  directed. Passed to \`influential::ivi()\`.

- d:

  Integer (default \`3\`). Distance horizon used by the
  collective-influence component. Larger \`d\` looks farther from each
  node when scoring it but increases run-time roughly linearly. The IVI
  paper recommends \`d = 3\` as a good trade-off.

- scale:

  Character (default \`"range"\`). How to scale the IVI values returned
  by \`influential::ivi()\`:

  \`"range"\`

  :   Normalise to \`\[1, 100\]\`. Use this when exploring a single
      network; lets you see the full spread of node influences.

  \`"z-scale"\`

  :   z-score standardisation. Use this when comparing IVI across
      multiple networks, or when you want a defensible numeric threshold
      (\`z \> 1.645\` is a common cutoff for "significantly
      influential").

  \`"none"\`

  :   Return the raw IVI scores with no scaling.

  On older versions of \`influential\` that exposed the boolean
  \`scaled\` argument instead of \`scale\`, this wrapper automatically
  translates: \`"none"\` -\> \`scaled = FALSE\`, anything else -\>
  \`scaled = TRUE\`.

- ncores:

  Integer (default \`1\`). Number of parallel cores for
  \`influential::ivi()\`'s internal cluster-rank step. The default
  (\`1\`) forces serial execution – this is deliberately conservative so
  the wrapper works inside \`R CMD check\` (which caps tests at 2 cores
  via \`\_R_CHECK_LIMIT_CORES\_\`) and inside vignette builds (where
  leftover parallel clusters are a common source of "9 simultaneous
  processes spawned" / "invalid connection" errors). For interactive
  production work on a real network, raising this to e.g. \`4\` or
  \`parallel::detectCores() - 1\` gives a meaningful speed-up.

- overwrite:

  Logical (default \`TRUE\`). If an \`IVI\` column already exists on the
  input graph, controls whether to overwrite it (silent overwrite when
  \`TRUE\`; warning + return unchanged when \`FALSE\`).

## Value

A \`tbl_graph\` whose node table is augmented with a single new column,
\`IVI\`. Other vertex / edge columns are preserved verbatim. Larger
\`IVI\` = more influential.

## Details

This function is a thin wrapper over \`influential::ivi()\`. The
\`influential\` package is the reference implementation maintained by
the IVI paper's first author and is required at runtime; install with
\`install.packages("influential")\` if not already available. Routing
through the canonical implementation guarantees the IVI values match the
published algorithm exactly.

## Forwarded arguments and version compatibility

Different \`influential\` releases have exposed slightly different
argument lists for \`ivi()\`. Recent releases (the one this package is
developed against) accept a \`scale\` string with values \`"range"\` /
\`"z-scale"\` / \`"none"\`; older releases used a \`scaled\` logical
instead. To stay forward- and backward-compatible, this wrapper inspects
\`formals(influential::ivi)\` at call time and translates the
user-facing \`scale\` argument to whichever signature the installed
version expects. Arguments the installed version does not recognise are
silently dropped from the call.

## References

Salavaty, A., Ramialison, M., & Currie, P. D. (2020). Integrated Value
of Influence: An Integrative Method for the Identification of the Most
Influential Nodes within Networks. \*Patterns\*, 1(5), 100052.

## See also

\[get_node_centrality()\] for the underlying per-node centrality panel;
\[ggnetview_zipi()\] for the complementary Zi-Pi role classification
(Guimera & Amaral 2005).

## Examples

``` r
if (FALSE) { # \dontrun{
# `obj` is a tbl_graph from any build_graph_from_*() constructor.
obj_aug <- get_node_ivi(obj)

# Top 10 most influential nodes:
obj_aug %>%
  tidygraph::activate(nodes) %>%
  tidygraph::as_tibble() %>%
  dplyr::arrange(dplyr::desc(IVI)) %>%
  utils::head(10)

# Map IVI to point fill in ggNetView(): the fill aesthetic uses a
# discrete scale, so bin the continuous IVI into ordered quartile
# factor levels first, then pass the bin column as `node_fill`.
obj_plot <- obj_aug %>%
  tidygraph::activate(nodes) %>%
  tidygraph::mutate(IVI_bin = cut(
    IVI,
    breaks = stats::quantile(IVI, probs = seq(0, 1, 0.25), na.rm = TRUE),
    labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
    include.lowest = TRUE,
    ordered_result = TRUE
  ))
ggNetView(obj_plot, layout = "fr", node_fill = "IVI_bin")

# Use z-score scaling when comparing across networks or applying a
# threshold (e.g. "significantly influential" = z > 1.645):
obj_z <- get_node_ivi(obj, scale = "z-scale")
} # }
```
