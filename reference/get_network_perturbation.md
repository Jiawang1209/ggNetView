# Structural perturbation ("virtual attack") of a network

Repeatedly removes nodes from a network – at random, by targeted attack
on a centrality, or by knocking out a named module / node set –
recomputes a panel of connectivity-sensitive topology metrics after each
removal, and returns the resulting \*perturbation curve\* together with
a single-number robustness index (Schneider R). This is the structural
("type 1") virtual-perturbation analysis: it asks how the network falls
apart as nodes are progressively lost.

## Usage

``` r
get_network_perturbation(
  graph_obj,
  strategy = c("random", "targeted", "module", "manual"),
  centrality = c("degree", "strength", "betweenness", "closeness", "eigenvector", "ivi"),
  target = NULL,
  module_col = "Modularity",
  fractions = seq(0.05, 1, by = 0.05),
  decreasing = TRUE,
  bootstrap = 100,
  seed = 123,
  plot = TRUE
)
```

## Arguments

- graph_obj:

  A \`tbl_graph\` from any \`build_graph_from\_\*()\` constructor.

- strategy:

  Character. One of \`"random"\`, \`"targeted"\`, \`"module"\`,
  \`"manual"\`.

- centrality:

  Character. Used when \`strategy = "targeted"\`. One of \`"degree"\`,
  \`"strength"\`, \`"betweenness"\`, \`"closeness"\`, \`"eigenvector"\`,
  \`"ivi"\`. \`"strength"\` uses the \`weight\` edge attribute;
  \`"ivi"\` requires the \`influential\` package.

- target:

  Character vector. The module label(s) (when \`strategy = "module"\`)
  or node \`name\`s (when \`strategy = "manual"\`) to remove.

- module_col:

  Character (default \`"Modularity"\`). Node column holding module
  labels, used when \`strategy = "module"\`.

- fractions:

  Numeric vector in \`(0, 1\]\`. Removal fractions for \`"random"\` /
  \`"targeted"\`. Default \`seq(0.05, 1, by = 0.05)\`.

- decreasing:

  Logical (default \`TRUE\`). For \`"targeted"\`, remove the
  most-central nodes first (\`TRUE\`) or least-central first
  (\`FALSE\`).

- bootstrap:

  Integer (default \`100\`). Number of random repetitions for \`strategy
  = "random"\`.

- seed:

  Integer (default \`123\`). Seed for the random strategy, so results
  are reproducible.

- plot:

  Logical (default \`TRUE\`). Attach a ready-made attack-curve ggplot of
  \`LCC_fraction\` (only for \`"random"\` / \`"targeted"\`).

## Value

A list with:

- \`curve\`: long data frame (\`strategy\`, \`fraction\`, \`metric\`,
  \`value\`, and for random \`value_sd\` / \`value_se\`).

- \`robustness_index\`: data frame with the Schneider R-index (area
  under the \`LCC_fraction\` curve) per strategy.

- \`plot\`: ggplot of the LCC attack curve, or \`NULL\`.

## Details

The random strategy is the multi-metric generalisation of the
node-removal robustness already computed inside
\[get_network_topology()\]; the targeted / module / manual strategies
answer "which nodes (or whole modules) hold the network together?".

## Strategies

- \`"random"\`:

  Remove a random subset at each fraction, repeated \`bootstrap\` times;
  mean / sd / se are reported. Seeded for reproducibility.

- \`"targeted"\`:

  Rank nodes by \`centrality\` and remove them in order (most-central
  first by default). The classic intentional attack – usually far more
  damaging than random failure.

- \`"module"\`:

  Knock out every node whose \`module_col\` label is in \`target\`;
  reported as a before/after comparison.

- \`"manual"\`:

  Knock out the exact node \`name\`s given in \`target\`; reported as a
  before/after comparison.

## Metrics tracked after each removal

- \`LCC_fraction\`:

  Size of the largest connected component as a fraction of the
  \*original\* node count (the main attack curve).

- \`N_components\`:

  Number of connected components.

- \`Natural_connectivity\`:

  \`log(mean(exp(eigenvalues(A))))\`; a spectral robustness measure that
  varies smoothly and does not jump discretely the way component counts
  do.

- \`Efficiency\`:

  Mean of \`1 / shortest-path-distance\`.

- \`Mean_degree\`, \`Density\`, \`Transitivity_global\`, \`Modularity\`:

  Standard summaries recomputed on the survivor subgraph.

## References

Albert R, Jeong H, Barabasi AL (2000). "Error and attack tolerance of
complex networks." *Nature* 406:378-382. Schneider CM et al. (2011).
"Mitigation of malicious attacks on networks." *PNAS* 108(10):3838-3841.

## See also

\[get_node_influence()\] for abundance-influence propagation;
\[press_perturbation()\] for the press-perturbation approximation;
\[get_network_topology()\] for static network-level metrics.

## Examples

``` r
# \donttest{
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
res <- get_network_perturbation(obj, strategy = "targeted",
                                centrality = "degree")
head(res$curve)
#>   strategy fraction               metric        value value_sd value_se
#> 1 targeted        0         LCC_fraction 2.000000e-02       NA       NA
#> 2 targeted        0         N_components 5.000000e+01       NA       NA
#> 3 targeted        0 Natural_connectivity 4.337808e-01       NA       NA
#> 4 targeted        0           Efficiency 1.882805e-04       NA       NA
#> 5 targeted        0          Mean_degree 1.000000e+00       NA       NA
#> 6 targeted        0              Density 1.010101e-02       NA       NA
res$robustness_index
#>   strategy    R_index
#> 1 targeted 0.01571429
# }
```
