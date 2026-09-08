# Build a consensus graph object from multiple adjacency matrices

Aggregates two or more adjacency matrices – typically produced by
different network inference methods (e.g. SparCC, SpiecEasi, WGCNA,
Hmisc, plain correlation) on the same set of features – into a single
consensus network and returns a \`tbl_graph\` consistent with the rest
of the \`build_graph_from\_\*()\` family. Combining multiple methods is
a standard strategy for reducing method-specific bias in biological
network inference; see Aghayeva et al. (2024, CMiNet) and Chowdhury et
al. (2024, hybrid Bayesian / ML / network framework) for related
approaches.

## Usage

``` r
build_graph_from_consensus(
  adj_list,
  method = c("rank_fusion", "intersection", "weighted_average", "majority_vote"),
  rank_fusion_algorithm = c("borda", "rra", "rrf"),
  threshold = NULL,
  binarize = c("none", "threshold", "topk"),
  binarize_threshold = 0,
  binarize_topk = NULL,
  weights = NULL,
  min_methods = NULL,
  rrf_k = 60,
  node_handling = c("intersect", "union"),
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  node_annotation = NULL,
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- adj_list:

  Named or unnamed list of at least two square numeric adjacency
  matrices. All matrices must have row and column names identifying
  features; matrices need not share identical feature sets (see
  \`node_handling\`). Pre-thresholded inputs (zeros for non-significant
  pairs) are supported and recommended.

- method:

  Character. Top-level consensus strategy. One of \`"rank_fusion"\`
  (default), \`"intersection"\`, \`"weighted_average"\`,
  \`"majority_vote"\`. See **Method details**.

- rank_fusion_algorithm:

  Character. Only used when \`method = "rank_fusion"\`. One of
  \`"borda"\` (default), \`"rra"\`, \`"rrf"\`. See **Rank fusion
  algorithms**.

- threshold:

  Numeric or \`NULL\` (default). Optional final absolute-weight cutoff
  applied to the consensus matrix; entries with \`abs(w) \< threshold\`
  are zeroed out.

- binarize:

  Character. Only used by \`"intersection"\` and \`"majority_vote"\`.
  One of \`"none"\` (default; treat any non-zero entry as an edge),
  \`"threshold"\` (apply \`binarize_threshold\`), \`"topk"\` (keep the
  global top-\`binarize_topk\` strongest edges per matrix).

- binarize_threshold:

  Numeric (default \`0\`). Per-method \`\|w\| \>= binarize_threshold\`
  edge cut-off when \`binarize = "threshold"\`.

- binarize_topk:

  Integer or \`NULL\`. Per-method top-K cut-off when \`binarize =
  "topk"\`.

- weights:

  Numeric vector or \`NULL\` (default). Per-method weights for
  \`"weighted_average"\`. Length must match \`length( adj_list)\`;
  values are renormalised to sum to 1. \`NULL\` uses uniform weights.

- min_methods:

  Integer or \`NULL\` (default). Minimum number of methods that must
  support an edge for \`"majority_vote"\`. \`NULL\` defaults to a strict
  majority (\`floor(M/2) + 1\` where \`M\` is the number of methods).

- rrf_k:

  Numeric (default \`60\`). Smoothing constant for reciprocal rank
  fusion (\`method = "rank_fusion", rank_fusion_algorithm = "rrf"\`).

- node_handling:

  Character. How to align features that are not present in every input
  matrix. \`"intersect"\` (default) keeps only features common to all
  inputs; \`"union"\` keeps every feature appearing in any input and
  treats missing entries as zero.

- module.method:

  Character. Module detection method passed through to
  \[build_graph_from_adj_mat()\] when constructing the final
  \`tbl_graph\`. One of \`"Fast_greedy"\`, \`"Walktrap"\`,
  \`"Edge_betweenness"\`, \`"Spinglass"\`.

- node_annotation:

  Optional data frame attached as vertex metadata; first column must
  match feature names of the consensus matrix.

- top_modules:

  Integer (default \`15\`). Number of top-ranked modules to retain in
  the final \`tbl_graph\`; smaller modules are collapsed into
  \`"Others"\`.

- seed:

  Integer (default \`1115\`). Random seed.

## Value

A \`tbl_graph\` object whose schema matches the rest of the
\`build_graph_from\_\*()\` family (vertex columns \`name\`,
\`modularity\`, \`modularity2\`, \`modularity3\`, \`Modularity\`,
\`Degree\`, \`Strength\`; edge columns \`weight\`, \`correlation\`,
\`corr_direction\`). Plug straight into \[ggNetView()\] for
visualisation.

## Details

Four consensus strategies are supported through \`method\`. Three of
them (\`"intersection"\`, \`"weighted_average"\`, \`"majority_vote"\`)
operate directly on the input weights; the fourth (\`"rank_fusion"\`)
converts each method's adjacency to ranks first, which makes the
aggregation invariant to the very different weight scales produced by
different inference algorithms (correlation in \[-1, 1\], partial
correlation, TOM in \[0, 1\], ...). When \`method = "rank_fusion"\`,
three classical rank-aggregation algorithms are available via
\`rank_fusion_algorithm\`: Borda count (mean rank), Robust Rank
Aggregation (RRA, Kolde et al. 2012), and Reciprocal Rank Fusion (RRF,
Cormack et al. 2009).

## Method details

- \`"intersection"\`:

  An edge survives in the consensus iff it is non-zero (after
  \`binarize\`) in \*\*every\*\* input matrix. Consensus weight is the
  mean of the input weights for surviving edges. Strictest of the four
  strategies.

- \`"weighted_average"\`:

  Each input matrix is min-max normalised to \[0, 1\] in absolute value
  (sign preserved), then combined as a weighted sum using \`weights\`
  (default uniform). Consensus weight is in \[-1, 1\] and inherits the
  sign of the methods that dominate the average.

- \`"majority_vote"\`:

  Each input matrix is binarised, edges are summed, and an edge survives
  in the consensus iff the count is at least \`min_methods\`. Consensus
  weight is the mean of the input weights restricted to the methods that
  voted yes.

- \`"rank_fusion"\`:

  Each input matrix's pairs are ranked by \`\|w\|\` (descending). Ranks
  are aggregated according to \`rank_fusion_algorithm\`, then min-max
  normalised to \[-1, 1\] (sign taken from the mean of the input
  weights). Most robust to across-method scale differences.

## Rank fusion algorithms

- \`"borda"\`:

  Mean rank across methods. Lower mean rank = more supported edge.
  Simple, no assumptions, good baseline.

- \`"rra"\`:

  Robust Rank Aggregation (Kolde et al. 2012, Bioinformatics). Computes
  a rho-score for each pair under the null hypothesis of uniformly
  distributed ranks; edges with very small rho consistently outrank
  random expectation. Uses the \`RobustRankAggreg\` package when
  available; otherwise falls back to an inline Beta-order-statistic
  implementation. Output score is \`-log10(rho)\` so larger = more
  confident.

- \`"rrf"\`:

  Reciprocal Rank Fusion (Cormack et al. 2009). Per-edge score is
  \`sum_m 1 / (k + rank_m)\` where \`k = rrf_k\` (default 60, the
  information-retrieval convention). Robust to extreme rank values, no
  parameter tuning required.

## References

Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). Robust rank
aggregation for gene list integration and meta-analysis.
\*Bioinformatics\*, 28(4), 573-580.

Cormack, G. V., Clarke, C. L. A., & Buettcher, S. (2009). Reciprocal
rank fusion outperforms Condorcet and individual rank learning methods.
\*Proceedings of the 32nd International ACM SIGIR Conference\*, 758-759.

Aghayeva, R., et al. (2024). CMiNet: An R package and user-friendly
Shiny App for constructing consensus microbiome networks.

Chowdhury, S., et al. (2024). A hybrid framework for disease biomarker
discovery in microbiome research combining Bayesian networks, machine
learning, and network-based methods.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 30
nm <- paste0("g", seq_len(n))

# Three "methods": shared signal + per-method noise.
core <- matrix(0, n, n, dimnames = list(nm, nm))
for (i in 1:6) for (j in (i + 1):8) core[i, j] <- core[j, i] <- 0.8
jitter <- function() core + matrix(stats::rnorm(n * n, sd = 0.05), n, n,
                                   dimnames = list(nm, nm))

adj_list <- list(method_A = jitter(), method_B = jitter(), method_C = jitter())

# Borda rank fusion with a final |w| >= 0.3 cutoff:
obj <- build_graph_from_consensus(
  adj_list = adj_list,
  method = "rank_fusion",
  rank_fusion_algorithm = "borda",
  threshold = 0.3
)
#> The max module in network is 2 we use the 2  modules for next analysis
obj
#> # A tbl_graph: 30 nodes and 330 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 30 × 7 (active)
#>    name  modularity modularity2 modularity3 Modularity Degree Strength
#>    <chr> <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl>
#>  1 g9    2          2           2           2              27     13.6
#>  2 g28   2          2           2           2              27     13.2
#>  3 g20   2          2           2           2              26     11.7
#>  4 g16   2          2           2           2              25     12.4
#>  5 g23   2          2           2           2              24     11.8
#>  6 g22   2          2           2           2              23     11.2
#>  7 g29   2          2           2           2              23     11.0
#>  8 g13   2          2           2           2              22     10.6
#>  9 g24   2          2           2           2              22     10.8
#> 10 g27   2          2           2           2              22     11.6
#> # ℹ 20 more rows
#> #
#> # Edge Data: 330 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1    23    24  0.991       0.991 Positive      
#> 2    20    23  0.979       0.979 Positive      
#> 3    21    23  0.955       0.955 Positive      
#> # ℹ 327 more rows
# }
```
