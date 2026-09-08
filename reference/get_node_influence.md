# Propagate a virtual perturbation from source node(s) across the network

Injects a virtual perturbation at one or more \`source\` nodes and lets
it spread along the (optionally signed) weighted edges, returning how
strongly every other node is affected. This is the abundance-influence
("type 2") virtual-perturbation analysis: it treats edge weights as
interaction strengths and asks "if I nudge species A, how far and how
strongly does the ripple reach?".

## Usage

``` r
get_node_influence(
  graph_obj,
  source,
  delta = 1,
  alpha = 0.5,
  signed = TRUE,
  drop_source = TRUE,
  overwrite = TRUE
)
```

## Arguments

- graph_obj:

  A \`tbl_graph\` from any \`build_graph_from\_\*()\` constructor.

- source:

  Character vector of node \`name\`s to perturb.

- delta:

  Numeric (default \`1\`). Magnitude of the injected perturbation placed
  on each source node.

- alpha:

  Numeric in \`(0, 1)\` (default \`0.5\`). Diffusion decay; the function
  automatically caps it just below \`1 / spectral-radius(W)\` so the
  series converges.

- signed:

  Logical (default \`TRUE\`). Use the signed \`correlation\` edge
  attribute (so anticorrelated neighbours receive negative influence).
  \`FALSE\` uses \`\|weight\|\` only.

- drop_source:

  Logical (default \`TRUE\`). Zero out the source node(s)' own influence
  in the returned column so the ranking reflects downstream spread only.

- overwrite:

  Logical (default \`TRUE\`). If an \`Influence\` column already exists
  on the input graph, controls whether to overwrite it (silent overwrite
  when \`TRUE\`; warning + return unchanged when \`FALSE\`).

## Value

The input \`tbl_graph\` with one new node column, \`Influence\` (signed
when \`signed = TRUE\`). Larger magnitude = more strongly affected. Map
it straight onto a figure with \`ggNetView(..., node_fill =
"Influence")\`.

## Details

Propagation uses a Katz / random-walk-with-restart diffusion,
\\influence = (I - \alpha W)^{-1} s\\, where \`W\` is the
column-normalised (signed) adjacency, \`s\` places \`delta\` on the
source node(s), and \`alpha\` is the decay. This always converges
(\`alpha\` is capped below the inverse spectral radius) and is
well-defined even on disconnected graphs.

## Important interpretation note

Correlation / co-occurrence networks encode \*\*association, not
causation\*\*, and give no edge direction. The score returned here is a
\*structural influence estimate\* – a weighted measure of how reachable
each node is from the source – and should be read as a qualitative
ranking, \*\*not\*\* as a quantitative ecological-dynamics prediction.
For a perturbation read with a (still approximate) mechanistic flavour,
see \[press_perturbation()\].

## See also

\[get_network_perturbation()\] for structural attacks;
\[press_perturbation()\] for the press-perturbation approximation.

## Examples

``` r
# \donttest{
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
src <- get_graph_nodes(obj)$name[1]
obj2 <- get_node_influence(obj, source = src)
obj2 %>%
  tidygraph::activate(nodes) %>%
  tidygraph::as_tibble() %>%
  dplyr::arrange(dplyr::desc(abs(Influence))) %>%
  utils::head(5)
#> # A tibble: 5 × 10
#>   name  group modularity modularity2 modularity3 Modularity Degree Segree
#>   <chr> <chr> <fct>      <fct>       <chr>       <fct>       <dbl>  <dbl>
#> 1 C28   C     1          1           1           1               1      1
#> 2 C13   C     1          1           1           1               1      1
#> 3 C2    C     10         10          10          10              1      1
#> 4 D9    D     10         10          10          10              1      1
#> 5 A3    A     11         11          11          11              1      1
#> # ℹ 2 more variables: Strength <dbl>, Influence <dbl>
# }
```
