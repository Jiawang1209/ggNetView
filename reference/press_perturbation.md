# Press-perturbation (sustained-disturbance) analysis

Approximates the classic ecological \*press perturbation\*: it treats
the signed correlation matrix as a proxy for the community (interaction)
matrix \\A\\, adds negative self-regulation on the diagonal, and inverts
it to obtain the net-effect matrix \\N = -A^{-1}\\. Entry \\N\_{ij}\\ is
the long-run net response of node \\i\\ when node \\j\\ is held under
sustained pressure (continuously elevated or suppressed). This is the
closest one can get to a "type 3" dynamical perturbation from a static
correlation network.

## Usage

``` r
press_perturbation(
  graph_obj = NULL,
  cor_mat = NULL,
  self_regulation = NULL,
  source = NULL
)
```

## Arguments

- graph_obj:

  A \`tbl_graph\` from any \`build_graph_from\_\*()\` constructor.
  Ignored if \`cor_mat\` is supplied.

- cor_mat:

  Optional numeric matrix. A signed correlation / interaction matrix
  with matching row/col names, used directly instead of extracting one
  from \`graph_obj\`.

- self_regulation:

  Numeric scalar, or \`NULL\` (default). The diagonal of \`A\`
  (intraspecific density dependence; must be negative). When \`NULL\`,
  it is set automatically to \`-(max Re eigenvalue(off-diagonal A) +
  1)\`, which guarantees a stable matrix. Supply your own (e.g. \`-1\`)
  to encode a specific assumption.

- source:

  Optional character vector of node \`name\`s. When given, \`response\`
  returns only the net effect on every node of pressing these source
  node(s).

## Value

A list with:

- \`net_effect\`: the \\N = -A^{-1}\\ matrix (columns = pressed node,
  rows = responding node).

- \`stable\`: logical; \`TRUE\` if \`A\` is dynamically stable.

- \`eigen_real_max\`: largest real part of \`A\`'s eigenvalues (must be
  \`\< 0\` for stability).

- \`self_regulation\`: the diagonal value actually used.

- \`response\`: if \`source\` was given, a data frame of the net
  response of each node to pressing the source(s); else \`NULL\`.

## Assumptions and honest limits

This is a deliberately approximate, \*qualitative\* method. Correlation
is \*\*not\*\* causation and carries no direction, so the interaction
signs and magnitudes are proxies, not measured coefficients. The
framework also requires the community matrix to be \*\*dynamically
stable\*\* (all eigenvalues of \`A\` have negative real part);
\`press_perturbation()\` checks this and warns when it fails. Treat the
output as a defensible qualitative scenario ("if I keep suppressing
taxon A, the community tends to shift this way"), never as a
quantitative prediction. For a purely structural read with no stability
assumption, use \[get_network_perturbation()\] or
\[get_node_influence()\].

## References

Bender EA, Case TJ, Gilpin ME (1984). "Perturbation experiments in
community ecology: theory and practice." *Ecology* 65(1):1-13. May RM
(1972). "Will a large complex system be stable?" *Nature* 238:413-414.
Novak M et al. (2016). "Characterizing species interactions to
understand press perturbations." *Annu. Rev. Ecol. Evol. Syst.*
47:409-432.

## See also

\[get_network_perturbation()\], \[get_node_influence()\].

## Examples

``` r
# \donttest{
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
pp <- press_perturbation(obj)
pp$stable
#> [1] TRUE
pp$net_effect[1:3, 1:3]
#>             C13         C28         C2
#> C13 0.009625942 0.002331755 0.00000000
#> C28 0.002331755 0.009625942 0.00000000
#> C2  0.000000000 0.000000000 0.01023473
# }
```
