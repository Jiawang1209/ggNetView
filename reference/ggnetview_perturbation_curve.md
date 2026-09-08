# Plot a network perturbation (attack) curve

Renders the perturbation curve produced by
\[get_network_perturbation()\] – a chosen topology metric against the
fraction of nodes removed, coloured by strategy. For the random strategy
a mean +/- standard-error ribbon is drawn. The styling matches the other
\`ggnetview\_\*\` plots (\`theme_classic\`, square aspect, black axes).

## Usage

``` r
ggnetview_perturbation_curve(curve, metric = "LCC_fraction")
```

## Arguments

- curve:

  A data frame as returned in the \`curve\` element of
  \[get_network_perturbation()\], or that element bound from several
  runs (e.g. \`rbind(random\$curve, targeted\$curve)\`) to overlay
  strategies.

- metric:

  Character (default \`"LCC_fraction"\`). Which metric to plot; must be
  present in \`curve\$metric\`.

## Value

A ggplot object.

## See also

\[get_network_perturbation()\].

## Examples

``` r
# \donttest{
data(ppi_example)
obj <- build_graph_from_df(
  df              = ppi_example$ppi,
  node_annotation = ppi_example$annotation
)
rnd <- get_network_perturbation(obj, strategy = "random", bootstrap = 20)
tgt <- get_network_perturbation(obj, strategy = "targeted",
                                centrality = "degree")
ggnetview_perturbation_curve(rbind(rnd$curve, tgt$curve))
#> Warning: Removed 21 rows containing missing values or values outside the scale range
#> (`geom_ribbon()`).

# }
```
