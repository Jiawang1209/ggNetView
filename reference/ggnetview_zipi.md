# Compute Zi-Pi (within-module connectivity and participation coefficient)

Calculates the within-module degree z-score (Zi) and among-module
connectivity (participation coefficient, Pi) for each node in a modular
network. These metrics classify nodes into roles such as module hubs,
connectors, and peripherals.

## Usage

``` r
ggnetview_zipi(
  nodes_bulk,
  z_bulk_mat,
  modularity_col,
  degree_col,
  zi_threshold = 2.5,
  pi_threshold = 0.62,
  na.rm = FALSE
)
```

## Arguments

- nodes_bulk:

  Data frame or tibble. Node table with modularity and degree
  information. Node IDs must be in `rownames` or in a `name` column
  (compatible with
  [`tidygraph::as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)
  output).

- z_bulk_mat:

  Numeric matrix. Adjacency or correlation matrix; rows and columns must
  correspond to nodes. Non-zero entries are treated as edges. `NA`/`Inf`
  are replaced with 0.

- modularity_col:

  Character. Column name in `nodes_bulk` containing module labels.

- degree_col:

  Character. Column name in `nodes_bulk` containing node degree (number
  of edges).

- zi_threshold:

  Numeric (default = 2.5). Threshold for within-module connectivity (Zi)
  in role classification.

- pi_threshold:

  Numeric (default = 0.62). Threshold for among-module connectivity (Pi)
  in role classification.

- na.rm:

  Logical (default = `FALSE`). If `TRUE`, remove rows with NA in Zi or
  Pi from the output. If `FALSE`, keep all rows; NA in Zi/Pi results in
  `type = NA`.

## Value

A list with two elements:

- `data`: Data frame merging `nodes_bulk` with
  `within_module_connectivities`, `among_module_connectivities`, and
  `type` (node role).

- `plot`: ggplot object of the Zi-Pi scatter plot with quadrant labels
  and background shading.

## Details

**Zi (within-module connectivity):** Reflects how strongly a node is
connected within its own module. Higher values indicate the node has
more connections within the module and may play a core role inside it.

**Pi (among-module connectivity):** Measures how much a node connects to
other modules. Higher values indicate the node acts as a bridge between
modules, facilitating information, material or energy flow across the
network.

**Node roles (by default thresholds Zi=2.5, Pi=0.62):**

- `Module hubs`: High Zi, low Pi. Core members within their module,
  important for module stability and function, but weakly connected to
  other modules.

- `Connectors`: Low Zi, high Pi. Not prominent within their module, but
  strongly connect across modules, acting as bridges.

- `Network hubs`: High Zi, high Pi. Core nodes both within and across
  modules, critical for overall network structure and stability.

- `Peripherals`: Low Zi, low Pi. Peripheral or satellite nodes with few
  connections within and across modules.

## References

Guimera R, Amaral LAN (2005). "Functional cartography of complex
metabolic networks." *Nature* 433(7028):895–900.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- build_graph_from_mat(otu_rare_relative, method = "WGCNA",
  transfrom.method = "none", cor.method = "pearson", proc = "Bonferroni",
  module.method = "Fast_greedy")
nodes_bulk <- get_graph_nodes(g)
adj_mat <- get_graph_adjacency(g)
res <- ggnetview_zipi(nodes_bulk, adj_mat, "Modularity", "Degree")
res$data   # Zi-Pi 计算结果
res$plot   # Zi-Pi 散点图
} # }
```
