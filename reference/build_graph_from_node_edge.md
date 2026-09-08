# Build a graph object from separate node and edge tables

Provides an explicit two-table interface where the user supplies a node
table and an edge table separately. Compared with
\[build_graph_from_df()\], which takes the edge list as the primary
input and treats the node annotation as optional, here the node table is
\*\*required and authoritative\*\*: nodes that have no incident edges
are kept in the resulting graph (no orphan-node removal). This is the
right entry point when the user already has a curated vertex set and
wants every node to appear in the visualization, even if it ends up
isolated after edge filtering.

## Usage

``` r
build_graph_from_node_edge(
  node,
  edge,
  directed = FALSE,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- node:

  Data frame. Node table whose \*\*first column is the node
  identifier\*\* that matches the \`from\` / \`to\` values in \`edge\`.
  Remaining columns are attached as vertex attributes (taxonomy, gene
  metadata, ...).

- edge:

  Data frame. Edge list with at least two columns; the \*\*first two
  columns are interpreted as \`from\` and \`to\`\*\* (column names are
  not enforced). An optional \`weight\` column (any case) is used as
  edge weight; if no \`weight\`-named column is found and \`ncol(edge)
  \>= 3\`, the third column is used.

- directed:

  Logical (default: \`FALSE\`). Whether edges between nodes are
  directed.

- module.method:

  Character. Network community detection method. Options:
  \`"Fast_greedy"\`, \`"Walktrap"\`, \`"Edge_betweenness"\`,
  \`"Spinglass"\`. \`"Fast_greedy"\` and \`"Spinglass"\` require an
  undirected graph.

- top_modules:

  Integer. Number of top-ranked modules (by node count) to retain;
  smaller modules are collapsed into \`"Others"\`.

- seed:

  Integer (default: \`1115\`). Random seed for reproducibility.

## Value

A \`tbl_graph\` object with vertex columns \`name\`, \`modularity\`,
\`modularity2\`, \`modularity3\`, \`Modularity\`, \`Degree\`,
\`Segree\`, \`Strength\` plus any extra columns supplied via \`node\`,
and edge columns \`weight\` (positive), \`correlation\` (signed) plus
any extra columns supplied via \`edge\`.

## Examples

``` r
if (FALSE) { # \dontrun{
node <- data.frame(
  name  = c("A", "B", "C", "D", "E"),
  group = c("g1", "g1", "g2", "g2", "g3")
)
edge <- data.frame(
  from   = c("A", "B", "C"),
  to     = c("B", "C", "D"),
  weight = c(0.8, 0.6, 0.7)
)
# Note: node "E" is isolated. Unlike build_graph_from_df(), the result
# of build_graph_from_node_edge() retains "E".
obj <- build_graph_from_node_edge(node = node, edge = edge)
obj
} # }
```
