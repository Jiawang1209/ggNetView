# Build a graph object from a STRING database export

Loads a STRING-formatted protein-protein interaction file (the TSV
export with columns \`node1, node2, ..., combined_score\`) or a
user-supplied data frame following the same schema, and converts it into
a \`tbl_graph\` consistent with the rest of the
\`build_graph_from\_\*()\` family. All STRING evidence channels in the
input – \`homology\`, \`coexpression\`,
\`experimentally_determined_interaction\`, \`automated_textmining\`,
etc. – are preserved as edge attributes so downstream visualizations can
colour or filter edges by evidence type.

## Usage

``` r
build_graph_from_stringdb(
  stringdb,
  node1_col = "node1",
  node2_col = "node2",
  score_col = "combined_score",
  score_threshold = NULL,
  node_annotation = NULL,
  directed = FALSE,
  module.method = c("Fast_greedy", "Walktrap", "Edge_betweenness", "Spinglass"),
  top_modules = 15,
  seed = 1115
)
```

## Arguments

- stringdb:

  Either a file path (\`.tsv\` / \`.tsv.gz\`) or a data.frame already
  loaded in memory. STRING-style headers (a leading \`#\` on the first
  column name) are handled automatically.

- node1_col, node2_col, score_col:

  Character. Column names to use as the edge \`from\`, \`to\`, and
  weight columns. Defaults match STRING's TSV layout (\`"node1"\`,
  \`"node2"\`, \`"combined_score"\`). Override if you want to weight by
  an individual evidence channel – e.g. \`score_col =
  "experimentally_determined_interaction"\`.

- score_threshold:

  Numeric or \`NULL\` (default). If non-NULL, edges with \`score_col \<
  score_threshold\` are dropped. STRING confidence cut-offs are
  typically \`0.4\` (medium), \`0.7\` (high), \`0.9\` (highest). Pass
  \`NULL\` to keep every edge in the input.

- node_annotation:

  Optional data frame; first column must match the protein/gene names
  appearing in \`node1_col\` / \`node2_col\`. If provided, these
  annotations are attached as vertex attributes via
  \[build_graph_from_df()\].

- directed:

  Logical (default: \`FALSE\`). STRING interactions are undirected by
  convention.

- module.method:

  Character. Network community detection method. Options:
  \`"Fast_greedy"\`, \`"Walktrap"\`, \`"Edge_betweenness"\`,
  \`"Spinglass"\`.

- top_modules:

  Integer. Number of top-ranked modules to retain; smaller modules are
  collapsed into \`"Others"\`.

- seed:

  Integer (default: \`1115\`). Random seed for reproducibility.

## Value

A \`tbl_graph\` object compatible with \[ggNetView()\]. Edge attributes
include \`weight\` (the chosen score), \`correlation\` (signed copy of
the score, equal to \`weight\` for STRING since scores are
non-negative), and every other column present in the input (e.g.
\`combined_score\`, \`coexpression\`, \`homology\`,
\`experimentally_determined_interaction\`).

## Details

STRING's TSV exports use a header line that starts with \`#node1\`. The
leading \`#\` is stripped automatically; all other column names are kept
verbatim.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a STRING TSV file:
obj <- build_graph_from_stringdb(
  stringdb        = "string_interactions.tsv",
  score_threshold = 0.7
)

# From an in-memory data.frame, no filtering:
df <- data.frame(
  node1          = c("AANAT", "AANAT", "ABCA1"),
  node2          = c("CRY1",  "TPH1",  "SIRT1"),
  combined_score = c(0.608,   0.675,   0.520)
)
obj <- build_graph_from_stringdb(stringdb = df)
} # }
```
