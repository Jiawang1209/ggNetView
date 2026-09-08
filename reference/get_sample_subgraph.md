# Extract sample-level subgraphs from a graph object

For each sample (column in `mat`), OTUs whose abundance exceeds
`min_abundance` are treated as present in that sample. The induced
subgraph of `graph_obj` on those OTUs is returned. Optionally, multiple
selected samples can be combined into a single merged subgraph via union
or intersection of their present OTUs.

## Usage

``` r
get_sample_subgraph(
  graph_obj,
  mat,
  min_abundance = 0,
  select_sample = NULL,
  combine = c("union", "intersect")
)
```

## Arguments

- graph_obj:

  A `tbl_graph` object from
  [`build_graph_from_mat`](https://jiawang1209.github.io/ggNetView/reference/build_graph_from_mat.md)
  or
  [`build_graph_from_df`](https://jiawang1209.github.io/ggNetView/reference/build_graph_from_df.md).
  Its node table must contain a `name` column matching `rownames(mat)`.

- mat:

  Numeric matrix. Rows are OTUs / features (must have rownames matching
  graph node names); columns are samples (must have colnames as sample
  IDs).

- min_abundance:

  Numeric (default = 0). An OTU is considered present in a sample when
  `mat[OTU, sample] > min_abundance`. The appropriate value depends on
  the scale of `mat`: e.g. `0` for raw or rarefied counts, `0.001` for
  relative abundance. The function does not infer data type; the user is
  responsible for choosing a meaningful threshold for their data.

- select_sample:

  Character vector (default = `NULL`). Sample IDs to extract into a
  single merged subgraph. Must be a subset of `colnames(mat)`. When
  `NULL`, only per-sample subgraphs are returned and `sub_graph_select`
  is `NULL`.

- combine:

  Character. One of `"union"` (default) or `"intersect"`. Controls how
  OTUs from the selected samples are combined when building
  `sub_graph_select`:

  - `"union"`: keep OTUs present in any of the selected samples.

  - `"intersect"`: keep only OTUs present in all of the selected
    samples.

  Edges are always the induced edges from `graph_obj` between the
  surviving nodes.

## Value

A list with three elements:

- `sub_graph_all`: named list of per-sample `tbl_graph` subgraphs.
  Samples with no present OTUs in the graph are dropped from this list
  (but are still recorded in `stat_sample`).

- `stat_sample`: data frame with columns `Sample`, `Node`, `Edge`,
  `Status`. One row per sample in `mat`, including samples with empty
  subgraphs.

- `sub_graph_select`: a single merged `tbl_graph` of the nodes combined
  from `select_sample` according to `combine`, or `NULL` if
  `select_sample` is `NULL` or the combination produces an empty node
  set. The node table of this graph carries two additional columns:

  - `n_present_samples` (integer): how many of the selected samples this
    node appears in.

  - `present_in_samples` (character): comma-separated sample IDs the
    node appears in, ordered following `select_sample`.

## Details

This function is the sample-wise analogue of
[`get_subgraph`](https://jiawang1209.github.io/ggNetView/reference/get_subgraph.md),
which splits a graph by its `Modularity` node attribute. Unlike module
membership (which is a partition), an OTU can belong to multiple
samples, so a `combine` switch is exposed.

## Examples

``` r
# \donttest{
data("otu_rare_relative")
data("tax_tab")
obj <- build_graph_from_mat(
  mat              = otu_rare_relative,
  transfrom.method = "none",
  r.threshold      = 0.7,
  p.threshold      = 0.05,
  method           = "WGCNA",
  cor.method       = "pearson",
  proc             = "bonferroni",
  module.method    = "Fast_greedy",
  node_annotation  = tax_tab,
  top_modules      = 15,
  seed             = 1115
)

# All per-sample subgraphs + per-sample stats
res <- get_sample_subgraph(graph_obj = obj, mat = otu_rare_relative)
head(res$stat_sample)
#>   Sample Node Edge Status
#> 1    KO1  131  187     OK
#> 2    KO2  141  194     OK
#> 3    KO3  116  404     OK
#> 4    KO4  110  126     OK
#> 5    KO5   85   73     OK
#> 6    KO6  114  146     OK

# Merged subgraph for 3 samples, union of present OTUs
res_u <- get_sample_subgraph(
  graph_obj     = obj,
  mat           = otu_rare_relative,
  select_sample = colnames(otu_rare_relative)[1:3],
  combine       = "union"
)
res_u$sub_graph_select
#> # A tbl_graph: 201 nodes and 660 edges
#> #
#> # An undirected simple graph with 36 components
#> #
#> # Node Data: 201 × 16 (active)
#>    name    modularity modularity2 modularity3 Modularity Degree Strength Kingdom
#>    <chr>   <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl> <chr>  
#>  1 ASV_12… 6          6           6           6              29     28.3 Bacter…
#>  2 ASV_14… 6          6           6           6              28     27.4 Bacter…
#>  3 ASV_649 6          6           6           6              27     26.5 Bacter…
#>  4 ASV_705 6          6           6           6              27     26.5 Bacter…
#>  5 ASV_913 6          6           6           6              27     26.2 Bacter…
#>  6 ASV_13… 6          6           6           6              27     26.5 Bacter…
#>  7 ASV_14… 6          6           6           6              27     26.5 Bacter…
#>  8 ASV_17… 6          6           6           6              27     26.3 Bacter…
#>  9 ASV_24… 6          6           6           6              27     26.5 Bacter…
#> 10 ASV_25… 6          6           6           6              27     26.4 Bacter…
#> # ℹ 191 more rows
#> # ℹ 8 more variables: Phylum <chr>, Class <chr>, Order <chr>, Family <chr>,
#> #   Genus <chr>, Species <chr>, n_present_samples <int>,
#> #   present_in_samples <chr>
#> #
#> # Edge Data: 660 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1   151   152  0.959       0.959 Positive      
#> 2   151   153  0.941       0.941 Positive      
#> 3   151   154  0.940       0.940 Positive      
#> # ℹ 657 more rows

# Same selection but intersection (core OTUs across the 3 samples)
res_i <- get_sample_subgraph(
  graph_obj     = obj,
  mat           = otu_rare_relative,
  select_sample = colnames(otu_rare_relative)[1:3],
  combine       = "intersect"
)
res_i$sub_graph_select
#> # A tbl_graph: 63 nodes and 43 edges
#> #
#> # An undirected simple graph with 30 components
#> #
#> # Node Data: 63 × 16 (active)
#>    name    modularity modularity2 modularity3 Modularity Degree Strength Kingdom
#>    <chr>   <fct>      <ord>       <chr>       <ord>       <dbl>    <dbl> <chr>  
#>  1 ASV_12… 6          6           6           6              29   28.3   Bacter…
#>  2 ASV_913 6          6           6           6              27   26.2   Bacter…
#>  3 ASV_767 6          6           6           6              22   21.1   Bacter…
#>  4 ASV_322 6          6           6           6              19   18.1   Bacter…
#>  5 ASV_18… 6          6           6           6              19   18.2   Bacter…
#>  6 ASV_277 6          6           6           6               4    3.76  Bacter…
#>  7 ASV_244 6          6           6           6               3    2.87  Bacter…
#>  8 ASV_367 6          6           6           6               3    2.83  Bacter…
#>  9 ASV_132 6          6           6           6               1    0.941 Bacter…
#> 10 ASV_927 10         10          10          10             16   15.6   Bacter…
#> # ℹ 53 more rows
#> # ℹ 8 more variables: Phylum <chr>, Class <chr>, Order <chr>, Family <chr>,
#> #   Genus <chr>, Species <chr>, n_present_samples <int>,
#> #   present_in_samples <chr>
#> #
#> # Edge Data: 43 × 5
#>    from    to weight correlation corr_direction
#>   <int> <int>  <dbl>       <dbl> <chr>         
#> 1     1     2  0.986       0.986 Positive      
#> 2     1     3  0.965       0.965 Positive      
#> 3     2     3  0.967       0.967 Positive      
#> # ℹ 40 more rows
# }
```
