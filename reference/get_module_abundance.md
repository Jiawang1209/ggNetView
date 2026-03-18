# Compute module relative abundance index

For each module, computes the sum or mean of member OTU abundances per
sample. Simpler than eigengene but more "abundance-focused".

## Usage

``` r
get_module_abundance(
  otu_mat,
  graph_obj,
  module_col = "Modularity",
  type = c("sum", "mean"),
  exclude_others = TRUE
)
```

## Arguments

- otu_mat:

  Numeric matrix. Rows = OTUs, columns = samples.

- graph_obj:

  A tbl_graph with `name` and module column.

- module_col:

  Character. Module column name.

- type:

  Character. `"sum"` (total abundance) or `"mean"` (average abundance
  per OTU).

- exclude_others:

  Logical. Exclude "Others" module.

## Value

A data frame with samples as rows and modules as columns.
