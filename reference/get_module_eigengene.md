# Compute module eigengene (PC1) from OTU abundance matrix

For each module, performs PCA on the module members' abundance matrix
(samples x OTUs) and returns the first principal component as the module
eigengene. This represents the "overall activity/state" of each module
across samples.

## Usage

``` r
get_module_eigengene(
  otu_mat,
  graph_obj,
  module_col = "Modularity",
  exclude_others = TRUE,
  scale_pca = TRUE
)
```

## Arguments

- otu_mat:

  Numeric matrix. Rows = OTUs/ASVs, columns = samples. Must have
  rownames (OTU IDs) and colnames (sample IDs).

- graph_obj:

  A tbl_graph from build_graph_from_mat or build_graph_from_df. Must
  have node attribute `name` and a module column.

- module_col:

  Character. Name of the module column in graph_obj nodes. One of
  `"Modularity"`, `"modularity3"`, `"modularity2"`.

- exclude_others:

  Logical. If TRUE, exclude "Others" module from output.

- scale_pca:

  Logical. If TRUE, scale variables before PCA (recommended).

## Value

A data frame with samples as rows and module eigengenes as columns.
Column names are module names.
