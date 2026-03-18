# Update module names in an existing graph object

Update module names in an existing graph object

## Usage

``` r
update_graph_modules(
  graph_obj,
  modules,
  old_col = NULL,
  new_col = NULL,
  levels = NULL,
  allow_partial = TRUE
)
```

## Arguments

- graph_obj:

  A graph object returned by `build_graph_*()`.

- modules:

  A named vector or a data frame used to rename modules. The mapping is
  `old_module -> new_module`.

- old_col:

  Character. Old-module column name in `modules` when `modules` is a
  data frame. Default `NULL` means auto-detect.

- new_col:

  Character. New-module column name in `modules` when `modules` is a
  data frame. Default `NULL` means auto-detect.

- levels:

  Character vector. Optional explicit module order. If `NULL`, order by
  module size (descending) and place `"Others"` at the end.

- allow_partial:

  Logical (default = `TRUE`). If `TRUE`, modules without rename rules
  keep original names. If `FALSE`, all existing modules must be covered
  by `modules`.

## Value

A new graph object with updated module-related node attributes:
`Modularity`, `modularity2`, and `modularity3` (and `modularity` if that
column exists).

## Examples

``` r
NULL
#> NULL
```
