# Update module assignments from an existing column in graph object

Uses the values of an existing node column as the new module assignment,
updating `Modularity`, `modularity2`, and `modularity3` for seamless
integration with
[`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md).
All other node columns are preserved.

## Usage

``` r
update_graph_modules2(graph_obj, modules_new, levels = NULL)
```

## Arguments

- graph_obj:

  A graph object returned by `build_graph_*()`.

- modules_new:

  Character. Name of an existing column in `graph_obj` nodes whose
  values define the new module assignment. The column may or may not
  contain `"Others"`.

- levels:

  Character vector. Optional explicit module order. If `NULL`, order by
  module size (descending) and place `"Others"` at the end when present.

## Value

A new graph object with updated `Modularity`, `modularity2`, and
`modularity3` (and `modularity` if that column exists). All other node
columns are preserved.

## Examples

``` r
NULL
#> NULL
```
