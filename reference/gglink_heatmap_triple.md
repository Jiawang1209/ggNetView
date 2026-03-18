# Visualize multi-orientation environmental–species correlation heatmaps2

Visualize multi-orientation environmental–species correlation heatmaps2

## Usage

``` r
gglink_heatmap_triple(
  Environment,
  Experiment,
  edge,
  node,
  sample_col = "Sample",
  delim = ",",
  hub_n = NULL,
  r = 6
)
```

## Arguments

- Environment:

  character or data.frame File path or data frame of environment data.

- Experiment:

  character or data.frame File path or data frame of experiment data.

- edge:

  character or data.frame File path or data frame of edge data.

- node:

  character or data.frame File path or data frame of node data.

- sample_col:

  Character (default = "Sample") Column name used as sample ID when
  input is a data frame or file.

- delim:

  Character (default = ",") Delimiter for reading input files.

- hub_n:

  Integer (default = NULL) Number of hub nodes used in layout; if NULL,
  uses all nodes.

- r:

  numeric (default = 6)

## Value

a ggplot2 object

## Examples

``` r
NULL
#> NULL
```
