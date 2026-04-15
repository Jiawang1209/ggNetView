# Generate a named color palette for discrete classes

Generate a named color palette for discrete classes

## Usage

``` r
get_palette(classes, others_label = "Others")
```

## Arguments

- classes:

  Character string. The discrete class names or factor levels to map to
  colors.

- others_label:

  Character, (default = "Others").

## Value

A named character vector whose names are the discrete classes and whose
values are hex colors. If \`others_label\` appears in \`classes\`, it is
always mapped to a neutral grey so "Others" has a stable color.

## Examples

``` r
get_palette(c("Module1", "Module2", "Module3"))
#>   Module1   Module2   Module3 
#> "#8dd3c7" "#ffffb3" "#bebada" 
get_palette(c("Module1", "Module2", "Others"))
#>   Module1   Module2    Others 
#> "#8dd3c7" "#ffffb3" "#bdbdbd" 
```
