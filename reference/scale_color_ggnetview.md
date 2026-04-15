# Custom discrete color scale for ggNetView

Custom discrete color scale for ggNetView

## Usage

``` r
scale_color_ggnetview(
  classes,
  ...,
  others_label = "Others",
  na_value = "#e0e0e0",
  drop = FALSE
)
```

## Arguments

- classes:

  Character string. The discrete class names or factor levels to map to
  colors.

- ...:

  Additional arguments passed to \`ggplot2::scale_color_manual()\`.

- others_label:

  Character, (default = "Others").

- na_value:

  Color for missing values. Default \`"#e0e0e0"\`.

- drop:

  Logical, passed to \`ggplot2::scale_color_manual()\`.

## Value

A \`ggplot2\` scale object.

## Examples

``` r
library(ggplot2)
df <- data.frame(
  x = 1:4, y = 1:4,
  group = c("Module1", "Module2", "Module3", "Others")
)
ggplot(df, aes(x, y, color = group)) +
  geom_point(size = 6) +
  scale_color_ggnetview(df$group)

```
