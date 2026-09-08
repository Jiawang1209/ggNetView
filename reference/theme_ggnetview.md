# ggNetView Theme

ggNetView Theme

## Usage

``` r
theme_ggnetview(
  base_size = 12,
  base_family = NULL,
  title_face = "bold",
  title_size = 18,
  title_hjust = 0.5,
  subtitle_size = 12,
  caption_size = 9,
  background = "white",
  foreground = NULL,
  border = FALSE,
  plot_margin = ggplot2::margin(20, 20, 20, 20),
  grid = c("none", "x", "y", "both")
)
```

## Arguments

- base_size:

  Base text size. Default 12.

- base_family:

  Base font family. Default NULL.

- title_face:

  Title font face (e.g., "bold"). Default "bold".

- title_size:

  Title font size. Default 18.

- title_hjust:

  Title horizontal justification (0-1). Default 0.5 (center).

- subtitle_size:

  Subtitle size. Default 12.

- caption_size:

  Caption size. Default 9.

- background:

  Plot background color; NULL for transparent. Default "white".

- foreground:

  Foreground color for strip background/border; NULL to skip.

- border:

  Logical; draw panel border with \`foreground\` color if TRUE. Default
  FALSE.

- plot_margin:

  Outer plot margin (ggplot2::margin). Default `margin(20, 20, 20, 20)`
  (in pt). Bumped from `10` to give labels rendered outside the plot
  panel (`coord_equal(clip = "off")` in
  [`ggNetView()`](https://jiawang1209.github.io/ggNetView/reference/ggNetView.md))
  enough room to render without being cropped by the device boundary.

- grid:

  "none","x","y","both" to toggle major grid lines. Default "none".

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)
library(ggNetView)
ggplot(mtcars, aes(wt, mpg, color = hp)) +
geom_point(size = 2) +
ggtitle("Example of ggNetView Theme") +
theme_ggnetview()

```
