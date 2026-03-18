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
  plot_margin = ggplot2::margin(10, 10, 10, 10),
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

  Title horizontal justification (0–1). Default 0.5 (center).

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

  Outer plot margin (ggplot2::margin). Default margin(10,10,10,10).

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
