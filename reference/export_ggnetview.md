# Save a ggNetView with sensible defaults

Save a ggNetView with sensible defaults

## Usage

``` r
export_ggnetview(p, filename, height, width, limitsize = F, dpi = 300)
```

## Arguments

- p:

  Plot to save, defaults to last plot displayed.

- filename:

  File name to create on disk.

- height:

  Plot size in units expressed by the units argument. If not supplied,
  uses the size of the current graphics device.

- width:

  Plot size in units expressed by the units argument. If not supplied,
  uses the size of the current graphics device.

- limitsize:

  When TRUE (the default), export_ggnetview() will not save images
  larger than 50x50 inches, to prevent the common error of specifying
  dimensions in pixels.

- dpi:

  Plot resolution. Also accepts a string input: "retina" (320), "print"
  (300), or "screen" (72). Only applies when converting pixel units, as
  is typical for raster output types.

## Value

a graph output

## Examples

``` r
NULL
#> NULL
```
