# Convert degrees to radians

Tiny helper, used to disambiguate when the magnitude rule below would
otherwise misclassify a small degree value as radians (e.g. if you
actually want `5` degrees, write `deg(5)`).

## Usage

``` r
deg(d)
```

## Arguments

- d:

  Numeric. Degrees.

## Value

Numeric. Radians (`d * pi / 180`).

## Examples

``` r
deg(45)   # 0.7853982
#> [1] 0.7853982
deg(180)  # 3.141593
#> [1] 3.141593
```
