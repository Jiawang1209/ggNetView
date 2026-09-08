# Auto-detect angle unit and normalise to radians

Internal helper. Accepts a single numeric and returns the equivalent
angle in radians. The rule is:

- `|x| <= 2*pi` -\> treated as radians (e.g. `pi/2`, `pi`, `2*pi`).

- `|x| > 2*pi` -\> treated as degrees and converted (e.g. `45`, `90`,
  `180`).

If the value falls in the ambiguous "small-but-bigger-than-pi" zone
(i.e. `(pi, 2*pi]`), a
[`message()`](https://rdrr.io/r/base/message.html) is emitted reporting
the interpretation and pointing the user at
[`deg()`](https://jiawang1209.github.io/ggNetView/reference/deg.md) as a
way to be explicit.

## Usage

``` r
.normalize_angle(x, name = "angle")
```

## Arguments

- x:

  Single finite numeric.

- name:

  Parameter name used in error / message text.

## Value

Numeric. Radians.
