# One-shot, non-interactive RMT threshold scan with auto selection

One-shot, non-interactive RMT threshold scan with auto selection

## Usage

``` r
ggNetView_RMT(
  mat,
  transfrom.method = c("none", "scale", "center", "log2", "log10", "ln", "rrarefy",
    "rrarefy_relative"),
  method = c("WGCNA", "SpiecEasi", "SPARCC", "cor"),
  cor.method = c("pearson", "kendall", "spearman"),
  SpiecEasi.method = c("mb", "glasso"),
  nr_thresholds = 51,
  interval = NULL,
  unfold.method = c("gaussian", "spline"),
  bandwidth = "nrd0",
  nr.fit.points = 51,
  discard.outliers = TRUE,
  discard.zeros = TRUE,
  min.mat.dim = 40,
  max.ev.spacing = 3,
  save_plots = FALSE,
  out_dir = "RMT_plots",
  verbose = TRUE,
  seed = 1115
)
```

## Arguments

- mat:

  Numeric matrix. A numeric matrix with samples in colums and variables
  in rows \#' @param transfrom.method Character. Data transformation
  methods applied before correlation analysis. Options include: "none"
  (raw data), "scale" (z-score standardization), "center" (mean
  centering only), "log2" (log2 transfrom), "log10" (log10 transfrom),
  "ln" (natural transfrom ), "rrarefy" (random rarefaction using
  [`vegan::rrarefy`](https://vegandevs.github.io/vegan/reference/rarefy.html)),
  "rrarefy_relative" (rarefy then convert to relative abundance).

- method:

  Character. Relationship analysis methods. Options include: "WGCNA",
  "SpiecEasi", "SPARCC" and "cor".

- cor.method:

  Character. Correlation analysis method. Options include "pearson",
  "kendall", and "spearman".

- SpiecEasi.method:

  Character. Method used in `SpiecEasi` network inference; options
  include "mb" and "glasso".

- nr_thresholds:

  Integer. Number of thresholds to scan (default 51).

- interval:

  Optional numeric length-2 vector \[min, max\] for scanning range over
  \|mat\| upper triangle (default auto).

- unfold.method:

  "gaussian" or "spline" (default "gaussian").

- bandwidth:

  Bandwidth for gaussian unfolding (passed to stats::density), default
  "nrd0".

- nr.fit.points:

  Integer, number of support points for spline unfolding (default 51).

- discard.outliers:

  Logical. Remove eigenvalue outliers via IQR before unfolding (default
  TRUE).

- discard.zeros:

  Logical. Drop all-zero rows/cols after thresholding (default TRUE).

- min.mat.dim:

  Integer. Early stop if effective dim \< this (default 40).

- max.ev.spacing:

  Numeric. Cutoff for spacing tail when computing NNSD metrics (default
  3).

- save_plots:

  Logical. Save PNGs (non-interactive) into out_dir (default FALSE).

- out_dir:

  Character. Output directory for plots (default "RMT_plots").

- verbose:

  Logical. Print progress (default TRUE).

- seed:

  Integer (default = 1115). Random seed for reproducibility.

## Value

A list with: chosen_threshold, chosen_reason, tested_thresholds, scores
(data.frame), unfolded (last-step unfolding), meta (matrix info &
params), plots (file paths if saved).

## Examples

``` r
NULL
#> NULL


# m <- cor(scale(matrix(rnorm(40000), 200, 200)))
# res <- ggNetView_RMT(m, save_plots=TRUE)
# res$chosen_threshold
```
