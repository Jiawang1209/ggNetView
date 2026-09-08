## Test environments

* local macOS 15 (aarch64-apple-darwin23), R 4.6.1 -- `R CMD check --as-cran`
  via `devtools::check()`, 2026-09-08: **0 errors | 0 warnings | 0 notes**
  (examples, `--run-donttest`, tests and vignette re-building all OK).
* win-builder (devel and release) -- to be run immediately before submission.
* R-hub (ubuntu-latest, windows-latest, macos-latest) -- to be run immediately
  before submission.

## R CMD check results

Local `--as-cran` run on the version in this repository (0.2.1):

    0 errors | 0 warnings | 0 notes

## Submission notes

This is the first CRAN submission of `ggNetView`.

`ggNetView` provides a unified, reproducible framework for analyzing and
visualizing complex biological, ecological, and microbial association
networks, with deterministic layout generators built on `ggraph` and
`ggplot2`.

* Dependencies are all on CRAN. `WGCNA` is a hard dependency (`Imports`) and is
  itself a CRAN package. `SpiecEasi` and `SparCC` appear in the Description and
  in the `method` arguments as *algorithm* names: both are implemented inside
  this package in C++ (`src/`, via `Rcpp` / `RcppArmadillo`), so neither
  requires an external Bioconductor or GitHub package.
* `LinkingTo: Rcpp, RcppArmadillo`; compiled code builds without warnings under
  the check's compilation flags.
* Exported functions are documented with runnable `\examples` using the bundled
  datasets, with two exceptions: `angle_utils`, a concept page describing the
  degree/radian convention rather than a function, and `mantel_utils`, the
  shared page for `mantel_pairwise()`, `mantel_between_blocks()` and
  `mantel_block_vs_col()`, which currently has none.
* `R/*.R` contains no non-ASCII characters.
* `LazyData: true`; the bundled datasets are small and stored compressed.
* The random-number generator is seeded locally: functions taking a `seed`
  argument restore the caller's RNG kind and state on exit, so the package does
  not change the state of the user's session.

## Downstream dependencies

There are currently no downstream dependencies on CRAN.
