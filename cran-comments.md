# cran-comments

## Submission: SciDataReportR 20.10.0

This is a new submission (first release of this package on CRAN).

## Test environments

* local: macOS (Apple Silicon), R 4.5.3, and R-devel (2026-04-24 build) via Docker
* win-builder: R-devel and R-release (checked 2026-07-05 / 07-06; issues found and fixed, see below)
* R-hub: Windows / Ubuntu / macOS      <!-- update after running -->

## Round-trip notes (for CRAN reviewer context; safe to ignore)

An initial win-builder R-devel check surfaced 2 ERRORs / 3 NOTEs. Root causes
and fixes:

* `MakeUnivariateRegressionTable()` calls `gtsummary::tbl_regression()`, which
  requires `broom.helpers` (>= 1.20.0) for lm/glm models. `broom.helpers` was
  only a transitive Suggests (via gtsummary), so it was never installed on a
  fresh check machine. Added `broom.helpers` to Imports with a
  `requireNamespace()` guard.
* `MultivariableRegressionTable()` checked for the optional `pROC` package
  before validating outcome factor levels post missing-data handling; reordered
  so data validation happens first.
* Fixed an internal `stop(paste(..., e$message))` that silently dropped the
  real error text when `e$message` was `NULL`/`character(0)` (a common
  rlang/cli condition shape) — now uses `conditionMessage(e)`. This is what
  let the two bugs above hide behind an uninformative "Error processing X and
  Y :" message.
* Fixed a `<prefix>` placeholder in `ApplyNormativeTScores.Rd` that rendered as
  invalid raw HTML; fixed `CodebookMergeApp.Rd` figure width to pixels; fixed
  an invalid file URI in README (linked to `LICENSE.md`, which is
  `.Rbuildignore`'d).

All 236 package tests pass under R-devel (verified via a local Docker
`rocker/r-devel` container) and under R 4.5.3 release.

A subsequent win-builder round (R-devel and R-release) then failed only in
vignette re-building: the `PlotCorrelationsHeatmap` vignette rendered
`ggstatsplot` scatterplots that require graphics fonts unavailable on the
Windows check machines ("unable to set or substitute a suitable font"). Those
illustrative chunks are now shown but not evaluated at build time; the code
remains visible for users to run interactively.

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE

* NEW submission NOTE: this is the package's first CRAN release.

* The installed package size (~7 Mb) is dominated by the rendered vignette
  documentation (doc/ ~5 Mb). The package ships ten worked vignettes with
  figure-heavy statistical visualizations; figure resolution has already
  been reduced to keep the source tarball at ~4 Mb.

## Comments

* Version number (20.7.0): the package has been developed and versioned
  publicly on GitHub for several years prior to this first CRAN submission;
  the version continues that established sequence rather than restarting at
  0.1.0, so existing users' upgrade paths remain intact.

* Vignettes use the Quarto engine (`VignetteBuilder: quarto`) with the
  `quarto` package in Suggests. Quarto-based reporting is the core purpose
  of the package. R Markdown vignettes are also included and build with
  knitr/rmarkdown.

* A small number of examples remain in \dontrun{}:
  - `CodebookMergeApp()` / `MergeCodebooks()` launch or depend on an
    interactive Shiny application.
  - The SOM clustering pipeline examples (`Pipeline_SOMClust()`,
    `Project_SOMClust()`) involve long-running model training.

* Suggested packages are used conditionally via `requireNamespace()` guards
  with informative error messages.
