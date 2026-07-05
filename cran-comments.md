# cran-comments

## Submission: SciDataReportR 20.7.0

This is a new submission (first release of this package on CRAN).

## Test environments

* local: macOS (Apple Silicon), R 4.x
* win-builder: R-devel and R-release  <!-- update after running -->
* R-hub: Windows / Ubuntu / macOS      <!-- update after running -->

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
