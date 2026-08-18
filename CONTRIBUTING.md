# Development checks

Use the smallest validation loop that covers the change you made.

## Routine development

For code and test changes, run the test suite first:

``` r

devtools::test()
```

For a package check that skips vignette building and vignette
re-building, run:

``` r

devtools::check(vignettes = FALSE, args = "--timings")
```

This is the normal local check. It still builds, installs, exercises
examples, and runs the test suite.

## Release and vignette validation

After changing a vignette, and before publishing a release, run the full
check:

``` r

devtools::check(vignettes = TRUE, args = "--timings")
```

This intentionally builds the package vignettes and asks `R CMD check`
to re-build their outputs. The clustering vignette is computationally
substantial, so this check is expected to take longer than the routine
workflow.

When running from a macOS terminal rather than RStudio, make sure
`rmarkdown` can find Pandoc before the full check. For the Quarto
installation used on Apple silicon Macs, run:

``` r

Sys.setenv(RSTUDIO_PANDOC = "/Applications/quarto/bin/tools/aarch64")
```

Then confirm `rmarkdown::find_pandoc()$version` is not `"0"`. On another
platform, set `RSTUDIO_PANDOC` to the directory that contains the local
Pandoc executable.
