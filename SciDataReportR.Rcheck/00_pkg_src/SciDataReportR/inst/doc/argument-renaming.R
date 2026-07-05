## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# make deprecation warnings fire deterministically for this document
options(lifecycle_verbosity = "warning")

## ----setup, message = FALSE---------------------------------------------------
library(SciDataReportR)

## ----warning = TRUE-----------------------------------------------------------
old_style <- MakeDataDictionary(DataFrame = mtcars, numdecimals = 2)

## -----------------------------------------------------------------------------
new_style <- MakeDataDictionary(mtcars, digits = 2)
identical(old_style, new_style)

## -----------------------------------------------------------------------------
identical(MakeDataDictionary(mtcars, 2), new_style)

## -----------------------------------------------------------------------------
# a directional example: xVars/yVars became predictor_vars/outcome_vars
res <- PlotCorrelationsHeatmap(
  mtcars,
  predictor_vars = c("hp", "wt"),
  outcome_vars   = "mpg"
)
round(res$Unadjusted$p, 4)

## ----eval = FALSE-------------------------------------------------------------
# options(lifecycle_verbosity = "error")

