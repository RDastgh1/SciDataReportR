# Extracted from test-scidata-colors.R:102

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
labels <- .SciDataContrastText(SciDataPalette(names = FALSE))
