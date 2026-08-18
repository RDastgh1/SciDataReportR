# Extracted from test-scidata-colors.R:65

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
values <- .SciDataNamedValues(
    c("1", "2", "3", "Noise"),
    fixed = c(Noise = "grey60")
  )
