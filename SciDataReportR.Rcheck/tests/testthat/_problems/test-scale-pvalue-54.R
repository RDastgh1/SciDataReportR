# Extracted from test-scale-pvalue.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(
    scale_color_pvalue(palette = "magma"),
    "should be one of"
  )
