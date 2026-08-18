# Extracted from test-scidata-colors.R:51

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_identical(.SciDataColorValues(50), .SciDataColorValues(50))
