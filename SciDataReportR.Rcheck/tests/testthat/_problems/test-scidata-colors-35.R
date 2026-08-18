# Extracted from test-scidata-colors.R:35

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
for (n in c(1, 2, 33, 34, 35, 68, 100, 170)) {
    expect_length(.SciDataColorValues(n), n)
  }
