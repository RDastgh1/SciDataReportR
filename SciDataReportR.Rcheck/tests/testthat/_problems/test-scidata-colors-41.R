# Extracted from test-scidata-colors.R:41

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
for (n in c(34, 35, 68, 102, 136, 170)) {
    expect_equal(anyDuplicated(.SciDataColorValues(n)), 0, info = paste("n =", n))
  }
