# Extracted from test-multivariable-regression-table.R:501

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(
    cohort = factor(rep(c("A", "B", "C"), each = 15)),
    x1 = rnorm(45)
  )
expect_error(
    MultivariableRegressionTable(df, "cohort", "x1", reference_levels = c(cohort = "D")),
    "not retained"
  )
