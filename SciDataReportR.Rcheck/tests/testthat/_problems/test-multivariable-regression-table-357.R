# Extracted from test-multivariable-regression-table.R:357

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_three <- data.frame(
    cohort = factor(rep(c("A", "B", "C"), each = 10)),
    x1 = rnorm(30)
  )
res <- MultivariableRegressionTable(df_three, "cohort", "x1")
