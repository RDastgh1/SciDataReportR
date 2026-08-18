# Extracted from test-multivariable-regression-table.R:467

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(643)
df <- data.frame(
    cohort = factor(sample(c("A", "B", "C"), 120, replace = TRUE)),
    y = rnorm(120),
    x1 = rnorm(120)
  )
ovr <- MultivariableRegressionTable(
    df, c("cohort", "y"), "x1",
    outcome_modes = c(cohort = "one_vs_rest")
  )
