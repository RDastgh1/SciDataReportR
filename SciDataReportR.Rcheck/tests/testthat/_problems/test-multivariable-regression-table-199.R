# Extracted from test-multivariable-regression-table.R:199

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(70)
df <- data.frame(
    grp = factor(c(rep("No", 20), rep("Yes", 20))),
    x1 = c(rep(-10, 20), rep(10, 20)),
    x2 = rnorm(40)
  )
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "grp",
    predictor_vars = c("x1", "x2"),
    Method = "lm"
  )
expect_true(res$ModelSummary$SeparationDetected)
