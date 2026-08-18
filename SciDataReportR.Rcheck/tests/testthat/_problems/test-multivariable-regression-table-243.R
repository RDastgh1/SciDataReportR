# Extracted from test-multivariable-regression-table.R:243

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(131)
n <- 90
df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
df$y <- 1.2 * df$x1 + rnorm(n)
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    Method = "lm"
  )
expect_equal(res$ModelSummary$ModelStatType, "F")
expect_true(is.finite(res$ModelSummary$ModelPValue))
