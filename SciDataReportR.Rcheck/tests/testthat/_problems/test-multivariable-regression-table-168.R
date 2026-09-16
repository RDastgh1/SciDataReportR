# Extracted from test-multivariable-regression-table.R:168

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("gt")
skip_if_not_installed("pROC")
set.seed(618)
df <- data.frame(
    y = rnorm(80),
    x1 = rnorm(80),
    x2 = rnorm(80),
    cov = rnorm(80)
  )
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    covariates = "cov",
    Method = "lm"
  )
