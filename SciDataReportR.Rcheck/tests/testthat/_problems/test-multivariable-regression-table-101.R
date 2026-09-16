# Extracted from test-multivariable-regression-table.R:101

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("glmnet")
set.seed(306)
df <- data.frame(
    y = rnorm(70),
    x1 = rnorm(70),
    x2 = rnorm(70),
    cov = rnorm(70)
  )
res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "y",
    PredictorVars = c("x1", "x2"),
    Covars = "cov",
    Method = "elasticnet",
    CVFolds = 5
  )
