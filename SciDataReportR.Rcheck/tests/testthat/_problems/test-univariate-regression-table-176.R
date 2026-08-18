# Extracted from test-univariate-regression-table.R:176

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(918)
df <- data.frame(
    y = rnorm(60),
    x1 = rnorm(60)
  )
lifecycle::expect_deprecated(
    res_old <- UnivariateRegressionTable(
      data = df,
      outcome_vars = "y",
      predictor_vars = "x1",
      ReturnModels = TRUE
    )
  )
