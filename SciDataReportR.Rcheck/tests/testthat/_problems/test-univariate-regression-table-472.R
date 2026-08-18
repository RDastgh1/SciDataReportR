# Extracted from test-univariate-regression-table.R:472

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(922)
df <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )
deprecation_warnings <- testthat::capture_warnings(
    res <- MakeUnivariateRegressionTable(
      Data = df,
      OutcomeVars = "y",
      PredictorVars = "x1",
      ReturnModels = TRUE
    )
  )
