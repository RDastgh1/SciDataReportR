# Extracted from test-univariate-regression-table.R:204

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(910)
df <- data.frame(
    ychar = sample(c("No", "Yes"), 70, replace = TRUE),
    x1 = rnorm(70)
  )
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "ychar",
    predictor_vars = "x1",
    ReturnModels = TRUE
  )
