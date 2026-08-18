# Extracted from test-univariate-regression-table.R:105

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(908)
df <- data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = "x1"
  )
expect_null(res$ModelSummaries)
expect_false(res$Metadata$AnalysisSettings$ReturnModels)
