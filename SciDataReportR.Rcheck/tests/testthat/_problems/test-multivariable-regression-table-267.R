# Extracted from test-multivariable-regression-table.R:267

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("glmnet")
set.seed(202)
n <- 90
df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
df$y <- 1.0 * df$x1 + rnorm(n)
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    Method = "lasso",
    CVFolds = 5
  )
expect_true(is.na(res$ModelSummary$ModelPValue))
ann <- SciDataReportR:::ScidrColumnAnnotations(res$ModelSummary)
