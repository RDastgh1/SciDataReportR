# Extracted from test-multivariable-regression-table.R:447

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("glmnet")
set.seed(642)
df <- data.frame(
    cohort = factor(sample(c("A", "B", "C"), 120, replace = TRUE)),
    x1 = rnorm(120),
    x2 = rnorm(120),
    age = rnorm(120, 50, 10)
  )
res <- MultivariableRegressionTable(
    df, "cohort", c("x1", "x2"), covariates = "age",
    Method = "lasso", CVFolds = 3, reference_levels = c(cohort = "B")
  )
