# Extracted from test-multivariable-regression-table.R:135

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("glmnet")
set.seed(512)
n <- 80
df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    cov = rnorm(n),
    cov_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
df$y <- 0.8 * df$x1 + 0.5 * df$cov + rnorm(n)
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    covariates = c("cov", "cov_cat"),
    Method = "lasso",
    CVFolds = 5
  )
tuning <- res$Metadata$AnalysisSettings$Tuning$y
expect_true(all(c("PenaltyFactors", "UnpenalizedTerms") %in% names(tuning)))
expect_setequal(tuning$UnpenalizedTerms, c("cov", "cov_cat"))
