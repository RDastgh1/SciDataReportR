# Extracted from test-multivariable-regression-table.R:222

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(88)
n <- 60
df <- data.frame(y = rnorm(n), a = rnorm(n))
df$b <- df$a
df$c <- 2 * df$a
res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("a", "b", "c"),
    Method = "lm"
  )
expect_gt(res$ModelSummary$AliasedTermCount, 0)
