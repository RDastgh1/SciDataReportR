# Extracted from test-univariate-regression-table.R:246

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(911)
x1 <- rnorm(80)
probability <- stats::plogis(-0.3 + 0.8 * x1)
df <- data.frame(
    ybin = factor(
      ifelse(stats::runif(80) < probability, "Case", "Control"),
      levels = c("Control", "Case")
    ),
    x1 = x1
  )
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "ybin",
    predictor_vars = "x1",
    Standardize = TRUE,
    ReturnModels = TRUE
  )
