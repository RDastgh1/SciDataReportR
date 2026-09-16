# Extracted from test-multivariable-regression-table.R:72

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(205)
x1 <- rnorm(120)
x2 <- rnorm(120)
probability <- stats::plogis(-0.2 + 0.9 * x1 - 0.4 * x2)
df <- data.frame(
    ybin = factor(ifelse(stats::runif(120) < probability, "Yes", "No")),
    x1 = x1,
    x2 = x2
  )
res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "ybin",
    PredictorVars = c("x1", "x2"),
    Method = "lm"
  )
