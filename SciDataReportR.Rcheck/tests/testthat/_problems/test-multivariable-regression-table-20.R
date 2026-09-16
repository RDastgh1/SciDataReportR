# Extracted from test-multivariable-regression-table.R:20

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("pROC")
set.seed(104)
df <- data.frame(
    y = rnorm(80),
    ybin = factor(sample(c("No", "Yes"), 80, replace = TRUE)),
    x1 = rnorm(80),
    x2 = rnorm(80),
    cov = rnorm(80)
  )
attr(df$x1, "label") <- "Marker one"
res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = c("y", "ybin"),
    PredictorVars = c("x1", "x2"),
    Covars = "cov",
    Method = "lm"
  )
