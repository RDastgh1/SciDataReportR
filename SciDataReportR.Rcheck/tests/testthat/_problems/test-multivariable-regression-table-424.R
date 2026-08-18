# Extracted from test-multivariable-regression-table.R:424

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(641)
df <- data.frame(
    severity = ordered(
      sample(c("Mild", "Moderate", "Severe"), 90, replace = TRUE),
      levels = c("Mild", "Moderate", "Severe")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )
res <- MultivariableRegressionTable(df, "severity", c("x1", "x2"), Method = "lm")
