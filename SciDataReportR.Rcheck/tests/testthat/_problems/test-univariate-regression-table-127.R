# Extracted from test-univariate-regression-table.R:127

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(917)
df <- data.frame(
    y = rnorm(90),
    ybin = factor(
      sample(c("Control", "Case"), 90, replace = TRUE),
      levels = c("Control", "Case")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )
attr(df$x1, "label") <- "Marker one"
attr(df$y, "label") <- "Outcome one"
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("y", "ybin"),
    predictor_vars = c("x1", "x2"),
    ReturnModels = TRUE
  )
