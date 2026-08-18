# Extracted from test-univariate-regression-table.R:84

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(924)
df <- data.frame(
    "Global Mean T" = rnorm(80),
    "Treatment Arm/Group" = factor(
      rep(c("Placebo", "Active"), each = 40),
      levels = c("Placebo", "Active")
    ),
    check.names = FALSE
  )
attr(df[["Treatment Arm/Group"]], "label") <- "Treatment group"
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = "Global Mean T",
    predictor_vars = "Treatment Arm/Group"
  )
