# Extracted from test-univariate-regression-table.R:52

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(923)
df <- data.frame(
    "Global Mean T" = rnorm(90),
    "Global Mean T change" = rnorm(90),
    "log AB-38 CSF" = rnorm(90),
    "log IL-12/IL-23p40 Plasma" = rnorm(90),
    "Age-years" = rnorm(90),
    check.names = FALSE
  )
res <- MakeUnivariateRegressionTable(
    data = df,
    outcome_vars = c("Global Mean T", "Global Mean T change"),
    predictor_vars = c("log AB-38 CSF", "log IL-12/IL-23p40 Plasma"),
    covariates = "Age-years",
    ReturnModels = TRUE
  )
