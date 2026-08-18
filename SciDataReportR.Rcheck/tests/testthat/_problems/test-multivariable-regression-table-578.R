# Extracted from test-multivariable-regression-table.R:578

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(645)
n <- 120
df <- data.frame(
    check.names = FALSE,
    "Gender at birth" = factor(sample(c("Female", "Male", "Other"), n, replace = TRUE)),
    "Education level" = factor(sample(c("High school", "College", "Graduate"), n, replace = TRUE)),
    "Continuous outcome" = rnorm(n),
    "Callosum Anterior Frontal_sfr" = rnorm(n),
    "Left Inferior Fronto-occipital_sfr" = rnorm(n),
    "Left Inferior Longitudinal" = rnorm(n),
    "Left Corticospinal_sfr" = rnorm(n),
    "Age" = rnorm(n, 50, 10)
  )
ordinary <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "Continuous outcome",
    predictor_vars = c(
      "Left Inferior Fronto-occipital_sfr",
      "Left Inferior Longitudinal",
      "Left Corticospinal_sfr"
    ),
    covariates = "Age",
    Method = "lm"
  )
