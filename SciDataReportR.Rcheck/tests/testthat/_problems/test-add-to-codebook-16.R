# Extracted from test-add-to-codebook.R:16

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_codebook <- data.frame(
    Variable = "sex",
    Label = "Sex assigned at birth",
    Type = "Double",
    Domain = "Clinical",
    stringsAsFactors = FALSE
  )
actual <- AddToCodebook(
    df_codebook,
    "age",
    "Age at enrollment",
    "Double",
    Domain = "Clinical"
  )
