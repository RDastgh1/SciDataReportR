# Extracted from test-label-ordinal-contract.R:14

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(severity = c(1, 2, 3, NA_real_))
df_Codebook <- data.frame(
    Variable = "severity",
    Recode = "Yes",
    Code = "1=Mild;2=Moderate;3=Severe",
    Type = "Ordinal",
    Label = "Clinical severity"
  )
out <- RevalueData(df_Test, df_Codebook)$RevaluedData
expect_true(is.ordered(out$severity))
expect_identical(attr(out$severity, "scidr_type"), "ordinal")
