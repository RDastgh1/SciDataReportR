# Extracted from test-label-ordinal-contract.R:42

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    group = factor(c("A", "B", "A")),
    severity = ordered(c("Mild", "Moderate", "Severe"),
                       levels = c("Mild", "Moderate", "Severe"))
  )
sjlabelled::set_label(df_Test$severity) <- "Clinical severity"
categorical <- ConvertOrdinalToNumeric(
    df_Test, TreatOrdinalAs = "Categorical", ReturnMetadata = TRUE
  )
