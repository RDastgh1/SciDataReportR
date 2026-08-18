# Extracted from test-label-ordinal-contract.R:113

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 8)),
    severity = rep(0:3, 4),
    marker = stats::rnorm(16)
  )
df_Codebook <- data.frame(
    Variable = "severity",
    Recode = "Yes",
    Code = "0=None;1=Mild;2=Moderate;3=Severe",
    Type = "Ordinal",
    Label = "Clinical severity"
  )
df_Revalued <- RevalueData(df_Test, df_Codebook)$RevaluedData
out <- PlotMiningMatrix(
    df_Revalued,
    outcome_vars = c("group", "severity", "marker"),
    TreatOrdinalAs = "Both"
  )
