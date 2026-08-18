# Extracted from test-ordinal-treatment-api.R:12

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 6)),
    severity = ordered(rep(c("None", "Mild", "Moderate"), 4),
                       levels = c("None", "Mild", "Moderate")),
    marker = seq_len(12)
  )
continuous <- PlotZScore(
    df_Test, TargetVar = "group", variables = c("severity", "marker"),
    TreatOrdinalAs = "Continuous"
  )
