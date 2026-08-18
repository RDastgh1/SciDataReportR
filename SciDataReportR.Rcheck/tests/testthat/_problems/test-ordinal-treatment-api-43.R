# Extracted from test-ordinal-treatment-api.R:43

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(1)
df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 12)),
    outcome = stats::rnorm(24),
    severity = ordered(rep(c("None", "Mild", "Moderate"), 8),
                       levels = c("None", "Mild", "Moderate")),
    marker = stats::rnorm(24)
  )
continuous <- PlotInteractionEffectsMatrix(
    df_Test, interVar = "group", outcome_vars = "outcome",
    predictor_vars = c("severity", "marker"), TreatOrdinalAs = "Continuous"
  )
