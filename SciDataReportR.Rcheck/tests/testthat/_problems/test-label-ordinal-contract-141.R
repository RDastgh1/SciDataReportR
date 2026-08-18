# Extracted from test-label-ordinal-contract.R:141

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 8)),
    severity = ordered(
      rep(c("None", "Mild", "Moderate", "Severe"), 4),
      levels = c("None", "Mild", "Moderate", "Severe")
    ),
    marker = seq_len(16)
  )
anova_out <- PlotAnovaRelationshipsMatrix(
    df_Test,
    CatVars = c("group", "severity"),
    ContVars = "marker",
    Relabel = FALSE
  )
