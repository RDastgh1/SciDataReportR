# Extracted from test-label-ordinal-contract.R:78

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("gtsummary")
df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 4)),
    severity = ordered(
      c("Mild", "Moderate", "Severe", "Moderate", "Mild", "Severe", "Moderate", "Mild"),
      levels = c("Mild", "Moderate", "Severe")
    )
  )
sjlabelled::set_label(df_Test$severity) <- "Clinical severity"
tbl <- MakeTable1(df_Test, variables = "severity", TreatOrdinalAs = "Both")
labels <- tbl$table_body$label[tbl$table_body$row_type == "label"]
expect_true(any(labels == "Clinical severity (categorical)"))
