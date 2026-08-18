# Extracted from test-plot-missing-data.R:23

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    FACS_date = as.Date("2026-01-01") + 0:3,
    Group = c("A", "B", NA, "A"),
    Marker = c(1, NA, 3, 4)
  )
sjlabelled::set_label(df_Test$FACS_date) <- "FACS date"
p_date <- PlotMissingData(df_Test, x_var = "FACS_date", Relabel = FALSE)
