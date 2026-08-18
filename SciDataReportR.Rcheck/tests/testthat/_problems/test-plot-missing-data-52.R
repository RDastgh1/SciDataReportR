# Extracted from test-plot-missing-data.R:52

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(
    Visit = 1:4,
    Site = c("A", "A", NA, NA),
    Marker = c(1, NA, NA, NA)
  )
p <- PlotMissingData(
    df_Test,
    x_var = "Visit",
    facet_by = "Site",
    Relabel = FALSE
  )
