# Extracted from test-clustering-pipelines.R:338

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("kohonen")
skip_if_not_installed("aweSOM")
skip_if_not_installed("tidyLPA")
skip_if_not_installed("dbscan")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training") %>%
    dplyr::slice(1:40)
