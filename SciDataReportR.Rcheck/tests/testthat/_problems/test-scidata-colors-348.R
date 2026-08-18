# Extracted from test-scidata-colors.R:348

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(SimulatedPhenotypeData, envir = environment())
training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
