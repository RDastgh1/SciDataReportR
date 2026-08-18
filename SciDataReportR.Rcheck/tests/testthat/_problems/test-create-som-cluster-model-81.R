# Extracted from test-create-som-cluster-model.R:81

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("aweSOM")
skip_if_not_installed("kohonen")
skip_if_not_installed("mclust")
skip_if_not_installed("R.utils")
suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
data("SampleData", package = "SciDataReportR")
data("SampleVariableTypes", package = "SciDataReportR")
df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory",
    k_range = 2,
    models = 1,
    lpa_timeout_seconds = 30
  )
