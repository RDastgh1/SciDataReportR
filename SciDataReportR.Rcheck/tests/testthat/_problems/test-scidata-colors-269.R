# Extracted from test-scidata-colors.R:269

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(SampleData, envir = environment())
data(SampleVariableTypes, envir = environment())
labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
variables <- getNumVars(labelled)[1:25]
categories <- paste0("Domain", seq_along(variables))
plot_obj <- suppressWarnings(CreateZScorePlot(
    data = labelled,
    TargetVar = "Diagnosis",
    variables = variables,
    VariableCategories = categories
  ))
built <- suppressWarnings(ggplot2::ggplot_build(plot_obj))
