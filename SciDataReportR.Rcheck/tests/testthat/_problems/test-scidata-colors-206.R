# Extracted from test-scidata-colors.R:206

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(SampleData, envir = environment())
data(SampleVariableTypes, envir = environment())
labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
plot_obj <- PlotCategoricalDistributions(
    labelled, variables = c("Diagnosis", "sex")
  )
plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot
fills <- unique(ggplot2::ggplot_build(plot_obj)$data[[1]]$fill)
expect_true(all(fills %in% .SciDataColorValues(170)))
