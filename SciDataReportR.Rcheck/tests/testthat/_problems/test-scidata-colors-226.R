# Extracted from test-scidata-colors.R:226

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(SampleData, envir = environment())
data(SampleVariableTypes, envir = environment())
labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
plot_obj <- PlotCategoricalDistributions(
    labelled, variables = c("Diagnosis", "sex", "Genotype")
  )
plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot
built <- ggplot2::ggplot_build(plot_obj)
text_layer <- built$data[[2]]
expect_setequal(unique(text_layer$colour), c("black", "white"))
