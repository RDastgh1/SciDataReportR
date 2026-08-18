# Extracted from test-semantic-colors.R:141

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
StatusPalette <- function() {
  c(PASS = "#2E7D32", WARNING = "#F9A825", FAIL = "#C62828")
}
MakeMergeValidation <- function() {
  set.seed(1)
  left <- data.frame(id = 1:50, x = rnorm(50))
  right <- data.frame(id = c(1:45, 100:104), y = rnorm(50))
  merged <- merge(left, right, by = "id")
  ValidateMerge(left, right, merged, keys = "id")
}

# test -------------------------------------------------------------------------
data(SampleData, envir = environment())
plot_obj <- PlotMissingData(SampleData[, c("age", "AXL", "tau")])
plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot
skip_if(is.null(plot_obj), "PlotMissingData returned no plot")
fills <- unique(ggplot2::ggplot_build(plot_obj)$data[[1]]$fill)
expect_false(any(fills %in% SciDataPalette(names = FALSE)))
