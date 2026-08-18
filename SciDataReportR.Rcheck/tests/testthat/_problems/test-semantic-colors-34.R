# Extracted from test-semantic-colors.R:34

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
plots <- PlotMergeValidation(
    MakeMergeValidation(), Plot = "All", interactive = FALSE
  )
fills <- unique(ggplot2::ggplot_build(plots$Checks)$data[[1]]$fill)
expect_true(all(fills %in% StatusPalette()))
expect_false(any(fills %in% SciDataPalette(names = FALSE)))
