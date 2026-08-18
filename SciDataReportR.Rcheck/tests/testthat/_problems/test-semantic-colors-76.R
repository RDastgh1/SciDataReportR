# Extracted from test-semantic-colors.R:76

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
set.seed(4)
old_data <- data.frame(id = 1:40, age = rnorm(40, 60, 8), grp = rep(c("A", "B"), 20))
new_data <- old_data
new_data$age[1:5] <- new_data$age[1:5] + 1
new_data <- new_data[-c(3, 7), ]
new_data$extra <- rnorm(nrow(new_data))
comparison <- CompareDatasets(old_data, new_data, keys = "id")
plots <- PlotDatasetComparison(comparison, interactive = FALSE)
fills <- unique(ggplot2::ggplot_build(plots$Checks)$data[[1]]$fill)
expect_true(all(fills %in% StatusPalette()))
expect_false(any(fills %in% SciDataPalette(names = FALSE)))
