# Extracted from test-semantic-colors.R:127

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
set.seed(9)
df <- data.frame(v = c(rnorm(60, 10, 2), 40, -15))
result <- IQROutliers(df, Variable = "v")
plot_obj <- result$p
skip_if(is.null(plot_obj), "IQROutliers returned no plot")
built <- ggplot2::ggplot_build(plot_obj)
colours <- unlist(lapply(built$data, function(layer) layer$colour))
colours <- unique(colours[!is.na(colours)])
expect_true("red" %in% colours || "black" %in% colours)
expect_false(any(colours %in% SciDataPalette(names = FALSE)))
