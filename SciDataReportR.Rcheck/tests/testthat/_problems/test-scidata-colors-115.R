# Extracted from test-scidata-colors.R:115

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(g = factor(paste0("L", seq_len(40))), y = seq_len(40))
fill_plot <- ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, fill = g)) +
    ggplot2::geom_col() +
    .SciDataFillScale()
