# Extracted from test-scidata-colors.R:144

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df <- data.frame(g = factor(c("a", "b", "c")), y = 1:3)
built_fill <- ggplot2::ggplot_build(
    ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, fill = g)) +
      ggplot2::geom_col() +
      scale_fill_SciData()
  )
