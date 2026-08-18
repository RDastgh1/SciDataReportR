# Extracted from test-scale-pvalue.R:135

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Guide <- data.frame(
    x = seq_len(7),
    y = seq_len(7),
    PValue = c(1, 0.1, 0.05, 0.01, 0.001, 1e-5, 1e-8)
  )
expected_positions <- c(0, 0.16, 0.32, 0.56, 0.75, 0.90, 1)
expected_labels <- c("1", "0.1", "0.05", "0.01", "0.001", "1e-05", "1e-08")
color_build <- ggplot2::ggplot_build(
    ggplot2::ggplot(
      df_Guide,
      ggplot2::aes(x = x, y = y, color = PValue)
    ) +
      ggplot2::geom_point() +
      scale_color_pvalue()
  )
