test_that("PlotClusterPrevalence summarizes binary features by cluster", {
  df <- data.frame(
    Cluster = factor(rep(c("B", "A"), each = 4), levels = c("B", "A")),
    high_bp = c("Yes", "Yes", "No", NA, "No", "No", "No", "Yes"),
    anemia = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
    ckd = c(1, 0, 0, 0, 1, 1, 1, 0)
  )
  attr(df$high_bp, "label") <- "Systolic BP >= 130 mmHg"

  plot <- PlotClusterPrevalence(df, "Cluster", c("high_bp", "anemia", "ckd"))
  summary <- attr(plot, "prevalence")

  expect_s3_class(plot, "ggplot")
  expect_identical(levels(summary$Cluster), c("All", "B", "A"))
  bp_b <- summary[summary$Variable == "high_bp" & summary$Cluster == "B", ]
  expect_equal(bp_b$Prevalence, 2 / 3)
  expect_equal(bp_b$n, 3)
  expect_equal(summary$Prevalence[summary$Variable == "ckd" & summary$Cluster == "All"], 0.5)
  expect_true("Systolic BP >= 130 mmHg" %in% as.character(summary$Label))
  expect_no_error(ggplot2::ggplot_build(plot))

  no_cohort <- attr(PlotClusterPrevalence(df, "Cluster", "anemia", show_cohort = FALSE), "prevalence")
  expect_false("All" %in% no_cohort$Cluster)
  expect_error(PlotClusterPrevalence(df, "Cluster", "missing_var"), "not found")
})
