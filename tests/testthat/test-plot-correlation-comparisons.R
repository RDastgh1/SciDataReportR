test_that("PlotCorrelationComparisons returns labelled, FDR-corrected static results", {
  set.seed(42)
  df_Test <- data.frame(
    Group = factor(rep(c("Reference", "Comparison"), each = 50),
      levels = c("Reference", "Comparison")),
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )
  df_Test$y[51:100] <- df_Test$x[51:100] + rnorm(50, sd = 0.2)
  df_Test$x <- sjlabelled::set_label(df_Test$x, "Exposure")
  df_Test$y <- sjlabelled::set_label(df_Test$y, "Outcome")

  out <- PlotCorrelationComparisons(
    data = df_Test,
    predictor_vars = c("x", "z"),
    outcome_vars = "y",
    group_var = "Group"
  )

  expect_s3_class(out$Unadjusted$plot, "ggplot")
  expect_s3_class(out$FDRCorrected$plot, "ggplot")
  expect_identical(out$Metadata$ReferenceGroup, "Reference")
  expect_identical(out$Metadata$ComparisonGroup, "Comparison")
  expect_true(all(c("PAdjusted", "InferenceStatus", "CellId") %in% names(out$Results)))
  expect_true(all(out$Results$PredictorLabel %in% c("Exposure", "z")))
  expect_match(out$Metadata$InferenceStatus, "Pearson")
  expect_equal(
    out$Unadjusted$delta_r,
    out$Unadjusted$r_comparison - out$Unadjusted$r_reference
  )
  expect_equal(
    out$FDRCorrected$p,
    ApplyFDRCorrection(out$Unadjusted$p, fdr_scope = "matrix", symmetric = FALSE)
  )
  expect_equal(out$p, out$Unadjusted)
  expect_equal(out$p_fdr, out$FDRCorrected)
})

test_that("PlotCorrelationComparisons reports approximation status", {
  set.seed(1)
  df_Test <- data.frame(
    Group = rep(c("A", "B"), each = 20),
    x = rnorm(40),
    y = rnorm(40),
    age = rnorm(40)
  )

  out_spearman <- PlotCorrelationComparisons(
    df_Test,
    predictor_vars = "x",
    outcome_vars = "y",
    group_var = "Group",
    method = "spearman"
  )
  out_partial <- PlotCorrelationComparisons(
    df_Test,
    predictor_vars = "x",
    outcome_vars = "y",
    group_var = "Group",
    covariates = "age"
  )

  expect_match(out_spearman$Results$InferenceStatus, "Approximate.*Spearman")
  expect_match(out_partial$Results$InferenceStatus, "Approximate.*partial")
})

test_that("PlotCorrelationComparisons validates group selection", {
  df_Test <- data.frame(
    Group = rep(c("A", "B", "C"), each = 5),
    x = seq_len(15),
    y = rev(seq_len(15))
  )

  expect_error(
    PlotCorrelationComparisons(df_Test, "x", "y", "Group"),
    "Specify `comparison_group` and `reference_group`"
  )
  expect_error(
    PlotCorrelationComparisons(
      df_Test, "x", "y", "Group",
      comparison_group = "A", reference_group = "Missing"
    ),
    "was not found"
  )

  expect_error(
    PlotCorrelationComparisons(
      df_Test[seq_len(7), ], "x", "y", "Group",
      comparison_group = "A", reference_group = "B"
    ),
    "contains only"
  )
})

test_that("PlotCorrelationComparisons provides Plotly and ggiraph widgets", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("ggiraph")

  set.seed(5)
  df_Test <- data.frame(
    Group = rep(c("A", "B"), each = 20),
    x = rnorm(40),
    y = rnorm(40)
  )

  out <- PlotCorrelationComparisons(
    df_Test,
    predictor_vars = "x",
    outcome_vars = "y",
    group_var = "Group",
    interactive = "both"
  )

  expect_s3_class(out$Interactive$Plotly$FDRCorrected, "plotly")
  expect_s3_class(out$Interactive$Girafe$FDRCorrected, "girafe")
  expect_no_error(plotly::plotly_build(out$Interactive$Plotly$FDRCorrected))
})
