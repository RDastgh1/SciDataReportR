test_that("DiagnosticLikelihoodRatioTable calculates binary diagnostic LRs", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    disease = factor(c("No", "No", "No", "No", "Yes", "Yes", "Yes", "Yes"), levels = c("No", "Yes")),
    test = factor(c("Negative", "Negative", "Negative", "Positive", "Negative", "Positive", "Positive", "Positive"), levels = c("Negative", "Positive"))
  )
  out <- DiagnosticLikelihoodRatioTable(df_Test, "disease", "test")
  expect_named(out, c("FormattedTable", "LargeTable", "BinaryFormattedTable", "BinaryLargeTable", "Results", "BinarySummary", "Metadata"))
  expect_equal(out$BinarySummary$Sensitivity, 3 / 4)
  expect_equal(out$BinarySummary$Specificity, 3 / 4)
  expect_equal(out$BinarySummary$LRPositive, 3)
  expect_equal(out$BinarySummary$LRNegative, 1 / 3)
  expect_true(all(is.finite(out$Results$LRLowerCI)))
})

test_that("DiagnosticLikelihoodRatioTable supports categorical levels and strata", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    disease = factor(rep(c("No", "Yes"), each = 6), levels = c("No", "Yes")),
    genotype = factor(rep(c("e2", "e3", "e4"), 4), levels = c("e2", "e3", "e4")),
    sex = rep(c("Female", "Male"), 6),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(DiagnosticLikelihoodRatioTable(df_Test, "disease", "genotype", stratify_by = "sex"))
  expect_equal(nrow(out$Results), 6)
  expect_true(all(c("sex", "Stratum") %in% names(out$Results)))
  expect_equal(sort(unique(out$Results$ResultLevel)), c("e2", "e3", "e4"))
  expect_equal(nrow(out$BinarySummary), 0)
})

test_that("DiagnosticLikelihoodRatioTable resolves levels and handles zero cells", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    outcome = c(FALSE, FALSE, TRUE, TRUE),
    numeric_test = c(0, 0, 1, 1),
    character_test = c("z", "z", "a", "a"),
    continuous = 1:4
  )
  out <- suppressWarnings(DiagnosticLikelihoodRatioTable(df_Test, "outcome", c("numeric_test", "character_test")))
  expect_equal(out$Metadata$OutcomeLevels$outcome$positive, "TRUE")
  expect_equal(out$Metadata$BinaryPredictorLevels$numeric_test$positive, "1")
  expect_equal(out$Metadata$BinaryPredictorLevels$character_test$positive, "z")
  expect_true(any(is.infinite(out$Results$LikelihoodRatio)))
  expect_true(all(is.na(out$Results$LRLowerCI)))
  corrected <- DiagnosticLikelihoodRatioTable(df_Test, "outcome", "numeric_test", continuity_correction = 0.5)
  expect_true(all(corrected$Results$Corrected))
  expect_true(all(is.finite(corrected$Results$LRLowerCI)))
  expect_error(DiagnosticLikelihoodRatioTable(df_Test, "outcome", "continuous"), "does not automatically dichotomize")
})

test_that("PlotDiagnosticLRHeatmap preserves static diagnostic data", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    disease = factor(rep(c("No", "Yes"), each = 4), levels = c("No", "Yes")),
    test = factor(c("Negative", "Negative", "Positive", "Positive", "Negative", "Positive", "Positive", "Positive"), levels = c("Negative", "Positive")),
    genotype = factor(rep(c("e2", "e2", "e3", "e3"), 2)),
    site = rep(c("A", "B"), 4)
  )
  out <- DiagnosticLikelihoodRatioTable(df_Test, "disease", c("test", "genotype"), stratify_by = "site", continuity_correction = 0.5)
  plot_all <- PlotDiagnosticLRHeatmap(out, result = "all", show_values = "auto", facet_strata = TRUE)
  expect_s3_class(plot_all, "ggplot")
  expect_true("HoverText" %in% names(attr(plot_all, "DiagnosticLRData")))
  expect_true(any(grepl("unadjusted", plot_all$labels$caption, ignore.case = TRUE)))
  plot_positive <- PlotDiagnosticLRHeatmap(out, result = "positive", predictor_order = "cluster", outcome_order = "cluster")
  expect_s3_class(plot_positive, "ggplot")
  expect_error(PlotDiagnosticLRHeatmap(out$Results, result = "positive"), "requires binary")
})
