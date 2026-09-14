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
  expect_named(plot_all, c("DiagnosticLR", "DiagnosticMatrices", "DiagnosticLRData", "DiagnosticMatrixData"))
  expect_s3_class(plot_all$DiagnosticLR, "ggplot")
  expect_s3_class(plot_all$DiagnosticMatrices, "ggplot")
  expect_true("HoverText" %in% names(attr(plot_all$DiagnosticLR, "DiagnosticLRData")))
  expect_true(any(grepl("unadjusted", plot_all$DiagnosticLR$labels$caption, ignore.case = TRUE)))
  expect_equal(nrow(plot_all$DiagnosticMatrixData), 16)
  expect_true(all(c("Outcome negative", "Outcome positive") %in% plot_all$DiagnosticMatrixData$OutcomeCondition))
  plot_positive <- PlotDiagnosticLRHeatmap(out, result = "positive", predictor_order = "cluster", outcome_order = "cluster")
  expect_s3_class(plot_positive$DiagnosticLR, "ggplot")
  expect_s3_class(plot_positive$DiagnosticMatrices, "ggplot")
  expect_error(PlotDiagnosticLRHeatmap(out$Results, result = "positive"), "requires binary")
})

test_that("PlotDiagnosticLRHeatmap splits multi-outcome diagnostic displays", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    outcome_a = factor(rep(c("No", "Yes"), each = 8), levels = c("No", "Yes")),
    outcome_b = factor(rep(c("No", "Yes"), 8), levels = c("No", "Yes")),
    test = factor(rep(c("Negative", "Positive"), 8), levels = c("Negative", "Positive")),
    genotype = factor(rep(c("e2", "e3", "e4", "e3"), 4))
  )
  out <- DiagnosticLikelihoodRatioTable(
    df_Test,
    outcome_vars = c("outcome_a", "outcome_b"),
    predictor_vars = c("test", "genotype"),
    continuity_correction = 0.5
  )
  split <- PlotDiagnosticLRHeatmap(out)
  expect_named(split, c("Overview", "ByOutcome", "DiagnosticMatrices", "OverviewData", "ByOutcomeData", "DiagnosticMatrixData"))
  expect_s3_class(split$Overview, "ggplot")
  expect_s3_class(split$ByOutcome, "ggplot")
  expect_s3_class(split$DiagnosticMatrices, "ggplot")
  expect_equal(nrow(split$OverviewData), 4)
  expect_true(all(c("test", "genotype") %in% split$OverviewData$Predictor))
  fill_scale <- split$Overview$scales$get_scales("fill")
  expect_true(any(grepl("LR +1", fill_scale$get_labels(c(-1, 0, 1)))))
  expect_named(PlotDiagnosticLRHeatmap(out, multi_outcome = "combined"), c("DiagnosticLR", "DiagnosticMatrices", "DiagnosticLRData", "DiagnosticMatrixData"))
  expect_named(PlotDiagnosticLRHeatmap(out, multi_outcome = "split"), c("Overview", "ByOutcome", "DiagnosticMatrices", "OverviewData", "ByOutcomeData", "DiagnosticMatrixData"))
})

test_that("PlotDiagnosticLRForest uses existing diagnostic LR results", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    outcome_a = factor(rep(c("No", "Yes"), each = 8), levels = c("No", "Yes")),
    outcome_b = factor(rep(c("No", "Yes"), 8), levels = c("No", "Yes")),
    test = factor(rep(c("Negative", "Positive"), 8), levels = c("Negative", "Positive")),
    genotype = factor(rep(c("e2", "e3", "e4", "e3"), 4))
  )
  out <- DiagnosticLikelihoodRatioTable(
    df_Test,
    outcome_vars = c("outcome_a", "outcome_b"),
    predictor_vars = c("test", "genotype"),
    continuity_correction = 0.5
  )
  p <- PlotDiagnosticLRForest(out)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("DiagnosticLRForestData", "DiagnosticLRForestLimits", "DiagnosticLRForestFacetBy") %in% names(attributes(p))))
  expect_equal(attr(p, "DiagnosticLRForestFacetBy"), "outcome")
  expect_true(all(attr(p, "DiagnosticLRForestLimits") > 0))
  expect_true(any(grepl("unadjusted", p$labels$caption, ignore.case = TRUE)))
  expect_silent(ggplot2::ggplot_build(p))

  p_clipped <- PlotDiagnosticLRForest(out, limits = c(0.5, 2))
  expect_true(any(attr(p_clipped, "DiagnosticLRForestData")$ArrowLeft | attr(p_clipped, "DiagnosticLRForestData")$ArrowRight))

  p_positive <- PlotDiagnosticLRForest(out, result = "positive", facet_by = "predictor")
  expect_true(all(attr(p_positive, "DiagnosticLRForestData")$Predictor == "test"))
  expect_equal(attr(p_positive, "DiagnosticLRForestFacetBy"), "predictor")
  expect_silent(ggplot2::ggplot_build(p_positive))

  df_Test$site <- rep(c("A", "B", "B", "A"), 4)
  out_stratified <- DiagnosticLikelihoodRatioTable(
    df_Test,
    outcome_vars = "outcome_a", predictor_vars = "test",
    stratify_by = "site", continuity_correction = 0.5
  )
  p_stratified <- PlotDiagnosticLRForest(out_stratified)
  expect_equal(dplyr::n_distinct(attr(p_stratified, "DiagnosticLRForestData")$Stratum), 2)
  expect_silent(ggplot2::ggplot_build(p_stratified))

  expect_error(PlotDiagnosticLRForest(out$Results, result = "positive"), "requires binary")
  expect_error(PlotDiagnosticLRForest(out, limits = c(1, 1)), "increasing positive")
})

test_that("PlotDiagnosticLRForest renders zero and infinite estimates at boundaries", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(
    disease = factor(c("No", "No", "Yes", "Yes"), levels = c("No", "Yes")),
    test = factor(c("Negative", "Negative", "Positive", "Positive"), levels = c("Negative", "Positive"))
  )
  out <- suppressWarnings(DiagnosticLikelihoodRatioTable(df_Test, "disease", "test"))
  p <- PlotDiagnosticLRForest(out, limits = c(0.25, 4))
  df_Plot <- attr(p, "DiagnosticLRForestData")
  expect_true(all(c("Zero", "Infinite") %in% df_Plot$EstimateStatus))
  expect_true(all(df_Plot$BoundaryLabel %in% c("0", "Inf")))
  expect_true(all(is.na(df_Plot$PlotLowerCI)))
  expect_true(all(is.na(df_Plot$PlotUpperCI)))
  expect_silent(ggplot2::ggplot_build(p))
})
