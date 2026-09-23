test_that("PlotPairwiseMiningMatrix returns non-directional separation magnitudes", {
  skip_if_not_installed("emmeans")
  set.seed(20260916)
  df_Test <- data.frame(
    Group = factor(rep(c("Referent", "A", "B"), each = 20), levels = c("Referent", "A", "B")),
    Marker = c(rnorm(20, 0, 1), rnorm(20, 1.5, 1), rnorm(20, -1.5, 1)),
    Exposure = factor(c(rep("No", 16), rep("Yes", 4), rep("No", 6), rep("Yes", 14), rep("No", 18), rep("Yes", 2)))
  )
  df_Test$Marker <- sjlabelled::set_label(df_Test$Marker, "Labelled marker")
  df_Test$Exposure <- sjlabelled::set_label(df_Test$Exposure, "Labelled exposure")

  res <- PlotPairwiseMiningMatrix(
    data = df_Test, group_var = "Group", variables = c("Marker", "Exposure"),
    Referent = "Referent", adjust_scope = "none", p_adjust_method = "none",
    variable_metadata = tibble::tibble(Variable = c("Marker", "Exposure"), AnchorN = c(60, 60), TimeExtended = c(FALSE, TRUE))
  )

  expect_s3_class(res, "SciDataReportRPairwiseMiningMatrix")
  expect_setequal(names(res$Plots), c("Continuous", "Categorical"))
  expect_true(all(c("A", "B") %in% res$Results$Group))
  expect_true(any(res$Results$VariableLabel == "Labelled marker"))
  expect_equal(nrow(dplyr::filter(res$Results, .data$Variable == "Exposure")), 2)
  expect_true(all(dplyr::filter(res$Results, .data$Variable == "Exposure")$EffectMagnitude >= 0))
  expect_true(all(dplyr::filter(res$Results, .data$Variable == "Exposure")$CramersV >= 0))
  expect_gt(dplyr::filter(res$Results, .data$Variable == "Marker", .data$Group == "A")$SignedEstimatedMeanDifference, 0)
  expect_lt(dplyr::filter(res$Results, .data$Variable == "Marker", .data$Group == "B")$SignedEstimatedMeanDifference, 0)
  expect_true(all(dplyr::filter(res$Results, .data$Variable == "Marker")$EffectMagnitude == abs(dplyr::filter(res$Results, .data$Variable == "Marker")$SignedEstimatedMeanDifference)))
  expect_true(all(c("AnchorN", "TimeExtended") %in% names(res$Results)))
  expect_true("IsAdjustedSignificant" %in% names(res$Results))
  expect_true(any(grepl("†", res$Results$DisplayRowLabel, fixed = TRUE)))
  expect_equal(res$Results$AdjustedPValue, res$Results$PValue)
  expect_equal(unname(res$Settings$colors), c("#F1F5F9", "#2166AC"))
  expect_true(nrow(res$CategoryLevelResults) > 0)
  expect_true(all(grepl(":", res$CategoryLevelResults$RowLabel, fixed = TRUE)))
})

test_that("PlotPairwiseMiningMatrix applies per-contrast FDR outlines", {
  skip_if_not_installed("emmeans")
  set.seed(20260916)
  df_Test <- data.frame(
    Group = factor(rep(c("Referent", "A", "B"), each = 30), levels = c("Referent", "A", "B")),
    Strong = c(rnorm(30, 0, 1), rnorm(30, 2.5, 1), rnorm(30, -2.5, 1)),
    Weak = rnorm(90, 0, 1)
  )
  res <- PlotPairwiseMiningMatrix(
    data = df_Test,
    group_var = "Group",
    variables = c("Strong", "Weak"),
    Referent = "Referent",
    adjust_scope = "per_group",
    p_adjust_method = "fdr",
    adjusted_outline = TRUE
  )
  expected <- rep(NA_real_, nrow(res$Results))
  for (idx in split(seq_len(nrow(res$Results)), res$Results$Group)) {
    expected[idx] <- stats::p.adjust(res$Results$PValue[idx], method = "fdr")
  }
  expect_equal(res$Results$AdjustedPValue, expected, tolerance = 1e-12)
  expect_equal(
    res$Results$IsAdjustedSignificant,
    res$Results$AdjustedPValue <= 0.05
  )
  expect_true(any(res$Results$IsAdjustedSignificant))
  expect_equal(res$Settings$adjusted_outline, TRUE)
  expect_gte(length(res$Plots$Continuous$layers), 3)

  res_no_outline <- PlotPairwiseMiningMatrix(
    data = df_Test,
    group_var = "Group",
    variables = c("Strong", "Weak"),
    Referent = "Referent",
    adjust_scope = "per_group",
    p_adjust_method = "fdr",
    adjusted_outline = FALSE
  )
  expect_equal(res_no_outline$Settings$adjusted_outline, FALSE)
  expect_lt(length(res_no_outline$Plots$Continuous$layers), length(res$Plots$Continuous$layers))
})

test_that("PlotPairwiseMiningMatrix requires an available referent", {
  df_Test <- data.frame(Group = factor(c("A", "B")), Marker = c(1, 2))
  expect_error(PlotPairwiseMiningMatrix(df_Test, "Group", "Marker", "Control"), "Referent")
})
