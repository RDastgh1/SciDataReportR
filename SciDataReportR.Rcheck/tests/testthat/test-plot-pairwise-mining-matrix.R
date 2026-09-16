test_that("PlotPairwiseMiningMatrix returns directed continuous and categorical contrasts", {
  skip_if_not_installed("emmeans")
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
    variable_metadata = tibble::tibble(Variable = c("Marker", "Exposure"), AnchorN = c(60, 60), TimeExtended = c(FALSE, FALSE))
  )

  expect_s3_class(res, "SciDataReportRPairwiseMiningMatrix")
  expect_true(all(c("Continuous", "Categorical", "Omnibus") %in% names(res$Plots)))
  expect_true(all(c("A", "B") %in% res$Results$Group))
  expect_true(any(res$Results$VariableLabel == "Labelled marker"))
  expect_true(any(res$Results$RowLabel == "Labelled exposure: Yes"))
  expect_gt(res$Results$EstimatedMeanDifference[res$Results$Variable == "Marker" & res$Results$Group == "A"], 0)
  expect_gt(res$Results$PercentagePointDifference[res$Results$Variable == "Exposure" & res$Results$Level == "Yes" & res$Results$Group == "A"], 0)
  expect_true(all(c("AnchorN", "TimeExtended") %in% names(res$Results)))
  expect_equal(res$Results$AdjustedPValue, res$Results$PValue)
})

test_that("PlotPairwiseMiningMatrix requires an available referent", {
  df_Test <- data.frame(Group = factor(c("A", "B")), Marker = c(1, 2))
  expect_error(PlotPairwiseMiningMatrix(df_Test, "Group", "Marker", "Control"), "Referent")
})
