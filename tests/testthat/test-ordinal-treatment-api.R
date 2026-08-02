test_that("PlotZScore uses TreatOrdinalAs and retains deprecated Ordinal support", {
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 6)),
    severity = ordered(rep(c("None", "Mild", "Moderate"), 4),
                       levels = c("None", "Mild", "Moderate")),
    marker = seq_len(12)
  )

  continuous <- PlotZScore(
    df_Test, TargetVar = "group", variables = c("severity", "marker"),
    TreatOrdinalAs = "Continuous"
  )
  excluded <- PlotZScore(
    df_Test, TargetVar = "group", variables = c("severity", "marker"),
    TreatOrdinalAs = "Exclude"
  )

  expect_s3_class(continuous, "ggplot")
  expect_s3_class(excluded, "ggplot")
  expect_error(
    PlotZScore(df_Test, TargetVar = "group", variables = "severity", TreatOrdinalAs = "Categorical"),
    "requires TreatOrdinalAs"
  )
  expect_warning(
    PlotZScore(df_Test, TargetVar = "group", variables = "severity", Ordinal = TRUE),
    "deprecated"
  )
})

test_that("PlotInteractionEffectsMatrix uses model-safe ordinal treatments", {
  set.seed(1)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 12)),
    outcome = stats::rnorm(24),
    severity = ordered(rep(c("None", "Mild", "Moderate"), 8),
                       levels = c("None", "Mild", "Moderate")),
    marker = stats::rnorm(24)
  )

  continuous <- PlotInteractionEffectsMatrix(
    df_Test, interVar = "group", outcome_vars = "outcome",
    predictor_vars = c("severity", "marker"), TreatOrdinalAs = "Continuous"
  )
  categorical <- PlotInteractionEffectsMatrix(
    df_Test, interVar = "group", outcome_vars = "outcome",
    predictor_vars = c("severity", "marker"), TreatOrdinalAs = "Categorical"
  )
  excluded <- PlotInteractionEffectsMatrix(
    df_Test, interVar = "group", outcome_vars = "outcome",
    predictor_vars = c("severity", "marker"), TreatOrdinalAs = "Exclude"
  )

  expect_true(is.list(continuous))
  expect_true(is.list(categorical))
  expect_true(is.list(excluded))
  expect_error(
    PlotInteractionEffectsMatrix(
      df_Test, interVar = "group", outcome_vars = "outcome",
      predictor_vars = "marker", TreatOrdinalAs = "Both"
    ),
    "not meaningful"
  )
  expect_error(
    PlotInteractionEffectsMatrix(
      df_Test, interVar = "group", outcome_vars = "severity",
      predictor_vars = "marker", TreatOrdinalAs = "Categorical"
    ),
    "Categorical ordinal outcomes"
  )
  expect_warning(
    PlotInteractionEffectsMatrix(
      df_Test, interVar = "group", outcome_vars = "outcome",
      predictor_vars = "marker", Ordinal = TRUE
    ),
    "deprecated"
  )
})
