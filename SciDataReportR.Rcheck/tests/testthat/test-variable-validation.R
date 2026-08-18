test_that("ScidrValidateVariables names the offending argument", {
  df_Test <- data.frame(a = 1:5, b = 6:10)

  expect_error(
    ScidrValidateVariables(df_Test, c("a", "NOPE"), "predictor_vars"),
    "Variables not found in `data`: NOPE \\(supplied to `predictor_vars`\\)",
    fixed = FALSE
  )
  expect_equal(ScidrValidateVariables(df_Test, c("a", "b")), c("a", "b"))
  expect_equal(ScidrValidateVariables(df_Test, NULL), character(0))
  expect_equal(ScidrValidateVariables(df_Test, c("a", "a")), "a")
  expect_equal(
    ScidrValidateVariables(df_Test, c("a", "a"), unique_only = FALSE),
    c("a", "a")
  )
  expect_error(
    ScidrValidateVariables(df_Test, NULL, "variables", allow_null = FALSE),
    "must name at least one column"
  )
  expect_error(ScidrValidateVariables(list(a = 1), "a"), "must be a data frame")
})

test_that("ScidrValidateVariable requires exactly one column", {
  df_Test <- data.frame(a = 1:5, b = 6:10)

  expect_equal(ScidrValidateVariable(df_Test, "a", "interVar"), "a")
  expect_error(
    ScidrValidateVariable(df_Test, c("a", "b"), "interVar"),
    "must name exactly one column"
  )
  expect_error(ScidrValidateVariable(df_Test, NULL, "interVar"), "must name a column")
  expect_null(ScidrValidateVariable(df_Test, NULL, "interVar", allow_null = TRUE))
})

# A misspelled variable used to behave five different ways across this family:
# silently dropped (shrinking the result matrix), or a bare base-R "undefined
# columns selected". A misspelled covariate was worse -- it silently produced
# unadjusted statistics reported as adjusted.
test_that("matrix family rejects unknown variables with one message", {
  skip_if_not_installed("rstatix")
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Test <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  Calls <- list(
    PlotCorrelationsHeatmap = function() PlotCorrelationsHeatmap(
      df_Test, predictor_vars = c("age", "NOPE"), outcome_vars = "age"),
    PlotAnovaRelationshipsMatrix = function() PlotAnovaRelationshipsMatrix(
      df_Test, CatVars = c("sex", "NOPE"), ContVars = "age"),
    PlotChiSqCovar = function() PlotChiSqCovar(
      df_Test, predictor_vars = c("sex", "NOPE"), outcome_vars = "sex"),
    PlotPointCorrelationsHeatmap = function() PlotPointCorrelationsHeatmap(
      df_Test, CatVars = c("sex", "NOPE"), ContVars = "age"),
    PlotPhiHeatmap = function() PlotPhiHeatmap(
      df_Test, CatVars = c("sex", "NOPE")),
    PlotMiningMatrix = function() PlotMiningMatrix(
      df_Test, outcome_vars = c("age", "NOPE"), predictor_vars = "AXL"),
    PlotDirectionalHeatmaps = function() PlotDirectionalHeatmaps(
      df_Test, variables = c("age", "NOPE")),
    PlotNumInteractionEffectsMatrix = function() PlotNumInteractionEffectsMatrix(
      df_Test, predictor_vars = c("AXL", "NOPE"), outcome_vars = "age",
      interVar = "Cortisol"),
    PlotCatInteractionEffectsMatrix = function() PlotCatInteractionEffectsMatrix(
      df_Test, predictor_vars = c("AXL", "NOPE"), outcome_vars = "age",
      interVar = "sex")
  )

  for (nm in names(Calls)) {
    expect_error(Calls[[nm]](), "Variables not found in `data`: NOPE",
                 fixed = TRUE, label = nm)
  }
})

test_that("matrix family rejects unknown covariates rather than silently unadjusting", {
  skip_if_not_installed("rstatix")
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Test <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  Calls <- list(
    PlotCorrelationsHeatmap = function() PlotCorrelationsHeatmap(
      df_Test, predictor_vars = "age", outcome_vars = "AXL", covariates = "NOPE"),
    PlotChiSqCovar = function() PlotChiSqCovar(
      df_Test, predictor_vars = "sex", outcome_vars = "Diagnosis", covariates = "NOPE"),
    PlotAnovaRelationshipsMatrix = function() PlotAnovaRelationshipsMatrix(
      df_Test, CatVars = "sex", ContVars = "age", covariates = "NOPE"),
    PlotPointCorrelationsHeatmap = function() PlotPointCorrelationsHeatmap(
      df_Test, CatVars = "sex", ContVars = "age", covariates = "NOPE"),
    PlotMiningMatrix = function() PlotMiningMatrix(
      df_Test, outcome_vars = "age", predictor_vars = "AXL", covariates = "NOPE"),
    PlotNumInteractionEffectsMatrix = function() PlotNumInteractionEffectsMatrix(
      df_Test, predictor_vars = "AXL", outcome_vars = "age",
      interVar = "Cortisol", covariates = "NOPE")
  )

  for (nm in names(Calls)) {
    expect_error(Calls[[nm]](),
                 "Variables not found in `data`: NOPE (supplied to `covariates`)",
                 fixed = TRUE, label = nm)
  }
})

test_that("valid variables and covariates are unaffected", {
  skip_if_not_installed("rstatix")
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Test <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  Vars <- c("age", "AXL", "Adiponectin")

  Square <- PlotCorrelationsHeatmap(df_Test, predictor_vars = Vars, outcome_vars = Vars)
  expect_equal(dim(Square$Unadjusted$p), c(3L, 3L))
  expect_equal(rownames(Square$Unadjusted$p), Vars)

  Adjusted <- PlotCorrelationsHeatmap(
    df_Test, predictor_vars = Vars, outcome_vars = Vars, covariates = "Cortisol")
  expect_equal(dim(Adjusted$Unadjusted$p), c(3L, 3L))

  # NULL variables still auto-detect rather than erroring
  expect_gt(nrow(PlotCorrelationsHeatmap(df_Test)$Unadjusted$p), 0)

  expect_s3_class(
    PlotPhiHeatmap(df_Test, CatVars = c("sex", "Diagnosis"))$Unadjusted$plot,
    "ggplot")
  expect_s3_class(
    PlotDirectionalHeatmaps(df_Test, variables = Vars)$Unadjusted$plot, "ggplot")
})

test_that("PlotDirectionalHeatmaps keeps deprecated formals trailing", {
  # yVars used to sit in formal slot 3, so PlotDirectionalHeatmaps(df, vars, TRUE)
  # bound TRUE to yVars instead of Relabel.
  Formals <- names(formals(PlotDirectionalHeatmaps))
  expect_equal(Formals[1:3], c("data", "variables", "Relabel"))
  expect_equal(tail(Formals, 3), c("Data", "xVars", "yVars"))
})
