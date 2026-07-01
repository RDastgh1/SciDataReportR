test_that("MultivariableRegressionTable returns a stable ordinary regression object", {
  skip_if_not_installed("pROC")

  set.seed(104)
  df <- data.frame(
    y = rnorm(80),
    ybin = factor(sample(c("No", "Yes"), 80, replace = TRUE)),
    x1 = rnorm(80),
    x2 = rnorm(80),
    cov = rnorm(80)
  )
  attr(df$x1, "label") <- "Marker one"

  res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = c("y", "ybin"),
    PredictorVars = c("x1", "x2"),
    Covars = "cov",
    Method = "lm"
  )

  expect_named(
    res,
    c(
      "Models", "FormattedTable", "LargeTable", "RegressionMatrix",
      "VariableImportanceMatrix", "Predictions", "Diagnostics",
      "ModelSummary", "Multicollinearity", "Plots", "Metadata"
    )
  )
  expect_s3_class(res$FormattedTable, "data.frame")
  expect_s3_class(res$LargeTable, "data.frame")
  expect_s3_class(res$RegressionMatrix, "data.frame")
  expect_s3_class(res$VariableImportanceMatrix, "data.frame")
  expect_s3_class(res$Predictions, "data.frame")
  expect_s3_class(res$Diagnostics, "data.frame")
  expect_s3_class(res$ModelSummary, "data.frame")

  expect_true(all(c("OutcomeIndex", "PredictorIndex", "HoverText") %in% names(res$RegressionMatrix)))
  expect_true(all(c("VariableImportance", "VariableImportanceType") %in% names(res$VariableImportanceMatrix)))
  expect_true(all(c("Leverage", "CooksDistance", "ClassificationThreshold") %in% names(res$Predictions)))
  expect_true(any(res$RegressionMatrix$PredictorLabel == "Marker one"))
  expect_true(any(grepl("<b>Outcome:</b>", res$RegressionMatrix$HoverText, fixed = TRUE)))

  expect_true(all(c("MissingRemoved", "PercentRemoved", "MaximumVIF", "MaximumCorrelation") %in% names(res$ModelSummary)))
  expect_true("CorrelationDataFrame" %in% names(res$Multicollinearity))
  expect_true("HighVIFVariables" %in% names(res$Multicollinearity))
})

test_that("MultivariableRegressionTable stores logistic diagnostics and predictions", {
  skip_if_not_installed("pROC")

  set.seed(205)
  x1 <- rnorm(120)
  x2 <- rnorm(120)
  probability <- stats::plogis(-0.2 + 0.9 * x1 - 0.4 * x2)
  df <- data.frame(
    ybin = factor(ifelse(stats::runif(120) < probability, "Yes", "No")),
    x1 = x1,
    x2 = x2
  )

  res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "ybin",
    PredictorVars = c("x1", "x2"),
    Method = "lm"
  )

  expect_equal(res$RegressionMatrix$OutcomeFamily, rep("logistic", 2))
  expect_true(all(c("AUC", "AUCLowerCI", "AUCUpperCI", "OptimalThreshold") %in% names(res$Diagnostics)))
  expect_false(is.na(res$Diagnostics$OptimalThreshold))
  expect_false(all(is.na(res$Predictions$PredictedProbability)))
  expect_false(all(is.na(res$Predictions$PredictedClass)))
  expect_equal(unique(stats::na.omit(res$Predictions$ClassificationThreshold)), res$Diagnostics$OptimalThreshold)
  expect_equal(unique(res$VariableImportanceMatrix$VariableImportanceType), "Likelihood Ratio Chi-square")
})

test_that("MultivariableRegressionTable stores penalized regression tuning", {
  skip_if_not_installed("glmnet")

  set.seed(306)
  df <- data.frame(
    y = rnorm(70),
    x1 = rnorm(70),
    x2 = rnorm(70),
    cov = rnorm(70)
  )

  res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "y",
    PredictorVars = c("x1", "x2"),
    Covars = "cov",
    Method = "elasticnet",
    CVFolds = 5
  )

  expect_true(all(is.na(res$RegressionMatrix$PValue)))
  expect_true(all(res$VariableImportanceMatrix$VariableImportanceType == "Absolute Standardized Beta"))
  expect_true(all(c("AlphaGrid", "CrossValidationErrors", "BestAlpha", "BestLambda") %in% names(res$Metadata$AnalysisSettings$Tuning$y)))
  expect_s3_class(res$Metadata$AnalysisSettings$Tuning$y$CrossValidationErrors, "data.frame")
  expect_true("RetainedPredictors" %in% names(res$Diagnostics))
  expect_type(res$Diagnostics$RetainedPredictors[[1]], "character")
})

test_that("MultivariableRegressionTable validates inputs", {
  df <- data.frame(
    y = rnorm(20),
    bad = letters[1:20],
    x1 = rnorm(20)
  )

  expect_error(
    MultivariableRegressionTable(df, "missing", "x1"),
    "not found"
  )
  expect_error(
    MultivariableRegressionTable(df, "bad", "x1"),
    "must be numeric, logical, or a two-level factor"
  )
})
