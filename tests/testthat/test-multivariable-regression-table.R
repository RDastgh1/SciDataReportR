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
  expect_s3_class(res$FormattedTable, "gt_tbl")
  expect_s3_class(res$LargeTable, "data.frame")
  expect_s3_class(res$RegressionMatrix, "data.frame")
  expect_s3_class(res$VariableImportanceMatrix, "data.frame")
  expect_s3_class(res$Predictions, "data.frame")
  expect_s3_class(res$Diagnostics, "data.frame")
  expect_s3_class(res$ModelSummary, "data.frame")
  expect_true(length(res$Plots) > 0)
  expect_s3_class(res$Plots$RegressionMatrix, "ggplot")
  expect_s3_class(res$Plots$VariableImportanceMatrix, "ggplot")
  expect_s3_class(res$Plots$ObservedVsPredicted, "ggplot")

  expect_true(all(c("OutcomeIndex", "PredictorIndex", "HoverText") %in% names(res$RegressionMatrix)))
  expect_true(all(c("VariableImportance", "VariableImportanceType") %in% names(res$VariableImportanceMatrix)))
  expect_true(all(c("Leverage", "CooksDistance", "ClassificationThreshold") %in% names(res$Predictions)))
  expect_true(any(res$RegressionMatrix$PredictorLabel == "Marker one"))
  expect_true(any(grepl("<b>Outcome:</b>", res$RegressionMatrix$HoverText, fixed = TRUE)))

  expect_true(all(c("MissingRemoved", "PercentRemoved", "MaximumVIF", "MaximumCorrelation") %in% names(res$ModelSummary)))
  expect_true("CorrelationDataFrame" %in% names(res$Multicollinearity))
  expect_true("HighVIFVariables" %in% names(res$Multicollinearity))
  expect_true(res$Metadata$AnalysisSettings$Standardize)
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

test_that("MultivariableRegressionTable treats covariates as mandatory (unpenalized) adjustments", {
  skip_if_not_installed("glmnet")

  set.seed(512)
  n <- 80
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    cov = rnorm(n),
    cov_cat = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  df$y <- 0.8 * df$x1 + 0.5 * df$cov + rnorm(n)

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    covariates = c("cov", "cov_cat"),
    Method = "lasso",
    CVFolds = 5
  )

  tuning <- res$Metadata$AnalysisSettings$Tuning$y
  expect_true(all(c("PenaltyFactors", "UnpenalizedTerms") %in% names(tuning)))
  expect_setequal(tuning$UnpenalizedTerms, c("cov", "cov_cat"))

  # Covariate columns (including factor dummies) get penalty factor 0,
  # predictor columns get 1.
  pf <- tuning$PenaltyFactors
  expect_equal(unname(pf[c("x1", "x2")]), c(1, 1))
  expect_true(all(pf[startsWith(names(pf), "cov")] == 0))

  # Unpenalized covariates are always retained by the lasso.
  # (Covariate rows live in LargeTable; RegressionMatrix reports predictors only.)
  cov_rows <- res$LargeTable$Predictor %in% c("cov", "cov_cat")
  expect_true(any(cov_rows))
  expect_true(all(res$LargeTable$Selected[cov_rows]))
})

test_that("MultivariableRegressionTable FormattedTable is a gt table matching the univariate style", {
  skip_if_not_installed("gt")
  skip_if_not_installed("pROC")

  set.seed(618)
  df <- data.frame(
    y = rnorm(80),
    x1 = rnorm(80),
    x2 = rnorm(80),
    cov = rnorm(80)
  )

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    covariates = "cov",
    Method = "lm"
  )

  expect_s3_class(res$FormattedTable, "gt_tbl")

  # The rendered table exposes a single combined estimate cell and a p-value
  # column, and reports predictor rows only (covariates are adjusted for but
  # not shown), mirroring MakeUnivariateRegressionTable's FormattedTable.
  gt_data <- res$FormattedTable[["_data"]]
  expect_true(all(c("Estimate_CI", "P") %in% names(gt_data)))
  expect_true(all(grepl("(", gt_data$Estimate_CI, fixed = TRUE)))
  expect_equal(nrow(gt_data), 2L)
})

test_that("MultivariableRegressionTable flags separation and blanks unreliable logistic estimates", {
  skip_if_not_installed("pROC")

  # x1 perfectly separates the outcome -> quasi-complete separation.
  set.seed(70)
  df <- data.frame(
    grp = factor(c(rep("No", 20), rep("Yes", 20))),
    x1 = c(rep(-10, 20), rep(10, 20)),
    x2 = rnorm(40)
  )

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "grp",
    predictor_vars = c("x1", "x2"),
    Method = "lm"
  )

  expect_true(res$ModelSummary$SeparationDetected)
  expect_false(res$ModelSummary$Converged)
  # Unreliable coefficients are blanked, not left to explode into the plot.
  expect_true(all(is.na(res$RegressionMatrix$StandardizedBeta)))
  expect_true(any(grepl("Separation", res$Diagnostics$Warnings, ignore.case = TRUE)))
})

test_that("MultivariableRegressionTable reports perfectly collinear (aliased) terms", {
  skip_if_not_installed("pROC")

  set.seed(88)
  n <- 60
  df <- data.frame(y = rnorm(n), a = rnorm(n))
  df$b <- df$a       # exact duplicate -> aliased
  df$c <- 2 * df$a   # linear combination -> aliased

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("a", "b", "c"),
    Method = "lm"
  )

  expect_gt(res$ModelSummary$AliasedTermCount, 0)
  expect_true(any(res$LargeTable$Aliased))
  expect_true(any(grepl("aliased", res$Diagnostics$Warnings, ignore.case = TRUE)))
})

test_that("MultivariableRegressionTable reports an omnibus model p-value with a matching column annotation", {
  skip_if_not_installed("pROC")

  set.seed(131)
  n <- 90
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- 1.2 * df$x1 + rnorm(n)  # strong linear signal -> significant omnibus F

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    Method = "lm"
  )

  expect_equal(res$ModelSummary$ModelStatType, "F")
  expect_true(is.finite(res$ModelSummary$ModelPValue))
  expect_lt(res$ModelSummary$ModelPValue, 0.05)

  ann <- SciDataReportR:::ScidrColumnAnnotations(res$ModelSummary)
  expect_true(any(grepl("^p", ann$Label)))
})

test_that("Penalized column annotation reports deviance explained, not a p-value", {
  skip_if_not_installed("glmnet")

  set.seed(202)
  n <- 90
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- 1.0 * df$x1 + rnorm(n)

  res <- MultivariableRegressionTable(
    data = df,
    outcome_vars = "y",
    predictor_vars = c("x1", "x2"),
    Method = "lasso",
    CVFolds = 5
  )

  expect_true(is.na(res$ModelSummary$ModelPValue))
  ann <- SciDataReportR:::ScidrColumnAnnotations(res$ModelSummary)
  expect_true(any(grepl("Dev. expl", ann$Label, fixed = TRUE)))
})

test_that("ScidrRobustFillLimits clamps extreme values to symmetric limits", {
  # A lone extreme should not blow out the range set by the bulk of the data.
  vals <- c(rnorm(100, sd = 0.3), 600)
  lims <- SciDataReportR:::ScidrRobustFillLimits(vals)
  expect_length(lims, 2)
  expect_equal(lims[1], -lims[2])
  expect_lt(lims[2], 600)

  expect_null(SciDataReportR:::ScidrRobustFillLimits(c(NA_real_, NA_real_)))
})

test_that("MultivariableRegressionTable drops sparse predictors and imputes by default", {
  skip_if_not_installed("pROC")

  set.seed(407)
  df <- data.frame(
    cohort = factor(sample(c("Control", "Long Covid"), 90, replace = TRUE)),
    x1 = rnorm(90),
    x2 = rnorm(90),
    x_sparse = rnorm(90)
  )
  df$x1[sample(seq_len(90), 12)] <- NA
  df$x2[sample(seq_len(90), 10)] <- NA
  df$x_sparse[sample(seq_len(90), 50)] <- NA

  res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "cohort",
    PredictorVars = c("x1", "x2", "x_sparse")
  )

  expect_equal(sort(res$RegressionMatrix$Predictor), c("x1", "x2"))
  expect_equal(res$ModelSummary$DroppedPredictorCount, 1)
  expect_equal(res$ModelSummary$ImputedPredictorCount, 2)
  expect_equal(res$Metadata$Missingness$Outcomes$DroppedVariables[[1]], "x_sparse")
  expect_equal(sort(res$Metadata$Missingness$Outcomes$ImputedVariables[[1]]), c("x1", "x2"))
  expect_equal(res$Metadata$AnalysisSettings$MissingDataStrategy, "drop_sparse_impute")
})

test_that("MultivariableRegressionTable preserves strict complete-case behavior when requested", {
  skip_if_not_installed("pROC")

  set.seed(508)
  df <- data.frame(
    cohort = factor(sample(c("Control", "Long Covid"), 80, replace = TRUE)),
    x1 = rnorm(80),
    x2 = rnorm(80)
  )
  df$x1[sample(seq_len(80), 20)] <- NA
  df$x2[sample(seq_len(80), 20)] <- NA

  res <- MultivariableRegressionTable(
    Data = df,
    OutcomeVars = "cohort",
    PredictorVars = c("x1", "x2"),
    MissingDataStrategy = "complete_cases"
  )

  expect_true(res$ModelSummary$SampleSize < nrow(df))
  expect_equal(res$ModelSummary$ImputedPredictorCount, 0)
  expect_equal(length(res$Metadata$Missingness$Outcomes$ImputedVariables[[1]]), 0)
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

test_that("MultivariableRegressionTable supports ordinary multinomial outcomes", {
  df_three <- data.frame(
    cohort = factor(rep(c("A", "B", "C"), each = 10)),
    x1 = rnorm(30)
  )

  res <- MultivariableRegressionTable(df_three, "cohort", "x1")

  expect_equal(unique(res$RegressionMatrix$OutcomeMode), "multinomial")
  expect_equal(unique(res$RegressionMatrix$ReferenceLevel), "A")
  expect_equal(unique(res$RegressionMatrix$OutcomeLevel), c("B", "C"))
  expect_equal(unique(res$RegressionMatrix$ComparisonLabel), c("cohort: B vs A", "cohort: C vs A"))
  expect_true(all(res$RegressionMatrix$EffectType == "Odds Ratio"))
  expect_s3_class(res$Plots$RegressionMatrix, "ggplot")
  expect_equal(res$Metadata$Outcomes$Engine, "nnet::multinom")
})

test_that("MultivariableRegressionTable validates missingness", {

  df_one_level <- data.frame(
    cohort = factor(c(rep("A", 10), rep("B", 10)), levels = c("A", "B")),
    x1 = c(rnorm(10), rep(NA_real_, 10))
  )

  expect_error(
    MultivariableRegressionTable(
      df_one_level,
      "cohort",
      "x1",
      MissingDataStrategy = "complete_cases"
    ),
    "only one level after missing-data handling"
  )

  df_all_dropped <- data.frame(
    cohort = factor(sample(c("A", "B"), 30, replace = TRUE)),
    x1 = c(rnorm(5), rep(NA_real_, 25))
  )

  expect_error(
    MultivariableRegressionTable(df_all_dropped, "cohort", "x1"),
    "All PredictorVars were dropped for missingness"
  )

  df_too_small <- data.frame(
    y = rnorm(10),
    x1 = c(rnorm(3), rep(NA_real_, 7)),
    x2 = c(rnorm(3), rep(NA_real_, 7))
  )

  expect_error(
    MultivariableRegressionTable(
      df_too_small,
      "y",
      c("x1", "x2"),
      MissingDataStrategy = "complete_cases",
      MinCompleteCases = 5
    ),
    "does not have enough rows after missing-data handling"
  )
})

test_that("MultivariableRegressionTable supports ordered proportional-odds outcomes", {
  set.seed(641)
  df <- data.frame(
    severity = ordered(
      sample(c("Mild", "Moderate", "Severe"), 90, replace = TRUE),
      levels = c("Mild", "Moderate", "Severe")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )

  res <- MultivariableRegressionTable(df, "severity", c("x1", "x2"), Method = "lm")

  expect_equal(unique(res$RegressionMatrix$OutcomeMode), "ordinal")
  expect_equal(unique(res$RegressionMatrix$Contrast), "higher vs lower")
  expect_equal(unique(res$RegressionMatrix$ComparisonLabel), "severity: higher vs lower")
  expect_equal(res$Metadata$Outcomes$Engine, "MASS::polr")
  expect_true(res$Metadata$Outcomes$Ordered)
  expect_true(isTRUE(res$Metadata$AnalysisSettings$Tuning$severity$ProportionalOdds))
})

test_that("MultivariableRegressionTable supports grouped penalized multinomial outcomes", {
  skip_if_not_installed("glmnet")
  set.seed(642)
  df <- data.frame(
    cohort = factor(sample(c("A", "B", "C"), 120, replace = TRUE)),
    x1 = rnorm(120),
    x2 = rnorm(120),
    age = rnorm(120, 50, 10)
  )

  res <- MultivariableRegressionTable(
    df, "cohort", c("x1", "x2"), covariates = "age",
    Method = "lasso", CVFolds = 3, reference_levels = c(cohort = "B")
  )

  expect_equal(unique(res$RegressionMatrix$ReferenceLevel), "B")
  expect_equal(sort(unique(res$RegressionMatrix$OutcomeLevel)), c("A", "C"))
  expect_true(all(is.na(res$RegressionMatrix$PValue)))
  expect_equal(res$Metadata$AnalysisSettings$Tuning$cohort$TypeMultinomial, "grouped")
  expect_equal(unname(res$Metadata$AnalysisSettings$Tuning$cohort$PenaltyFactors["age"]), 0)
})

test_that("MultivariableRegressionTable supports one-vs-rest, binary subsets, and skip", {
  set.seed(643)
  df <- data.frame(
    cohort = factor(sample(c("A", "B", "C"), 120, replace = TRUE)),
    y = rnorm(120),
    x1 = rnorm(120)
  )

  ovr <- MultivariableRegressionTable(
    df, c("cohort", "y"), "x1",
    outcome_modes = c(cohort = "one_vs_rest")
  )
  expect_equal(
    unique(ovr$RegressionMatrix$ComparisonLabel[ovr$RegressionMatrix$Outcome == "cohort"]),
    paste0("cohort: ", c("A", "B", "C"), " vs all others")
  )
  expect_true(all(c("nominal", "ordinal", "one_vs_rest", "binary_subset", "skip") %in%
    names(ovr$Metadata$ModelingAdvice)))

  subset <- MultivariableRegressionTable(
    df, "cohort", "x1", outcome_modes = c(cohort = "binary_subset"),
    binary_subsets = list(cohort = c("C", "A"))
  )
  expect_equal(unique(subset$RegressionMatrix$ComparisonLabel), "cohort: A vs C")
  expect_equal(unique(subset$RegressionMatrix$ReferenceLevel), "C")

  expect_warning(
    skipped <- MultivariableRegressionTable(
      df, c("cohort", "y"), "x1", outcome_modes = c(cohort = "skip")
    ),
    "Skipping multi-category outcome"
  )
  expect_equal(unique(skipped$RegressionMatrix$Outcome), "y")
  expect_equal(skipped$Metadata$SkippedOutcomes, "cohort")
})

test_that("MultivariableRegressionTable validates multi-category controls", {
  df <- data.frame(
    cohort = factor(rep(c("A", "B", "C"), each = 15)),
    x1 = rnorm(45)
  )

  expect_error(
    MultivariableRegressionTable(df, "cohort", "x1", reference_levels = c(cohort = "D")),
    "not retained"
  )
  expect_error(
    MultivariableRegressionTable(
      df, "cohort", "x1", outcome_modes = c(cohort = "binary_subset"),
      binary_subsets = list(cohort = c("A", "D"))
    ),
    "exactly two distinct retained levels"
  )
  expect_error(
    MultivariableRegressionTable(df, "cohort", "x1", Method = "lasso", CVFolds = 20),
    "exceeds the smallest class size"
  )
})

test_that("penalized ordinal outcomes require ordinalNet", {
  if (requireNamespace("ordinalNet", quietly = TRUE)) skip("ordinalNet is installed")
  df <- data.frame(
    severity = ordered(rep(c("Mild", "Moderate", "Severe"), each = 20)),
    x1 = rnorm(60)
  )

  expect_error(
    MultivariableRegressionTable(df, "severity", "x1", Method = "lasso", CVFolds = 3),
    "ordinalNet"
  )
})

test_that("MultivariableRegressionTable supports penalized proportional odds", {
  skip_if_not_installed("ordinalNet")
  set.seed(644)
  df <- data.frame(
    severity = ordered(
      sample(c("Mild", "Moderate", "Severe"), 90, replace = TRUE),
      levels = c("Mild", "Moderate", "Severe")
    ),
    x1 = rnorm(90),
    x2 = rnorm(90),
    age = rnorm(90, 50, 8)
  )

  res <- MultivariableRegressionTable(
    df, "severity", c("x1", "x2"), covariates = "age",
    Method = "lasso", CVFolds = 3
  )

  expect_equal(unique(res$RegressionMatrix$OutcomeMode), "ordinal")
  expect_true(all(is.na(res$RegressionMatrix$PValue)))
  expect_equal(res$Metadata$AnalysisSettings$Tuning$severity$Family, "cumulative")
  expect_equal(res$Metadata$AnalysisSettings$Tuning$severity$UnpenalizedTerms, "age")
  expect_true(all(res$Metadata$AnalysisSettings$Tuning$severity$PenaltyFactors["age"] == 0))
})
