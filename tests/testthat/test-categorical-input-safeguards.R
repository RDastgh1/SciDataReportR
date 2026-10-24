PrepareCategoricalForTest <- function(data, variables = "group", arg_name = "group_var") {
  ScidrPrepareCategorical(data, variables, arg_name)
}

CaptureCategoricalWarnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, messages = messages)
}

test_that("numeric category preparation is bounded and preserves input metadata", {
  for (x in list(c(0L, 1L, NA_integer_), c(1, 0, NA), 1:10)) {
    df_Test <- data.frame(group = x)
    attr(df_Test$group, "label") <- "Study cohort"
    df_Original <- df_Test
    result <- CaptureCategoricalWarnings(PrepareCategoricalForTest(df_Test))
    expect_length(result$messages, 1)
    expect_match(result$messages, "converting locally to a factor")
    expect_match(result$messages, "Few unique values")
    expect_true(is.factor(result$value$group))
    expect_equal(levels(result$value$group), as.character(sort(unique(x[!is.na(x)]))))
    expect_equal(attr(result$value$group, "label"), "Study cohort")
    expect_equal(is.na(result$value$group), is.na(x))
    expect_identical(df_Test, df_Original)
  }
  for (x in list(1:11, 1:1000, c(0.1, 0.2), c(0, Inf), c(-Inf, 1))) {
    expect_error(PrepareCategoricalForTest(data.frame(group = x)),
      "explicitly convert `group` to a factor")
  }
  expect_error(PrepareCategoricalForTest(data.frame(group = 1:1000)),
    "1000 observed non-missing values")
  expect_warning(PrepareCategoricalForTest(data.frame(group = c(1, NA, NaN))),
    "1 observed non-missing")
})

test_that("declared categories preserve ordering and count only observed levels", {
  df_Test <- data.frame(group = factor(c("b", "a", NA), levels = c("b", "a", letters[3:26]), ordered = TRUE))
  attr(df_Test$group, "label") <- "Ordered category"
  expect_identical(PrepareCategoricalForTest(df_Test), df_Test)
  for (x in list(letters[1:11], factor(letters[1:11]))) {
    df_Test <- data.frame(group = x)
    result <- CaptureCategoricalWarnings(PrepareCategoricalForTest(df_Test))
    expect_length(result$messages, 1)
    expect_match(result$messages, "11 observed categories")
    expect_identical(result$value, df_Test)
  }
  df_Test <- data.frame(group = c(TRUE, FALSE, NA))
  expect_identical(PrepareCategoricalForTest(df_Test), df_Test)
  expect_silent(PrepareCategoricalForTest(data.frame(group = factor(c(0.1, 0.2)))))
  expect_silent(PrepareCategoricalForTest(data.frame(group = factor(1:10))))
})

test_that("empty, unsupported, and misspelled category inputs fail clearly", {
  for (x in list(c(NA_real_, NA_real_), factor(c(NA, NA)), character())) {
    expect_error(PrepareCategoricalForTest(data.frame(group = x)), "no observed non-missing")
  }
  for (x in list(as.Date("2026-01-01"), as.POSIXct("2026-01-01"), 1 + 1i, I(list(1, 2)), I(matrix(1:4, nrow = 2)))) {
    df_Test <- data.frame(row = seq_along(x))
    if (!is.null(dim(x))) df_Test <- data.frame(row = seq_len(nrow(x)))
    df_Test$group <- x
    expect_error(PrepareCategoricalForTest(df_Test), "unsupported column type")
  }
  expect_error(PrepareCategoricalForTest(data.frame(group = 1), "typo"), "typo.*group_var")
})

test_that("warning state spans nested calls but does not leak across public calls", {
  Nested <- function(data) ScidrPrepareCategorical(data, "group", "nested_group")
  Outer <- function(data) {
    df_Prepared <- ScidrPrepareCategorical(data, c("group", "other"), "group_var")
    Nested(df_Prepared)
  }
  df_Test <- data.frame(group = factor(letters[1:11]), other = factor(letters[1:11]))
  expect_length(CaptureCategoricalWarnings(Outer(df_Test))$messages, 2)
  expect_length(CaptureCategoricalWarnings(Outer(df_Test))$messages, 2)
  expect_length(CaptureCategoricalWarnings(PrepareCategoricalForTest(data.frame(group = 0:1)))$messages, 1)
  expect_length(CaptureCategoricalWarnings(PrepareCategoricalForTest(data.frame(group = 0:1)))$messages, 1)
})

test_that("binary numeric Fill builds rainclouds and leaves data unchanged", {
  skip_if_not_installed("ggrain")
  set.seed(62)
  df_Test <- data.frame(marker = rnorm(120), cohort = rep(c(0, 1, NA), each = 40))
  attr(df_Test$cohort, "label") <- "Study cohort"
  df_Original <- df_Test
  result <- CaptureCategoricalWarnings(PlotContinuousDistributions(df_Test, "marker", Fill = "cohort"))
  expect_length(result$messages, 1)
  expect_s3_class(result$value, "ggplot")
  expect_true(is.factor(result$value$data$cohort))
  expect_equal(result$value$scales$get_scales("fill")$name, "Study cohort")
  built <- CaptureCategoricalWarnings(ggplot2::ggplot_build(result$value))
  expect_false(any(grepl("Continuous value supplied|Ignoring unknown|dropped.*group", built$messages)))
  expect_identical(df_Test, df_Original)
  automatic <- suppressWarnings(PlotContinuousDistributions(df_Test, Fill = "cohort"))
  expect_length(levels(automatic$data$name), 1)
  expect_match(levels(automatic$data$name), "^marker")
  expect_error(suppressWarnings(PlotContinuousDistributions(df_Test, c("marker", "cohort"), Fill = "cohort")), "must not also")
  expect_error(PlotContinuousDistributions(df_Test, "marker", Fill = c("cohort", "marker")), "exactly one")
})

test_that("categorical public entry points reject accidental continuous groups early", {
  df_Test <- data.frame(group = seq_len(1000), marker = seq_len(1000) / 7,
    other = rep(c("a", "b"), 500), count = rep(1:2, 500))
  Calls <- list(
    PlotContinuousDistributions = function() PlotContinuousDistributions(df_Test, "marker", Fill = "group"),
    PlotSplitViolin = function() PlotSplitViolin(df_Test, "marker", "group"),
    PlotSpiderChart = function() PlotSpiderChart(df_Test, "marker", group_var = "group"),
    PlotClusterBoxplot = function() PlotClusterBoxplot(df_Test, ClusterVar = "group", variables = "marker"),
    PlotClusterPrevalence = function() PlotClusterPrevalence(df_Test, "group", "other"),
    PlotSplitCorrelationHeatmap = function() PlotSplitCorrelationHeatmap(df_Test, "marker", "marker", split_var = "group"),
    PlotCorrelationComparisons = function() PlotCorrelationComparisons(df_Test, "marker", "marker", group_var = "group"),
    Plot2GroupStats = function() Plot2GroupStats(df_Test, variables = "marker", group_var = "group"),
    PlotPValueComparisons = function() PlotPValueComparisons(df_Test, "group", "marker"),
    PlotZScore = function() PlotZScore(df_Test, "group", "marker"),
    PlotPairwiseMiningMatrix = function() PlotPairwiseMiningMatrix(df_Test, "group", "marker", Referent = "1"),
    MakePairwiseHeatmap = function() MakePairwiseHeatmap(df_Test, "group", "marker", Referent = "1"),
    MakeComparisonTable = function() MakeComparisonTable(df_Test, "group", "marker"),
    CreateStatisticsTable = function() CreateStatisticsTable(df_Test, "group"),
    CreateNormativeTScoreModel = function() CreateNormativeTScoreModel(df_Test, "marker", "count", "marker", "group", 1),
    PlotCategoricalDistributions = function() PlotCategoricalDistributions(df_Test, "group"),
    PlotAnovaRelationshipsMatrix = function() PlotAnovaRelationshipsMatrix(df_Test, "group", "marker"),
    PlotChiSqCovarPredictor = function() PlotChiSqCovar(df_Test, "group", "other"),
    PlotChiSqCovarOutcome = function() PlotChiSqCovar(df_Test, "other", "group"),
    PlotChiSqCovarStratum = function() PlotChiSqCovar(df_Test, "other", "other", covariates = "group"),
    DiagnosticLikelihoodRatioTablePredictor = function() DiagnosticLikelihoodRatioTable(df_Test, "other", "group"),
    DiagnosticLikelihoodRatioTableStratum = function() DiagnosticLikelihoodRatioTable(df_Test, "other", "other", stratify_by = "group")
  )
  for (name in names(Calls)) {
    expect_error(Calls[[name]](), "explicitly convert `group` to a factor", info = name)
  }
})

test_that("group conversion retains binary requirements and explicit factor override", {
  skip_if_not_installed("emmeans")
  df_Test <- data.frame(group = rep(1:3, each = 20), marker = seq_len(60) / 7)
  expect_error(suppressWarnings(PlotSplitViolin(df_Test, "marker", "group")), "<=2 levels")
  expect_error(suppressWarnings(PlotSplitCorrelationHeatmap(df_Test, "marker", "marker", "group")), "exactly two observed")
  df_Test$group <- factor(rep(c(0.1, 0.2, 0.3), each = 20))
  expect_s3_class(PlotCategoricalDistributions(df_Test, "group"), "ggplot")
})

test_that("analysis group preparation agrees with explicit factors", {
  df_Test <- data.frame(group = rep(0:1, each = 30), marker = seq_len(60) / 7)
  result <- CaptureCategoricalWarnings(MakeComparisonTable(df_Test, "group", "marker"))
  expect_length(result$messages, 1)
  df_Test$group <- factor(df_Test$group)
  explicit <- MakeComparisonTable(df_Test, "group", "marker")
  expect_equal(result$value$table_body, explicit$table_body)
})

test_that("numeric diagnostic codes support multiple categories without changing outcomes", {
  skip_if_not_installed("gt")
  df_Test <- data.frame(outcome = rep(c(FALSE, TRUE), each = 30),
    code = rep(1:3, 20), stratum = rep(0:1, 30))
  result <- CaptureCategoricalWarnings(DiagnosticLikelihoodRatioTable(df_Test, "outcome", "code", stratify_by = "stratum"))
  expect_length(result$messages, 2)
  expect_equal(nrow(result$value$Results), 6)
  expect_equal(result$value$Metadata$OutcomeLevels$outcome$positive, "TRUE")
})

test_that("labelled numeric codes retain their variable label", {
  df_Test <- data.frame(group = sjlabelled::set_label(c(0, 1, NA), "Diagnosis cohort"))
  df_Original <- df_Test
  result <- CaptureCategoricalWarnings(PrepareCategoricalForTest(df_Test))
  expect_length(result$messages, 1)
  expect_true(is.factor(result$value$group))
  expect_equal(sjlabelled::get_label(result$value$group), "Diagnosis cohort")
  expect_identical(df_Test, df_Original)
})

test_that("faceted comparisons retain numeric-facet exclusion and warn once", {
  df_Test <- data.frame(group = factor(rep(letters[1:11], each = 3)),
    numeric_facet = rep(1:3, 11), marker = seq_len(33))
  Seen <- new.env(parent = emptyenv())
  Comparison <- function(data, group_var, variables, ...) {
    Seen$group <- group_var
    df_Prepared <- ScidrPrepareCategorical(data, group_var, "group_var")
    gtsummary::tbl_summary(df_Prepared[variables])
  }
  result <- CaptureCategoricalWarnings(MakeFacetCatComparisonTable(
    df_Test, FacetVariables = c("group", "numeric_facet"),
    variables = "marker", CompFun = Comparison
  ))
  expect_length(result$messages, 1)
  expect_equal(Seen$group, "group")
  expect_s3_class(result$value, "tbl_merge")
})

test_that("normative plots prepare reference groups without changing model inputs", {
  skip_if_not_installed("broom")
  set.seed(25)
  df_Test <- data.frame(score = rnorm(80, 50, 4), count = 1,
    age = runif(80, 20, 70), group = rep(0:1, each = 40))
  df_Original <- df_Test
  without_plots <- CreateNormativeTScoreModel(df_Test, "score", "count", "age",
    "group", 0, log_transform = FALSE, return_plots = FALSE)
  result <- CaptureCategoricalWarnings(CreateNormativeTScoreModel(
    df_Test, "score", "count", "age", "group", 0,
    log_transform = FALSE, return_plots = TRUE
  ))
  expect_length(result$messages, 1)
  expect_identical(result$value$data, without_plots$data)
  expect_equal(stats::coef(result$value$model), stats::coef(without_plots$model))
  expect_true(is.numeric(result$value$data$group))
  expect_true(is.factor(result$value$plots$raw$data$group))
  expect_identical(df_Test, df_Original)
})

test_that("PCA continuous colours and explicit categorical authorization coexist", {
  skip_if_not_installed("plotly")
  set.seed(90)
  df_Test <- as.data.frame(matrix(rnorm(240), nrow = 60))
  df_Test$colour <- seq(0.1, 6, by = 0.1)
  pca <- suppressWarnings(CreatePCAObject(df_Test, VarsToReduce = paste0("V", 1:4), numComponents = 2))
  continuous <- CaptureCategoricalWarnings(plotPCA(pca, Var = "colour", ColorType = "continuous", Mode = "2D", Components = 1:2))
  expect_false(any(grepl("categorical|observed categories|converting locally", continuous$messages)))
  expect_s3_class(continuous$value, "plotly")
  explicit <- CaptureCategoricalWarnings(plotPCA(pca, Var = "colour", ColorType = "factor", Mode = "2D", Components = 1:2))
  expect_equal(sum(grepl("60 observed categories", explicit$messages)), 1)
  expect_s3_class(explicit$value, "plotly")
  expect_no_error(suppressWarnings(plotly::plotly_build(explicit$value)))
})
