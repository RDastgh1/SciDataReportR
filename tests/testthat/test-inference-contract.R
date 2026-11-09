InferenceContractFixture <- function() {
  set.seed(43)
  df_Test <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80), g = factor(rep(c('A', 'B'), 40)), h = factor(rep(c('C', 'D'), each = 40)), time = rexp(80), event = sample(0:1, 80, TRUE))
  obj_Module <- structure(list(datExpr = matrix(df_Test$x, ncol = 1), eigengenes = data.frame(MEblue = df_Test$x)), class = 'feature_wgcna_obj')
  attr(obj_Module$datExpr, 'sample_ids') <- seq_len(nrow(df_Test))
  df_Test$id <- seq_len(nrow(df_Test))
  obj_Correlation <- suppressWarnings(PlotCorrelationsHeatmap(df_Test, predictor_vars = 'x', outcome_vars = 'y'))
  obj_Anova <- suppressWarnings(PlotAnovaRelationshipsMatrix(df_Test, 'g', 'y'))
  calls <- list(
   PlotMiningMatrix = quote(PlotMiningMatrix(df_Test, c('y', 'g'), c('x','h'))),
   PlotCorrelationsHeatmap = quote(PlotCorrelationsHeatmap(df_Test, 'x', 'y')),
   PlotSplitCorrelationHeatmap = quote(PlotSplitCorrelationHeatmap(df_Test, 'x', 'y', 'g')),
   PlotCorrelationComparisons = quote(PlotCorrelationComparisons(df_Test, 'x', 'y', 'g')),
   PlotAnovaRelationshipsMatrix = quote(PlotAnovaRelationshipsMatrix(df_Test, 'g', 'y')),
   PlotChiSqCovar = quote(PlotChiSqCovar(df_Test, 'g', 'h')),
   PlotPhiHeatmap = quote(PlotPhiHeatmap(df_Test, c('g','h'))),
   PlotPointCorrelationsHeatmap = quote(PlotPointCorrelationsHeatmap(df_Test, 'g', 'y')),
   PlotDirectionalHeatmaps = quote(PlotDirectionalHeatmaps(df_Test, c('g','h','x','y'))),
   PlotInteractionEffectsMatrix = quote(PlotInteractionEffectsMatrix(df_Test, 'g', 'y', 'x')),
   PlotNumInteractionEffectsMatrix = quote(PlotNumInteractionEffectsMatrix(df_Test,'x','y', interVar = 'z')),
   PlotCatInteractionEffectsMatrix = quote(PlotCatInteractionEffectsMatrix(df_Test,'x','y', interVar = 'g')),
   PlotInteractionEffectsContinuous = quote(PlotInteractionEffectsContinuous(df_Test, 'g','y','x')),
   MakeUnivariateRegressionTable = quote(MakeUnivariateRegressionTable(df_Test,'y','x')),
   MultivariableRegressionTable = quote(MultivariableRegressionTable(df_Test,'y',c('x','z'))),
   MakeCoxRegressionTable = quote(MakeCoxRegressionTable(df_Test,'time','event','x')),
   MakeComparisonTable = quote(MakeComparisonTable(df_Test,'g',c('x','y'))),
   MakeFacetCatComparisonTable = quote(MakeFacetCatComparisonTable(df_Test,c('g','h'), c('x','y'))),
   CreateStatisticsTable = quote(CreateStatisticsTable(df_Test[,c('g','x','y')],'g')),
   MakePairwiseHeatmap = quote(MakePairwiseHeatmap(df_Test, 'g', c('x','y'), 'A')),
   PlotPairwiseMiningMatrix = quote(PlotPairwiseMiningMatrix(df_Test, 'g', c('x','y'), 'A')),
   Plot2GroupStats = quote(Plot2GroupStats(df_Test,c('x','y'),impClust = 'B',normalClust = 'A',group_var = 'g')),
   PlotSplitViolin = quote(PlotSplitViolin(df_Test, 'y', 'g')),
   PlotPValueComparisons = quote(PlotPValueComparisons(df_Test, 'g', c('x','y'))),
   PlotZScore = quote(PlotZScore(df_Test,'g',c('x','y'))),
   PlotBetaProfile = quote(PlotBetaProfile(df_Test,c('x','z'),'y')),
   PlotVolcanoEffects = quote(PlotVolcanoEffects(df_Test,c('x','z'),'y')),
   EvaluateBiomarkerPerformance = quote(EvaluateBiomarkerPerformance(df_Test, 'g', 'x', CIBootstrapR = 0, Validation = 'none')),
   ScreenBiomarkerPerformance = quote(ScreenBiomarkerPerformance(df_Test, 'g', c('x','z'), CIBootstrapR = 0, Validation = 'none')),
   DiagnosticLikelihoodRatioTable = quote(DiagnosticLikelihoodRatioTable(df_Test, 'g', 'h')),
   PlotAssociations = quote(PlotAssociations(df_Test,'x','y')),
   PlotPartialRegressionScatter = quote(PlotPartialRegressionScatter(df_Test,'x','y','z')),
   PlotBlandAltman = quote(PlotBlandAltman(df_Test,'x','y')),
   AnalyzeFeatureModuleTraits = quote(AnalyzeFeatureModuleTraits(obj_Module, df_Test, 'id', outcome = 'y')),
   calculate_pathway_results = quote(calculate_pathway_results(df_Test, 'y', metabolites = c('x','z'))),
   plotSigAssociations = quote(plotSigAssociations(df_Test,obj_Anova,Pthresh=.99)),
   plotSigCorrelations = quote(plotSigCorrelations(df_Test,obj_Correlation,Pthresh=.99))
  )
  list(Data = df_Test, Calls = calls, Environment = environment())
}

CaptureInferenceWarnings <- function(expr) {
  captured <- new.env(parent = emptyenv())
  captured$Messages <- character()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    captured$Messages <- c(captured$Messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(Value = value, Warnings = captured$Messages)
}

CheckInferenceContract <- function(value, function_name) {
  expect_s3_class(value$Results, "tbl_df")
  expect_type(value$Metadata, "list")
  expect_true(all(c("Outcome", "Predictor", "OutcomeLabel", "PredictorLabel", "Group",
    "Term", "Estimate", "ConfLow", "ConfHigh", "PValue", "PValueFDR", "N",
    "Status", "Reason") %in% names(value$Results)), info = function_name)
  expect_true(all(c("Function", "PackageVersion", "Call", "Arguments", "FormulaTemplate",
    "Models", "ReferenceLevels", "Standardization", "FDRFamily", "ModelCounts", "Failed") %in%
    names(value$Metadata)), info = function_name)
  expect_identical(value$Metadata$Function, function_name)
  expect_identical(value$Metadata$PackageVersion, as.character(utils::packageVersion("SciDataReportR")))
}

test_that("all 37 inferential entrypoints expose numerical results and provenance", {
  fixture <- InferenceContractFixture()
  expect_length(fixture$Calls, 37)
  for (function_name in names(fixture$Calls)) {
    if (function_name == "PlotBlandAltman" && !requireNamespace("BlandAltmanLeh", quietly = TRUE)) next
    if (function_name == "CreateStatisticsTable" && !requireNamespace("arsenal", quietly = TRUE)) next
    captured <- CaptureInferenceWarnings(eval(fixture$Calls[[function_name]], fixture$Environment))
    CheckInferenceContract(captured$Value, function_name)
    expect_false(anyNA(captured$Value$Results$Outcome), info = function_name)
    expect_false(anyNA(captured$Value$Results$Predictor), info = function_name)
    expect_true(all(captured$Value$Results$Status == "ok"), info = function_name)
    expect_false(any(grepl("could not|failed|Insufficient", captured$Warnings)), info = function_name)
  }
})

test_that("canonical linear and logistic estimates retain their model scale", {
  set.seed(910)
  df_Test <- data.frame(x = rnorm(90), cov = rnorm(90), y = rnorm(90))
  df_Test$binary <- factor(rbinom(90, 1, plogis(0.6 * df_Test$x)))
  attr(df_Test$y, "label") <- "Outcome label"
  attr(df_Test$x, "label") <- "Predictor label"
  result <- MakeUnivariateRegressionTable(df_Test, "y", "x", covariates = "cov")
  model <- stats::lm(y ~ x + cov, data = df_Test)
  ci <- stats::confint(model)["x", ]
  expect_equal(result$Results$Estimate, unname(stats::coef(model)["x"]))
  expect_equal(result$Results$ConfLow, unname(ci[1]))
  expect_equal(result$Results$ConfHigh, unname(ci[2]))
  expect_equal(result$Results$PValue, summary(model)$coefficients["x", "Pr(>|t|)"])
  expect_identical(result$Results$OutcomeLabel, "Outcome label")
  expect_identical(result$Results$PredictorLabel, "Predictor label")
  expect_identical(result$Metadata$Arguments$data, "df_Test")
  expect_equal(result$Metadata$ModelCounts$Fit, 1L)

  result <- MultivariableRegressionTable(df_Test, "binary", c("x", "cov"), Standardize = FALSE)
  model <- stats::glm(binary ~ x + cov, data = df_Test, family = stats::binomial())
  expect_equal(result$Results$Estimate, exp(unname(stats::coef(model)[result$Results$Predictor])))
  expect_true(all(result$Results$ReferenceValue == 1))
  expect_equal(result$Metadata$ModelCounts$Fit, 1L)
})

test_that("correlation intervals and pair counts come from the numerical test", {
  set.seed(911)
  df_Test <- data.frame(x = rnorm(45), y = rnorm(45), z = rnorm(45))
  df_Test$x[c(2, 5)] <- NA_real_
  model <- stats::cor.test(df_Test$x, df_Test$y)
  result <- PlotCorrelationsHeatmap(df_Test, "x", "y")
  expect_equal(result$Results$Estimate, unname(model$estimate))
  expect_equal(result$Results$ConfLow, unname(model$conf.int[1]))
  expect_equal(result$Results$ConfHigh, unname(model$conf.int[2]))
  expect_equal(result$Results$N, 43)
  expect_equal(result$Results$PValue, model$p.value)
  expect_equal(result$Results$PValueFDR, model$p.value)
})

test_that("raw p-values control defaults while FDR selection is explicit", {
  set.seed(901)
  df_Test <- as.data.frame(matrix(rnorm(100 * 8), ncol = 8))
  names(df_Test) <- paste0("x", seq_len(8))
  df_Test$y <- 0.2 * df_Test$x1 + rnorm(100)
  vars_Predictors <- paste0("x", seq_len(8))
  raw <- MakeUnivariateRegressionTable(df_Test, "y", vars_Predictors)
  adjusted <- MakeUnivariateRegressionTable(df_Test, "y", vars_Predictors, PVar = "PValueFDR")
  expect_equal(raw$Results$PValueFDR, stats::p.adjust(raw$Results$PValue, "BH"))
  expect_identical(raw$Results$Significant, raw$Results$PValue < 0.05)
  expect_identical(adjusted$Results$Significant, adjusted$Results$PValueFDR < 0.05)
  expect_true(any(raw$Results$Significant != adjusted$Results$Significant))
  expect_equal(raw$Results$PValue, adjusted$Results$PValue)
  expect_equal(raw$Results$Estimate, adjusted$Results$Estimate)
  expect_false(formals(MakePairwiseHeatmap)$adjusted_outline)
  expect_false(formals(PlotPairwiseMiningMatrix)$adjusted_outline)
})

test_that("standardization metadata records actual SDs and primary model IDs", {
  set.seed(912)
  df_Test <- data.frame(x = rnorm(60, 10, 3), y = rnorm(60, 20, 5), z = rnorm(60))
  result <- MakeUnivariateRegressionTable(df_Test, "y", "x", Standardize = TRUE)
  df_Scale <- result$Metadata$Standardization
  expect_equal(df_Scale$SD[match(c("y", "x"), df_Scale$Variable)], c(stats::sd(df_Test$y), stats::sd(df_Test$x)))
  expect_true(all(df_Scale$ModelID %in% result$Metadata$Models$ModelID))
  result <- MultivariableRegressionTable(df_Test, "y", c("x", "z"))
  expect_equal(result$Metadata$ModelCounts$Fit, 1L)
  expect_true(all(result$Metadata$Standardization$ModelID %in% result$Metadata$Models$ModelID))
})

test_that("all interaction terms have intervals while the matrix remains available", {
  set.seed(913)
  df_Test <- data.frame(y = rnorm(90), x = rnorm(90), g = factor(rep(c("A", "B", "C"), 30)))
  result <- PlotInteractionEffectsMatrix(df_Test, "g", "y", "x")
  model <- stats::lm(y ~ x * g, data = df_Test)
  vars_Terms <- grep(":", names(stats::coef(model)), value = TRUE)
  expect_setequal(result$Results$Term, vars_Terms)
  expect_equal(result$Results$Estimate, unname(stats::coef(model)[result$Results$Term]))
  expect_true(all(!is.na(result$Results$ConfLow)))
  expect_true(all(!is.na(result$Results$ConfHigh)))
  expect_equal(result$Metadata$ModelCounts$Fit, 1L)
  expect_true(all(c("Unadjusted", "FDRCorrected", "raw_data") %in% names(result)))
})

test_that("pairwise Results keep raw p-values separate from legacy corrections", {
  set.seed(914)
  df_Test <- data.frame(y = rnorm(90), g = factor(rep(c("A", "B", "C"), 30)))
  result <- MakeComparisonTable(df_Test, "g", "y", AddPairwise = TRUE, PairwiseMethod = "bonferroni")
  df_Pairwise <- result$Results[result$Results$TestRole == "Pairwise", ]
  expect_equal(nrow(df_Pairwise), 3)
  expected <- stats::t.test(df_Test$y[df_Test$g == "A"], df_Test$y[df_Test$g == "B"])
  expect_equal(df_Pairwise$PValue[df_Pairwise$Term == "A - B"], expected$p.value)
  expect_equal(df_Pairwise$PValueFDR, stats::p.adjust(df_Pairwise$PValue, "BH"))
  expect_equal(result$Metadata$ModelCounts$Fit, 4L)
  expect_s3_class(result, "gtsummary")
  expect_true(any(grepl("^pw_", names(result$table_body))))
})

test_that("failed requests stay visible and independent requests continue", {
  set.seed(915)
  df_Test <- data.frame(x = rnorm(40), y = rnorm(40), constant = 1, missing = NA_real_)
  captured <- CaptureInferenceWarnings(MakeUnivariateRegressionTable(df_Test, "y", c("x", "constant", "missing")))
  result <- captured$Value
  expect_setequal(result$Results$Predictor, c("x", "constant", "missing"))
  expect_equal(sum(result$Results$Status == "ok"), 1)
  expect_equal(nrow(result$Metadata$Failed), 2)
  expect_true(all(c("Outcome", "Predictor", "Group", "Stage", "Reason") %in% names(result$Metadata$Failed)))
  expect_true(any(grepl("y.*constant", captured$Warnings)))
  expect_true(any(grepl("y.*missing", captured$Warnings)))
  for (call in list(
    quote(PlotCorrelationsHeatmap(df_Test, "constant", "y")),
    quote(MakeUnivariateRegressionTable(df_Test, "y", "missing")),
    quote(PlotBetaProfile(df_Test, "constant", "y")))) {
    captured <- CaptureInferenceWarnings(eval(call))
    expect_true(all(captured$Value$Results$Status == "failed"))
    expect_true(all(is.na(captured$Value$Results$Estimate)))
    expect_true(nrow(captured$Value$Metadata$Failed) > 0L)
    expect_true(length(captured$Warnings) > 0L)
  }
  expect_error(MakeUnivariateRegressionTable(df_Test, "absent", "x"), "absent")
})

test_that("plot and table classes, composition, indexing, and legacy elements survive", {
  fixture <- InferenceContractFixture()
  result <- eval(fixture$Calls$PlotSplitViolin, fixture$Environment)
  expect_s3_class(result, "ggplot")
  modified <- result + ggplot2::theme_bw()
  expect_identical(modified$Results, result$Results)
  expect_no_error(ggplot2::ggplot_build(modified))
  result <- eval(fixture$Calls$MakeComparisonTable, fixture$Environment)
  expect_s3_class(result, "gtsummary")
  expect_no_error(gtsummary::as_gt(result))
  result <- eval(fixture$Calls$calculate_pathway_results, fixture$Environment)
  expect_s3_class(result, "data.frame")
  expect_true("Metabolite" %in% names(result))
  expect_false("Results" %in% names(result))
  result <- eval(fixture$Calls$plotSigCorrelations, fixture$Environment)
  expect_length(result, nrow(result$Results))
  if (length(result)) expect_s3_class(result[[1]], "ggplot")
  result <- eval(fixture$Calls$MakeUnivariateRegressionTable, fixture$Environment)
  expect_no_error(ggplot2::ggplot_build(PlotForestFromTable(result$Results)))
})

test_that("empty requests return a typed schema and diagnostic boundaries remain valid", {
  df_Test <- data.frame(x = seq_len(10), y = seq_len(10))
  result <- PlotMiningMatrix(df_Test, character(), character())
  CheckInferenceContract(result, "PlotMiningMatrix")
  expect_equal(nrow(result$Results), 0L)
  expect_equal(result$Metadata$ModelCounts$Attempted, 0L)
  df_Diagnostic <- data.frame(outcome = factor(c("No", "No", "Yes", "Yes")), test = factor(c("Negative", "Negative", "Positive", "Positive")))
  result <- DiagnosticLikelihoodRatioTable(df_Diagnostic, "outcome", "test", positive_level = "Yes")
  expect_true(any(result$Results$Estimate == 0))
  expect_true(any(is.infinite(result$Results$Estimate)))
  expect_true(all(is.na(result$Results$PValue)))
  expect_equal(nrow(result$Metadata$Failed), 0L)
})

test_that("every estimation family returns an inspectable all-failed analysis", {
  fixture <- InferenceContractFixture()
  df_Test <- fixture$Data
  df_Test$y <- NA_real_
  df_Test$x <- 1
  df_Test$z <- 1
  df_Test$h <- factor(rep(c("C", "D"), 40))
  df_Test$g[seq_len(40)] <- NA
  df_Test$h[41:80] <- NA
  assign("df_Test", df_Test, envir = fixture$Environment)
  for (function_name in setdiff(names(fixture$Calls), c("plotSigAssociations", "plotSigCorrelations"))) {
    if (function_name == "PlotBlandAltman" && !requireNamespace("BlandAltmanLeh", quietly = TRUE)) next
    if (function_name == "CreateStatisticsTable" && !requireNamespace("arsenal", quietly = TRUE)) next
    captured <- CaptureInferenceWarnings(eval(fixture$Calls[[function_name]], fixture$Environment))
    CheckInferenceContract(captured$Value, function_name)
    expect_gt(nrow(captured$Value$Results), 0L)
    expect_true(all(captured$Value$Results$Status == "failed"), info = function_name)
    expect_true(all(is.na(captured$Value$Results$Estimate)), info = function_name)
    expect_gt(nrow(captured$Value$Metadata$Failed), 0L)
    expect_true(any(grepl("outcome.*predictor", captured$Warnings)), info = function_name)
  }
  obj_Correlation <- eval(fixture$Calls$PlotCorrelationsHeatmap, InferenceContractFixture()$Environment)
  obj_Anova <- eval(fixture$Calls$PlotAnovaRelationshipsMatrix, InferenceContractFixture()$Environment)
  for (call in list(quote(plotSigCorrelations(df_Test, obj_Correlation, Pthresh = 1e-100)),
                    quote(plotSigAssociations(df_Test, obj_Anova, Pthresh = 1e-100)))) {
    result <- eval(call)
    expect_length(result, 0L)
    expect_equal(nrow(result$Results), 0L)
    expect_equal(nrow(result$Metadata$Failed), 0L)
  }
})

test_that("factor regression results retain every coefficient and model count", {
  set.seed(916)
  df_Test <- data.frame(y = rnorm(90), g = factor(rep(c("A", "B", "C"), 30)), x = rnorm(90),
    outcome = factor(sample(c("Control", "Mild", "Severe"), 90, replace = TRUE)))
  result <- MultivariableRegressionTable(df_Test, "y", c("g", "x"), Standardize = FALSE)
  model <- stats::lm(y ~ g + x, data = df_Test)
  expect_setequal(result$Results$Term, c("gB", "gC", "x"))
  expect_equal(result$Results$Estimate, unname(stats::coef(model)[result$Results$Term]))
  expect_equal(result$Results$ConfLow, unname(stats::confint(model)[result$Results$Term, 1]))
  expect_equal(result$Results$PValueFDR, stats::p.adjust(result$Results$PValue, "BH"))
  expect_equal(result$Metadata$ModelCounts$Fit, 1)
  nominal <- MultivariableRegressionTable(df_Test, "outcome", c("g", "x"), outcome_modes = "multinomial")
  expect_equal(nrow(nominal$Results), 6)
  expect_equal(nominal$Metadata$ModelCounts$Fit, 1)
  ordinal <- MultivariableRegressionTable(df_Test, "outcome", c("g", "x"), outcome_modes = "ordinal")
  expect_equal(nrow(ordinal$Results), 3)
  expect_equal(ordinal$Metadata$ModelCounts$Fit, 1)
  split <- MultivariableRegressionTable(df_Test, "outcome", c("g", "x"), outcome_modes = "one_vs_rest")
  expect_equal(nrow(split$Results), 9)
  expect_equal(split$Metadata$ModelCounts$Fit, 3)
  expect_true(all(split$Metadata$Standardization$ModelID %in% split$Metadata$Models$ModelID))
})

test_that("categorical and mixed association plots expose library statistics", {
  fixture <- InferenceContractFixture()
  for (variables in list(c("g", "h"), c("g", "x"))) {
    result <- PlotAssociations(fixture$Data, variables[1], variables[2])
    expect_true(all(result$Results$Status == "ok"))
    expect_true(all(is.finite(result$Results$PValue)))
    expect_equal(nrow(result$Metadata$Failed), 0)
  }
})

test_that("entirely missing categorical data remain failed analysis requests", {
  set.seed(917)
  df_Test <- data.frame(y = rnorm(40), x = rnorm(40),
    g = factor(rep(c("A", "B"), 20)), h = factor(rep(c("C", "D"), each = 20)))
  df_Test$g[] <- NA
  calls <- list(
    quote(PlotAnovaRelationshipsMatrix(df_Test, "g", "y")),
    quote(PlotPointCorrelationsHeatmap(df_Test, "g", "y")),
    quote(PlotPhiHeatmap(df_Test, c("g", "h"))),
    quote(MakeComparisonTable(df_Test, "g", "y")),
    quote(MakeUnivariateRegressionTable(df_Test, "g", "x")),
    quote(MultivariableRegressionTable(df_Test, "g", "x")),
    quote(DiagnosticLikelihoodRatioTable(df_Test, "g", "h"))
  )
  for (call in calls) {
    captured <- CaptureInferenceWarnings(eval(call))
    expect_true(all(captured$Value$Results$Status == "failed"))
    expect_true(all(captured$Value$Results$N == 0))
    expect_gt(nrow(captured$Value$Metadata$Failed), 0)
    expect_true(any(grepl("outcome.*predictor", captured$Warnings)))
  }
  # A standalone categorical validator retains its strict input safeguards.
  expect_error(ScidrPrepareCategorical(df_Test, "g", "group_var"), "no observed")
})

test_that("missing grouping columns remain inspectable across grouped entrypoints", {
  fixture <- InferenceContractFixture()
  df_Test <- fixture$Data
  df_Test$g[] <- NA
  assign("df_Test", df_Test, envir = fixture$Environment)
  for (function_name in c("PlotSplitCorrelationHeatmap", "PlotCorrelationComparisons", "MakePairwiseHeatmap",
    "PlotPairwiseMiningMatrix", "Plot2GroupStats", "PlotPValueComparisons", "plotSigAssociations")) {
    captured <- CaptureInferenceWarnings(eval(fixture$Calls[[function_name]], fixture$Environment))
    CheckInferenceContract(captured$Value, function_name)
    expect_gt(nrow(captured$Value$Metadata$Failed), 0L)
    expect_true(any(grepl("outcome.*predictor", captured$Warnings)), info = function_name)
  }
})

test_that("unadjusted correlation fallbacks warn and record the failed adjustment", {
  set.seed(918)
  df_Test <- data.frame(x = rnorm(30), y = rnorm(30), cov = factor(rep("Only level", 30)))
  captured <- CaptureInferenceWarnings(PlotCorrelationsHeatmap(df_Test, "x", "y", covariates = "cov"))
  expect_true(any(grepl("y.*x.*unadjusted", captured$Warnings)))
  expect_equal(captured$Value$Metadata$Failed$Stage, "adjustment")
  expect_identical(captured$Value$Results$AdjustmentStatus, "unadjusted_fallback")
  expect_identical(captured$Value$Results$Status, "failed")
  expect_true(is.na(captured$Value$Results$PValue))
  expect_equal(captured$Value$Metadata$ModelCounts$Failed, 1L)
  expect_equal(captured$Value$Unadjusted$r[1, 1], stats::cor(df_Test$x, df_Test$y))
})

test_that("wrappers aggregate unavailable adjustments without duplicating warnings", {
  set.seed(919)
  df_Test <- data.frame(x = rnorm(40), y = rnorm(40),
    group = factor(rep(c("A", "B"), 20)), cov = factor(rep("Only level", 40)))
  captured <- CaptureInferenceWarnings(PlotSplitCorrelationHeatmap(df_Test, "x", "y", "group", covariates = "cov"))
  expect_equal(length(captured$Warnings), 2L)
  expect_true(all(captured$Value$Results$Status == "failed"))
  expect_true(all(is.na(captured$Value$Results$Estimate)))
  expect_setequal(captured$Value$Metadata$Failed$Group, c("A", "B"))
  expect_true(all(captured$Value$Metadata$Failed$Stage == "adjustment"))
})

test_that("pairwise primary stars use BH independently of legacy corrections", {
  set.seed(921)
  df_Test <- as.data.frame(matrix(rnorm(100 * 12), ncol = 12))
  df_Test$g <- factor(rep(c("A", "B"), each = 50))
  df_Test$V1 <- df_Test$V1 + 0.4 * (df_Test$g == "B")
  vars_Outcomes <- paste0("V", seq_len(12))
  for (Compute in list(MakePairwiseHeatmap, PlotPairwiseMiningMatrix)) {
    raw <- Compute(df_Test, "g", vars_Outcomes, "A", p_adjust_method = "none")
    fdr <- Compute(df_Test, "g", vars_Outcomes, "A", p_adjust_method = "none", PVar = "PValueFDR")
    expect_equal(fdr$Results$PValueFDR, stats::p.adjust(fdr$Results$PValue, "BH"))
    expect_equal(fdr$Results$AdjustedPValue, fdr$Results$PValue)
    expect_true(any(nzchar(raw$Results$SignificanceLabel)))
    expect_true(all(!nzchar(fdr$Results$SignificanceLabel)))
    plot <- if (!is.null(fdr$Plot)) fdr$Plot else fdr$Plots$Continuous
    expect_true(all(!nzchar(plot$data$SignificanceLabel)))
    none <- Compute(df_Test, "g", vars_Outcomes, "A", star_p = "none")
    expect_true(all(!nzchar(none$Results$SignificanceLabel)))
    legacy <- Compute(df_Test, "g", vars_Outcomes, "A", star_p = "adjusted", p_adjust_method = "holm")
    expect_equal(legacy$Metadata$Arguments$PVar, "AdjustedPValue")
    if (identical(Compute, MakePairwiseHeatmap)) {
      expect_equal(raw$Results$ConfLow, raw$Results$ConfidenceIntervalLower)
      expect_equal(raw$Results$ConfHigh, raw$Results$ConfidenceIntervalUpper)
    }
  }
  raw <- MakeComparisonTable(df_Test, "g", vars_Outcomes)
  fdr <- MakeComparisonTable(df_Test, "g", vars_Outcomes, PVar = "PValueFDR")
  strict <- MakeComparisonTable(df_Test, "g", vars_Outcomes, Pthresh = 0.001)
  expect_match(gtsummary::as_tibble(raw)[["**p-value**"]][1], "__")
  expect_false(grepl("__", gtsummary::as_tibble(fdr)[["**p-value**"]][1]))
  expect_false(grepl("__", gtsummary::as_tibble(strict)[["**p-value**"]][1]))
})

test_that("one-sided pairwise results retain the requested raw test", {
  set.seed(922)
  df_Test <- data.frame(y = rnorm(80), g = factor(rep(c("A", "B"), each = 40)))
  df_Test$y <- df_Test$y + 0.5 * (df_Test$g == "B")
  for (parametric in c(TRUE, FALSE)) {
    result <- MakeComparisonTable(df_Test, "g", "y", AddPairwise = TRUE,
      alternative = "greater", Parametric = parametric)
    expected <- if (parametric) stats::t.test(df_Test$y[df_Test$g == "B"], df_Test$y[df_Test$g == "A"], alternative = "greater") else stats::wilcox.test(df_Test$y[df_Test$g == "B"], df_Test$y[df_Test$g == "A"], alternative = "greater", exact = FALSE)
    df_Pairs <- result$Results[result$Results$TestRole == "Pairwise", ]
    expect_equal(nrow(df_Pairs), 1)
    expect_equal(df_Pairs$PValue, expected$p.value)
    expect_equal(df_Pairs$PValueFDR, expected$p.value)
    expect_equal(df_Pairs$N, nrow(df_Test))
  }
})
