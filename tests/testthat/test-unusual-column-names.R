NameSafetyFixture <- function(marker = NULL) {
  set.seed(1874)
  df_Test <- data.frame(x = rnorm(96), y = rnorm(96), z = rnorm(96),
    w = rnorm(96), g = factor(rep(c("A", "B"), 48)),
    h = factor(rep(c("C", "D", "E"), each = 32)),
    count = rep(1:3, 32), id = seq_len(96))
  df_Test$y <- 5 + df_Test$x + 0.4 * df_Test$z + rnorm(96)
  df_Test$z[2] <- NA_real_
  vars_Names <- names(df_Test)
  for (variable in vars_Names) {
    df_Test[[variable]] <- sjlabelled::set_label(df_Test[[variable]], paste("Label", variable))
  }
  vars_Map <- stats::setNames(vars_Names, vars_Names)
  if (!is.null(marker)) {
    vars_Map[] <- paste0(marker, vars_Names)
    names(df_Test) <- unname(vars_Map)
  }
  list(Data = df_Test, Names = vars_Map)
}

NameSafetyNumbers <- function(df_Results) {
  df_Results <- as.data.frame(df_Results)
  df_Results[vapply(df_Results, is.numeric, logical(1))]
}

CheckNameSafetyInference <- function(result, reference) {
  df_Results <- if (is.data.frame(result)) attr(result, "Results") else result$Results
  df_Reference <- if (is.data.frame(reference)) attr(reference, "Results") else reference$Results
  expect_gt(nrow(df_Results), 0)
  if ("Status" %in% names(df_Results)) expect_true(all(df_Results$Status == "ok"))
  expect_equal(NameSafetyNumbers(df_Results), NameSafetyNumbers(df_Reference), tolerance = 1e-8)
  expect_true(any(is.finite(df_Results$PValue)))
}

CheckNameSafetyPlots <- function(result) {
  if (inherits(result, "ggplot")) {
    expect_s3_class(suppressMessages(ggplot2::ggplot_build(result)), "ggplot_built")
  } else if (is.list(result)) {
    for (field in intersect(names(result), c("Plot", "plot", "Plots", "plots", "Unadjusted", "FDRCorrected", "Continuous", "Categorical", "p", "p_FDR"))) {
      CheckNameSafetyPlots(result[[field]])
    }
  }
}

test_that("shared formula helpers preserve exact column names", {
  vars_Names <- c("% marker", "marker value", "marker+value", "marker`value", "marker\\value")
  for (variable in vars_Names) {
    expect_identical(all.vars(ScidrRegressionFormula(variable, outcome = variable)), variable)
    expect_identical(all.vars(stats::as.formula(paste("~", ScidrQuoteFormulaNames(variable)))), variable)
    expect_identical(all.vars(.SdrFormula("outcome", variable)), c("outcome", variable))
    expect_identical(all.vars(.BiomarkerFormula("outcome", variable)), c("outcome", variable))
  }
})

# All roles are renamed together, with labels unchanged, to expose downstream lookup bugs.
name_safety_calls <- list(
  beta = function(df_Data, vars_Names, vars_Covariates) PlotBetaProfile(df_Data, vars_Names[["x"]], vars_Names[["y"]], covariates = vars_Covariates, InteractiveLabels = FALSE),
  volcano_continuous = function(df_Data, vars_Names, vars_Covariates) PlotVolcanoEffects(df_Data, vars_Names[["x"]], vars_Names[["y"]], covariates = vars_Covariates, InteractiveLabels = FALSE),
  volcano_categorical = function(df_Data, vars_Names, vars_Covariates) PlotVolcanoEffects(df_Data, vars_Names[["x"]], vars_Names[["g"]], covariates = vars_Covariates, InteractiveLabels = FALSE),
  anova = function(df_Data, vars_Names, vars_Covariates) PlotAnovaRelationshipsMatrix(df_Data, vars_Names[["g"]], vars_Names[["y"]], covariates = vars_Covariates),
  mining = function(df_Data, vars_Names, vars_Covariates) PlotPairwiseMiningMatrix(df_Data, vars_Names[["g"]], vars_Names[["y"]], "A", covariates = vars_Covariates),
  violin_parametric = function(df_Data, vars_Names, vars_Covariates) do.call(PlotSplitViolin, list(data = df_Data, Var = vars_Names[["y"]], group_var = vars_Names[["g"]], covariates = vars_Covariates)),
  violin_nonparametric = function(df_Data, vars_Names, vars_Covariates) do.call(PlotSplitViolin, list(data = df_Data, Var = vars_Names[["y"]], group_var = vars_Names[["g"]], covariates = vars_Covariates, nonparametric = TRUE)),
  comparison = function(df_Data, vars_Names, vars_Covariates) MakeComparisonTable(df_Data, vars_Names[["g"]], vars_Names[["y"]], covariates = vars_Covariates),
  heatmap = function(df_Data, vars_Names, vars_Covariates) MakePairwiseHeatmap(df_Data, vars_Names[["g"]], vars_Names[["y"]], "A", covariates = vars_Covariates),
  interaction = function(df_Data, vars_Names, vars_Covariates) PlotInteractionEffectsMatrix(df_Data, vars_Names[["g"]], vars_Names[["y"]], vars_Names[["x"]], covariates = vars_Covariates),
  interaction_covariates = function(df_Data, vars_Names, vars_Covariates) PlotInteractionEffectsMatrix(df_Data, vars_Names[["g"]], vars_Names[["y"]], vars_Names[["x"]], covariates = vars_Covariates, covariate_interactions = TRUE),
  interaction_numeric = function(df_Data, vars_Names, vars_Covariates) PlotNumInteractionEffectsMatrix(df_Data, vars_Names[["x"]], vars_Names[["y"]], interVar = vars_Names[["w"]], covariates = vars_Covariates),
  interaction_categorical = function(df_Data, vars_Names, vars_Covariates) PlotCatInteractionEffectsMatrix(df_Data, vars_Names[["g"]], vars_Names[["y"]], interVar = vars_Names[["h"]]),
  pathway_binary = function(df_Data, vars_Names, vars_Covariates) calculate_pathway_results(df_Data, vars_Names[["g"]], covariates = vars_Covariates, metabolites = vars_Names[["y"]]),
  pathway_point = function(df_Data, vars_Names, vars_Covariates) calculate_pathway_results(df_Data, vars_Names[["g"]], covariates = vars_Covariates, metabolites = vars_Names[["y"]], use_point_correlation = TRUE),
  pathway_continuous = function(df_Data, vars_Names, vars_Covariates) calculate_pathway_results(df_Data, vars_Names[["x"]], covariates = vars_Covariates, metabolites = vars_Names[["y"]]),
  biomarker = function(df_Data, vars_Names, vars_Covariates) EvaluateBiomarkerPerformance(df_Data, vars_Names[["g"]], vars_Names[["x"]], covariates = vars_Covariates, CIBootstrapR = 0, Validation = "none")
)

for (case_name in names(name_safety_calls)) {
  test_that(paste(case_name, "is invariant to unusual column names"), {
    skip_if_not_installed("emmeans")
    fixture <- NameSafetyFixture()
    Compute <- name_safety_calls[[case_name]]
    for (adjusted in c(FALSE, TRUE)) {
      reference <- Compute(fixture$Data, fixture$Names, if (adjusted) unname(fixture$Names[c("z", "w")]) else NULL)
      for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
        renamed <- NameSafetyFixture(marker)
        result <- Compute(renamed$Data, renamed$Names, if (adjusted) unname(renamed$Names[c("z", "w")]) else NULL)
        CheckNameSafetyInference(result, reference)
        CheckNameSafetyPlots(result)
        if (case_name == "beta") expect_equal(result$ResultsTable$Label, "Label x")
        if (grepl("^pathway", case_name)) {
          expect_equal(NameSafetyNumbers(result), NameSafetyNumbers(reference), tolerance = 1e-8)
          expect_identical(result$Metabolite, renamed$Names[["y"]])
        }
      }
    }
  })
}

test_that("normative predictions retain unusual covariate and practice-count names", {
  skip_if_not_installed("broom")
  Compute <- function(fixture, practice) {
    vars_Names <- fixture$Names
    CreateNormativeTScoreModel(fixture$Data, vars_Names[["y"]], vars_Names[["count"]],
      unname(vars_Names[c("z", "w")]), vars_Names[["g"]], "A", include_practice_effect = practice,
      log_transform = FALSE, return_plots = FALSE)
  }
  for (practice in c(FALSE, TRUE)) {
    reference <- Compute(NameSafetyFixture(), practice)
    for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
      result <- Compute(NameSafetyFixture(marker), practice)
      expect_equal(result$data$NormPredicted, reference$data$NormPredicted)
      expect_equal(result$data$NormT, reference$data$NormT)
      expect_equal(result$model_fit, reference$model_fit)
      expect_equal(NameSafetyNumbers(result$model_summary), NameSafetyNumbers(reference$model_summary))
    }
  }
})

test_that("module-trait models preserve unusual names for every outcome family", {
  Compute <- function(fixture) {
    df_Data <- fixture$Data
    vars_Names <- fixture$Names
    obj_Module <- structure(list(datExpr = matrix(df_Data[[vars_Names[["x"]]]], ncol = 1),
      eigengenes = data.frame(MEblue = df_Data[[vars_Names[["x"]]]])), class = "feature_wgcna_obj")
    attr(obj_Module$datExpr, "sample_ids") <- df_Data[[vars_Names[["id"]]]]
    AnalyzeFeatureModuleTraits(obj_Module, df_Data, vars_Names[["id"]],
      outcome = unname(vars_Names[c("y", "g", "h")]), covariates = unname(vars_Names[c("z", "w")]))
  }
  reference <- Compute(NameSafetyFixture())
  for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
    result <- Compute(NameSafetyFixture(marker))
    CheckNameSafetyInference(result, reference)
    expect_true(all(result$join_audit$Matched))
  }
})

test_that("all-missing association plots accept literal unusual names", {
  for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
    fixture <- NameSafetyFixture(marker)
    fixture$Data[[fixture$Names[["x"]]]] <- NA_real_
    result <- suppressWarnings(PlotAssociations(fixture$Data, fixture$Names[["x"]], fixture$Names[["y"]]))
    expect_s3_class(ggplot2::ggplot_build(result), "ggplot_built")
  }
})

test_that("statistics tables accept unusual target and measure names", {
  skip_if_not_installed("arsenal")
  skip_if_not_installed("data.table")
  skip_if_not_installed("kableExtra")
  Compute <- function(fixture) {
    vars_Names <- fixture$Names
    CreateStatisticsTable(fixture$Data[unname(vars_Names[c("g", "y")])], vars_Names[["g"]])
  }
  reference <- Compute(NameSafetyFixture())
  for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
    CheckNameSafetyInference(Compute(NameSafetyFixture(marker)), reference)
  }
})

test_that("latent class estimates are invariant to unusual indicator names", {
  skip_if_not_installed("poLCA")
  set.seed(810)
  df_Test <- as.data.frame(replicate(4, factor(sample(c("No", "Yes"), 100, TRUE))))
  Compute <- function(df_Data) CreateClusterModel_LatentClass(df_Data, method = "finalize", final_k = 2, nrep = 1, seed = 437)
  reference <- Compute(df_Test)
  for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
    df_Renamed <- df_Test
    names(df_Renamed) <- paste0(marker, names(df_Test))
    result <- Compute(df_Renamed)
    expect_equal(result$ModelInfo$lca_model$llik, reference$ModelInfo$lca_model$llik)
    expect_equal(result$ModelInfo$lca_model$posterior, reference$ModelInfo$lca_model$posterior)
    expect_equal(result$ProbFit$individual$Cluster, reference$ProbFit$individual$Cluster)
    expect_identical(names(result$ModelInfo$levels), names(df_Renamed))
  }
})

test_that("comparison term lookups handle unusual groups in robust and multi-group models", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("car")
  Compute <- function(fixture, group, parametric, adjusted, alternative) {
    vars_Names <- fixture$Names
    MakeComparisonTable(fixture$Data, vars_Names[[group]], vars_Names[["y"]],
      covariates = if (adjusted) unname(vars_Names[c("z", "w")]) else NULL,
      Parametric = parametric, alternative = alternative)
  }
  for (group in c("g", "h")) {
    for (parametric in c(FALSE, TRUE)) {
      for (adjusted in c(FALSE, TRUE)) {
        for (alternative in if (group == "g") c("two.sided", "greater") else "two.sided") {
          reference <- Compute(NameSafetyFixture(), group, parametric, adjusted, alternative)
          for (marker in c("% ", "with space ", "a+b", "tick`", "slash\\")) {
            result <- Compute(NameSafetyFixture(marker), group, parametric, adjusted, alternative)
            CheckNameSafetyInference(result, reference)
          }
        }
      }
    }
  }
})

test_that("raw names and codebook labels survive model and plotting output", {
  fixture <- NameSafetyFixture("tick`")
  vars_Names <- fixture$Names
  for (Compute in list(PlotBetaProfile, PlotVolcanoEffects)) {
    result <- Compute(fixture$Data, vars_Names[["x"]], vars_Names[["y"]],
      covariates = unname(vars_Names[c("z", "w")]), InteractiveLabels = FALSE)
    expect_identical(result$Results$Predictor, vars_Names[["x"]])
    expect_identical(result$Results$Outcome, vars_Names[["y"]])
    expect_identical(result$Results$PredictorLabel, "Label x")
    expect_identical(result$Results$OutcomeLabel, "Label y")
    expect_identical(result$Metadata$Arguments$covariates, unname(vars_Names[c("z", "w")]))
  }
  # Without display labels, the legacy interaction matrices use raw names in reshaping.
  for (variable in names(fixture$Data)) attr(fixture$Data[[variable]], "label") <- NULL
  result <- PlotNumInteractionEffectsMatrix(fixture$Data, vars_Names[["x"]], vars_Names[["y"]],
    interVar = vars_Names[["w"]], covariates = vars_Names[["z"]])
  expect_true(all(result$Results$Status == "ok"))
  expect_identical(result$Results$Outcome, vars_Names[["y"]])
  CheckNameSafetyPlots(result)
})
