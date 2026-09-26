test_that("CreateClusterModel_SOM_MClust returns the complete model-review output set", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory",
    k_range = 2:4,
    models = 1,
    lpa_timeout_seconds = 30
  )

  expect_s3_class(model$fit_plot, "ggplot")
  expect_s3_class(model$ModelInfo_SOM$plots$Circular, "htmlwidget")
  expect_s3_class(model$ModelInfo_SOM$plots$Line, "htmlwidget")
  expect_s3_class(model$ModelInfo_SOM$plots$Cloud, "htmlwidget")

  som_fit_plots <- model$ModelInfo_SOM$SOMFit$plots
  reference_plots <- model$ModelInfo_SOM$ProjectionReference$plots
  probability_plots <- model$ProbFit$plots

  expect_true(all(vapply(som_fit_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(reference_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(probability_plots, inherits, logical(1), what = "ggplot")))
  expect_false(any(grepl("MaxProbBoxplot", names(probability_plots))))
  expect_true(all(c(
    "node_ProbAssignedDensity", "individual_ProbAssignedDensity"
  ) %in% names(probability_plots)))
  expect_equal(
    probability_plots$node_ProbAssignedDensity$scales$get_scales("x")$limits,
    c(0, 1)
  )
  expect_equal(
    probability_plots$individual_ProbAssignedDensity$scales$get_scales("x")$limits,
    c(0, 1)
  )

  fit_table <- model$ModelInfo_MClust$fit_table
  expect_false("CAICCLC" %in% names(fit_table))
  expect_true(all(c("CAIC", "CLC") %in% names(fit_table)))
  expect_true(all(c("MinProfileNodeN", "MaxProfileNodeN",
    "MinProfileNodeProportion", "MaxProfileNodeProportion",
    "BLRTStatistic", "BLRTPValue") %in% names(fit_table)))
  expect_true(all(fit_table$MinProfileNodeN == as.integer(fit_table$MinProfileNodeN)))
  expect_true(all(fit_table$MaxProfileNodeN == as.integer(fit_table$MaxProfileNodeN)))
  expect_false(any(c("n_min", "n_max", "BLRT_val", "BLRT_p") %in% names(fit_table)))
  expect_true("Bootstrap likelihood-ratio test p-value (−log10 scale)" %in%
    as.character(model$fit_plot$data$name))
  expect_equal(
    model$fit_plot$data$PlotValue[model$fit_plot$data$name ==
      "Bootstrap likelihood-ratio test p-value (−log10 scale)"],
    -log10(pmax(fit_table$BLRTPValue, .Machine$double.xmin)))
  expected_best <- fit_table[which.max(fit_table$ahp_index), , drop = FALSE]
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Model, expected_best$Model)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, expected_best$Classes)
})

test_that("SOM Mclust lifecycle settings are mutually exclusive", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  expect_error(CreateClusterModel_SOM_MClust(
    df_Labelled, c("age", "AXL", "Adiponectin"), method = "finalize",
    k_range = 2:10, models = c(1, 3), final_k = 2, final_model = 1,
    lpa_timeout_seconds = NULL), "k_range, models")
})

test_that("CreateClusterModel_SOM_MClust handles a single exploratory candidate", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory",
    k_range = 2,
    models = 1,
    lpa_timeout_seconds = 30
  )

  expect_equal(nrow(model$ModelInfo_MClust$fit_table), 1)
  expect_equal(model$ModelInfo_MClust$fit_table$ahp_index, 0)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, 2)
})

test_that("SOM Mclust supports tidyLPA model 6", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled,
    c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "finalize", final_k = 2, final_model = 6,
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30, Relabel = FALSE
  ))

  expect_identical(model$ModelInfo_MClust$fit_table$Model, 6L)
  expect_true(any(grepl(
    "Model 6: varying variance, varying covariance",
    levels(model$fit_plot$data$Model), fixed = TRUE
  )))

  mixed_model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled,
    c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory", k_range = 2, models = c(1, 6),
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30, Relabel = FALSE
  ))
  model_six_label <- "Model 6: varying variance, varying covariance"
  colour_scale <- ggplot2::ggplot_build(mixed_model$fit_plot)$plot$scales$get_scales("colour")

  expect_true(all(c(
    "Model 1: equal variance, zero covariance", model_six_label
  ) %in% colour_scale$get_breaks()))
  expect_identical(
    unname(colour_scale$map(model_six_label)),
    unname(.SciDataColorValues(4)[[4]])
  )
})

test_that("CreateClusterModel_SOM_MClust adds deterministic subsample stability to fit review", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")

  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  fit_model <- function() {
    CreateClusterModel_SOM_MClust(
      data = df_Labelled,
      variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
      method = "exploratory",
      k_range = 2,
      models = 1,
      stability_resamples = 2,
      stability_seed = 1234,
      lpa_timeout_seconds = 30
    )
  }

  model_one <- fit_model()
  model_two <- fit_model()
  stability <- model_one$ModelInfo_MClust$Stability
  fit_table <- model_one$ModelInfo_MClust$fit_table

  expect_equal(stability$settings$resamples, 2)
  expect_equal(nrow(stability$replicates), 2)
  expect_true(all(stability$replicates$Status == "success"))
  expect_true(all(stability$cluster_recovery$Jaccard >= 0))
  expect_true(all(stability$cluster_recovery$Jaccard <= 1))
  expect_true(all(c(
    "StabilitySuccessRate", "StabilityARI_Mean", "StabilityARI_P05",
    "StabilityJaccard_Mean", "StabilityJaccard_Min",
    "ReproducibilityScore", "Reproducibility_scaled"
  ) %in% names(fit_table)))
  expect_false("CAICCLC" %in% names(fit_table))
  expect_true(all(c("CAIC", "CLC") %in% names(fit_table)))
  expect_false(any(grepl("MaxProbBoxplot", names(model_one$ProbFit$plots))))
  expect_true("ReproducibilityScore" %in% unique(model_one$fit_plot$data$name))
  expect_equal(
    model_one$ModelInfo_MClust$Stability$replicates$ARI,
    model_two$ModelInfo_MClust$Stability$replicates$ARI
  )
  expect_equal(mclust::adjustedRandIndex(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  negative_ari <- mclust::adjustedRandIndex(
    c(1, 1, 1, 2), c(1, 2, 2, 2)
  )
  expect_lt(negative_ari, 0)
  expect_equal(mean(c(negative_ari, 0.8)), (negative_ari + 0.8) / 2)
})

test_that("SOM Mclust finalize with stability records the user-specified model", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- CreateClusterModel_SOM_MClust(
    data = df_Labelled,
    variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "finalize", final_k = 2, final_model = 1,
    stability_resamples = 2, stability_cores = 1,
    lpa_timeout_seconds = 30, Relabel = FALSE
  )

  expect_identical(model$method, "finalize")
  expect_equal(model$Specification$selected, list(k = 2, model = 1))
  expect_match(model$ModelInfo_MClust$AHP$recommendation, "^User-specified Model 1")
  expect_true(all(is.na(model$ModelInfo_MClust$fit_table$ahp_index)))
  expect_true("StabilityARI_Mean" %in% names(model$ModelInfo_MClust$fit_table))
  expect_identical(model$ModelInfo, model$ModelInfo_MClust)

  # Individuals carry node assignments, not the node codebook vectors.
  individual <- model$ProbFit$individual
  expect_false(any(grepl("^Z_", names(individual))))
  expect_true(all(c("Cluster", "prob_1", "prob_assigned") %in% names(individual)))
  expect_true(any(grepl("^Z_", names(model$ProbFit$node))))
  expect_true(all(individual$Projection_Fit_Class %in% c(NA, "Good fit",
    "Uncertain membership", "Poor fit to training structure",
    "Potential novel phenotype")))

  projected <- ProjectCluster(model, df_Labelled)
  expect_false(any(grepl("^Z_", names(projected$ProbFit$individual))))
})

test_that("SOM Mclust records inestimable fits as failed instead of stopping", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  # Model 6 with 5 profiles needs more parameters than a 4 x 4 SOM has nodes,
  # so mclust returns fit indices without class assignments.
  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled,
    c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory", k_range = c(2, 5), models = c(1, 6),
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30, Relabel = FALSE
  ))

  diagnostics <- model$ModelInfo_MClust$diagnostics$lpa_fit_diagnostics
  failed <- diagnostics[diagnostics$Model == 6 & diagnostics$Classes == 5, ]
  expect_identical(failed$status, "failed")
  expect_match(failed$error, "no class assignments")
  expect_false(any(model$ModelInfo_MClust$fit_table$Model == 6 &
    model$ModelInfo_MClust$fit_table$Classes == 5))
  expect_true(all(c(2L, 5L) %in%
    model$ModelInfo_MClust$fit_table$Classes[model$ModelInfo_MClust$fit_table$Model == 1]))
})

test_that("SOM Mclust reports Euclidean distances with per-variable residuals", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  vars <- c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin")

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled, vars, method = "finalize", final_k = 2, final_model = 1,
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  individual <- model$ProbFit$individual
  mapped <- !is.na(individual$SOM_Distance)
  expect_equal(individual$SOM_Distance[mapped],
    sqrt(model$ModelInfo_SOM$som_model$distances))

  residuals <- model$ModelInfo_SOM$SOMFit$residuals
  expect_identical(names(residuals), c(".row_id", paste0("Resid_", vars)))
  resid_sq <- rowSums(as.matrix(residuals[paste0("Resid_", vars)])^2)
  expect_equal(resid_sq, individual$SOM_Distance[match(residuals$.row_id, individual$.row_id)]^2)
  expect_true(all(individual$Top_Distance_Variable[mapped] %in% vars))
  contribution <- model$ModelInfo_SOM$SOMFit$variable_contribution
  expect_equal(sum(contribution$share_of_distance), 1)

  # Projecting the training data onto its own map reproduces training fit.
  projected <- suppressWarnings(ProjectCluster(model, df_Labelled))
  drift <- projected$ProjectionFit$variable_drift
  expect_equal(drift$sq_residual_ratio, rep(1, length(vars)))
  summary <- projected$ProjectionFit$summary
  expect_equal(summary$value[summary$metric == "cluster_occupancy_js_divergence"], 0)
  expect_equal(summary$value[summary$metric == "excess_high_distance"],
    summary$value[summary$metric == "high_distance_burden"] - 0.05)
  expect_s3_class(projected$ProjectionFit$plots$variable_drift, "ggplot")
})

test_that("SOM Mclust can leave sparse nodes out of the mixture fit", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  skip_if_not_installed("R.utils")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Labelled, c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
    method = "exploratory", k_range = 2:3, models = 1,
    som_xdim = 8, som_ydim = 8, min_nodes_per_cluster = NULL,
    lpa_min_node_n = 3, min_cluster_prop = 0.2,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  node <- model$ProbFit$node
  expect_identical(node$InLPAFit, node$NodeN >= 3)
  expect_true(any(!node$InLPAFit))
  expect_false(anyNA(node$Cluster))
  expect_equal(unname(rowSums(as.matrix(node[grep("^prob_[0-9]+$", names(node))]))),
    rep(1, nrow(node)), tolerance = 1e-6)
  expect_equal(model$ModelInfo_MClust$diagnostics$lpa_preprocess$n_lpa_nodes_used,
    sum(node$NodeN >= 3))

  fit_table <- model$ModelInfo_MClust$fit_table
  expect_true(all(c("MinProfileParticipantN", "MinProfileParticipantProportion",
    "Eligible") %in% names(fit_table)))
  expect_identical(fit_table$Eligible, fit_table$MinProfileParticipantProportion >= 0.2)
  expect_match(model$ModelInfo_MClust$AHP$recommendation, "^Composite rank index")
})

test_that("SOM Mclust composite index accepts AHP pairwise weights", {
  pairwise <- matrix(c(1, 1/3, 1/2,
                       3, 1, 2,
                       2, 1/2, 1), 3, byrow = TRUE,
    dimnames = list(c("AIC", "BIC", "Entropy"), c("AIC", "BIC", "Entropy")))
  weights <- .ResolveCompositeWeights(pairwise)
  expect_identical(weights$method, "AHP pairwise")
  expect_equal(sum(weights$weights), 1)
  expect_identical(names(which.max(weights$weights)), "BIC")
  expect_lt(weights$consistency_ratio, 0.10)

  expect_warning(.ResolveCompositeWeights(matrix(c(1, 9, 1/9,
                                                   1/9, 1, 9,
                                                   9, 1/9, 1), 3, byrow = TRUE,
    dimnames = list(c("AIC", "BIC", "Entropy"), c("AIC", "BIC", "Entropy")))),
    "inconsistent")
  expect_error(.ResolveCompositeWeights(matrix(c(1, 2, 2, 1), 2,
    dimnames = list(c("AIC", "BIC"), c("AIC", "BIC")))), "reciprocal")
  expect_equal(.ResolveCompositeWeights(c(BIC = 2, Entropy = 2))$weights,
    c(BIC = 0.5, Entropy = 0.5))
})
