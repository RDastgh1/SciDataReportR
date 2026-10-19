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
    lpa_bootstrap_n = 20,
    stability_resamples = 0,
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
  expect_identical(.ClusterMetricDirection(c("BIC", "ICL"), fit_table),
    c("lower", "higher"))
  expect_true("BIC (lower is better)" %in%
    model$fit_plot$facet$params$labeller(data.frame(name = "BIC"))[[1]])
  expected_best <- fit_table[which.max(fit_table$ahp_index), , drop = FALSE]
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Model, expected_best$Model)
  expect_equal(model$ModelInfo_MClust$AHP$ahp_best_row$Classes, expected_best$Classes)
})

test_that("SOM stability fractions retain full-cohort scoring and expose partition diagnostics", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    assignment = "soft", stability_resamples = 2L,
    stability_fractions = c(.90, .80), stability_seed = 889L,
    stability_cores = 1L, Relabel = FALSE, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30
  ))
  stability <- model$Stability
  expect_setequal(unique(stability$replicates$RetainedFraction), c(.90, .80))
  expect_true(all(stability$replicates$ParticipantsCompared == nrow(df_Training)))
  expect_setequal(unique(stability$deletion_curve$RetainedFraction), c(.90, .80))
  expect_true(all(c("VI", "NMI", "FowlkesMallows") %in% names(stability$replicates)))
  expect_true(nrow(stability$co_membership_by_cluster) > 0)
  expect_true(nrow(stability$split_merge_reference) > 0)
  expect_true(nrow(stability$profile_recovery) > 0)
  expect_error(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, stability_fractions = c(.90, 1)),
    "strictly between 0 and 1")
})

test_that("co-membership and profile recovery are label invariant", {
  reference <- c(1L, 1L, 1L, 2L, 2L, 2L)
  refit <- c(7L, 7L, 7L, 9L, 9L, 9L)
  profiles <- cbind(
    x = c(0, 0, 0, 1, 1, 1),
    y = c(1, 1, 1, 0, 0, 0),
    z = c(2, 2, 2, 3, 3, 3)
  )
  co_membership <- SciDataReportR:::.ClusterCoMembershipSummary(reference, refit)
  split_merge <- SciDataReportR:::.ClusterSplitMergeDiagnostics(reference, refit)
  recovery <- SciDataReportR:::.ClusterProfileRecovery(reference, refit, profiles)
  expect_true(all(co_membership$by_cluster$Cohesion == 1))
  expect_true(all(co_membership$by_cluster$CrossMixing == 0))
  expect_true(all(split_merge$reference$DominantRefitShare == 1))
  expect_true(all(split_merge$reference$EffectiveRefitClusters == 1))
  expect_true(all(split_merge$reference$BestJaccard == 1))
  expect_equal(recovery$ProfileCorrelation, rep(1, 2), tolerance = 1e-12)
  expect_equal(recovery$ProfileRMSE, rep(0, 2), tolerance = 1e-12)
})

test_that("co-membership summaries equal pairwise hand calculations without a participant matrix", {
  reference <- c(1L, 1L, 2L, 2L)
  refit <- c(7L, 9L, 9L, 9L)
  summary <- SciDataReportR:::.ClusterCoMembershipSummary(reference, refit)
  expect_equal(summary$by_cluster$Cohesion, c(0, 1))
  expect_equal(summary$by_cluster$CrossMixing, c(.5, .5))
  expect_equal(summary$between_clusters$CoassignmentProbability, .5)

  large_summary <- SciDataReportR:::.ClusterCoMembershipSummary(
    rep(seq_len(5L), each = 1000L), rep(seq_len(5L), each = 1000L))
  expect_equal(nrow(large_summary$by_cluster), 5L)
  expect_false(any(vapply(large_summary, is.matrix, logical(1))))
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
    stability_resamples = 0,
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
    stability_resamples = 0,
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
    stability_resamples = 0,
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

  expect_identical(formals(CreateClusterModel_SOM_MClust)$stability_resamples, 0L)
  expect_equal(stability$settings$resamples, 2)
  expect_identical(stability$settings$comparison_scope, "all_reference_participants")
  expect_identical(stability$settings$som_grid$initialization_seed_policy,
    "fixed_to_reference_fit")
  expect_equal(nrow(stability$replicates), 2)
  expect_true(all(stability$replicates$Status == "success"))
  expect_true(all(stability$replicates$ParticipantsCompared ==
    sum(!is.na(model_one$DataWithClusters$Cluster))))
  expect_true(all(stability$cluster_recovery$Jaccard >= 0))
  expect_true(all(stability$cluster_recovery$Jaccard <= 1))
  expect_equal(fit_table$ReproducibilityScore, fit_table$StabilityARI_Mean)
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
    stability_resamples = 0,
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
    stability_resamples = 0,
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
    stability_resamples = 0,
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

test_that("SOM Mclust stability can be limited to a shortlist of candidates", {
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
    method = "exploratory", k_range = 2:4, models = 1,
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    stability_resamples = 2, stability_top_n = 2, stability_cores = 1,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  fit_table <- model$ModelInfo_MClust$fit_table
  assessed <- !is.na(fit_table$ReproducibilityScore)
  expect_equal(sum(assessed), 2)
  expect_length(model$Stability$settings$candidates_assessed, 2)
  expect_equal(nrow(model$Stability$replicates), 2 * 2)
  expect_true(all(is.na(fit_table$ahp_index[!assessed])))
  best <- model$ModelInfo_MClust$AHP$ahp_best_row
  expect_true(paste0("model_", best$Model, "_class_", best$Classes) %in%
    model$Stability$settings$candidates_assessed)
  expect_length(model$ModelInfo_MClust$node_assignments, 3)
})

test_that("hard SOM Mclust assignment equals the best-matching node class", {
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
    df_Labelled, vars, method = "finalize", final_k = 3, final_model = 1,
    som_xdim = 5, som_ydim = 5, stability_resamples = 0,
    assignment = "hard",
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30, Relabel = FALSE))
  individual <- model$ProbFit$individual
  mapped <- !is.na(individual$SOM_Node)

  expect_identical(individual$Cluster[mapped],
    model$ProbFit$node$Cluster[individual$SOM_Node[mapped]])
  expect_identical(individual$NodeCluster[mapped],
    model$ProbFit$node$Cluster[individual$SOM_Node[mapped]])
  expect_false(any(individual$Flag_NodeDisagreement))
  expect_identical(model$DataWithClusters$Cluster, individual$Cluster)

  # Projecting the training data preserves the node-inherited assignment.
  projected <- suppressWarnings(ProjectCluster(model, df_Labelled))
  expect_identical(projected$DataWithClusters$Cluster, model$DataWithClusters$Cluster)
  expect_true("Flag_NodeDisagreement" %in% names(projected$ProbFit$individual))
})

test_that("SOM Mclust stability can assess the best candidate at each k", {
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
    method = "exploratory", k_range = 2:3, models = c(1, 2),
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    stability_resamples = 2, stability_candidates = "best_per_k", stability_cores = 1,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  assessed <- model$ModelInfo_MClust$fit_table %>%
    dplyr::filter(!is.na(ReproducibilityScore))
  expect_setequal(assessed$Classes, 2:3)
  expect_equal(nrow(assessed), 2)
})

test_that("SOM Mclust summarizes posteriors and core stability per candidate", {
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
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    stability_resamples = 2, stability_cores = 1, low_prob_threshold = 0.7,
    criteria_weights = c(Entropy = 1, Posterior = 1, Reproducibility = 1),
    lpa_timeout_seconds = 30, Relabel = FALSE))

  fit_table <- model$ModelInfo_MClust$fit_table
  expect_true(all(c("PosteriorMean", "PosteriorMedian", "PropLowPosterior",
    "StabilityCoreProportion", "StabilityARI_Core", "StabilityJaccard_CoreMin",
    "ReproducibilityScore_All") %in% names(fit_table)))
  expect_true(all(fit_table$PosteriorMean > 0 & fit_table$PosteriorMean <= 1))

  # Posterior summaries inherit the assigned node's posterior.
  best <- model$ModelInfo_MClust$AHP$ahp_best_row
  individual <- model$ProbFit$individual
  row <- fit_table[fit_table$Model == best$Model & fit_table$Classes == best$Classes, ]
  expect_equal(row$PosteriorMean, mean(individual$max_prob, na.rm = TRUE))
  expect_equal(row$PropLowPosterior, mean(individual$max_prob < 0.7, na.rm = TRUE))
  expect_true(all(fit_table$StabilityCoreProportion >= 0 & fit_table$StabilityCoreProportion <= 1))
  expect_match(model$ModelInfo_MClust$AHP$recommendation, "Posterior")
})

test_that("SOM Mclust skips tidyLPA bootstraps by default", {
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
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    stability_resamples = 0,
    lpa_timeout_seconds = 30, Relabel = FALSE))

  fit_table <- model$ModelInfo_MClust$fit_table
  expect_true(all(is.na(fit_table$BLRTPValue)))
  expect_false(any(grepl("likelihood-ratio", as.character(model$fit_plot$data$name))))
  expect_error(CreateClusterModel_SOM_MClust(df_Labelled, "age", lpa_bootstrap_n = -1),
    "lpa_bootstrap_n")
})

test_that("SOM hard and soft projection retain participant confidence outputs", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    soft_n_nodes = 3L,
    stability_resamples = 0L, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))

  individual <- model$ProbFit$individual
  mapped <- !is.na(individual$SOM_Node)
  soft_probability_columns <- grep("^SoftProb_", names(individual), value = TRUE)
  expected_columns <- c("Cluster_Hard", "Cluster_Soft", "SoftEntropy",
    "SoftMargin", "Flag_BoundaryParticipant", "BoundaryReason")
  expect_true(all(expected_columns %in% names(individual)))
  expect_true(all(expected_columns %in% names(model$DataWithClusters)))
  expect_equal(rowSums(as.matrix(individual[mapped, soft_probability_columns])),
    rep(1, sum(mapped)), tolerance = 1e-8)
  expect_identical(individual$Cluster, individual$Cluster_Soft)
  expect_identical(individual$Cluster_Hard[mapped],
    model$ProbFit$node$Cluster[individual$SOM_Node[mapped]])
  expect_true(all(individual$SoftEntropy[mapped] >= 0 & individual$SoftEntropy[mapped] <= 1))
  expect_true(all(individual$SoftMargin[mapped] >= 0 & individual$SoftMargin[mapped] <= 1))
  expect_true(all(individual$BoundaryReason[!is.na(individual$BoundaryReason)] %in%
    c("low_margin", "high_entropy", "both")))

  projected <- suppressWarnings(ProjectCluster(model, df_Training[1:20, ]))
  expect_true(all(expected_columns %in% names(projected$ProbFit$individual)))
  expect_true(all(expected_columns %in% names(projected$DataWithClusters)))
  expect_identical(projected$DataWithClusters$Cluster,
    projected$ProbFit$individual$Cluster)

  # Earlier objects may retain now-unsupported training provenance. Projection
  # must use their frozen codebook and node mixture rather than inspect it.
  legacy_model <- model
  legacy_model$ModelInfo_SOM$consensus <- list(strategy = "median_only")
  legacy_model$Specification$projection$consensus_strategy <- "median_only"
  legacy_projection <- suppressWarnings(ProjectCluster(legacy_model, df_Training[1:20, ]))
  expect_identical(legacy_projection$ProbFit$individual$Cluster,
    projected$ProbFit$individual$Cluster)

  node_df <- tibble::tibble(NodeID = 1:3, Cluster = c(1L, 2L, 3L),
    prob_1 = c(0.8, 0.1, 0.2), prob_2 = c(0.1, 0.8, 0.2),
    prob_3 = c(0.1, 0.1, 0.6))
  soft_at_node <- .SOMSoftProjection(matrix(c(0, 0), nrow = 1),
    rbind(c(0, 0), c(1, 0), c(0, 1)), node_df, n_nodes = 3L)
  expect_equal(unname(soft_at_node$probability[1, ]), c(0.8, 0.1, 0.1))
  expect_equal(soft_at_node$cluster, 1L)
})

test_that("standard batch-SOM provenance and boundary stability use a fixed reference mask", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData", package = "SciDataReportR")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")

  standard_one <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    seed_som = 889L, stability_resamples = 0L,
    Relabel = FALSE, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  standard_two <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    seed_som = 889L, stability_resamples = 0L,
    Relabel = FALSE, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  expect_false("consensus_resamples" %in% names(formals(CreateClusterModel_SOM_MClust)))
  expect_false("consensus_strategy" %in% names(formals(CreateClusterModel_SOM_MClust)))
  expect_identical(eval(formals(CreateClusterModel_SOM_MClust)$assignment),
    c("soft", "hard"))
  expect_identical(formals(CreateClusterModel_SOM_MClust)$soft_n_nodes, 3L)
  expect_identical(formals(CreateClusterModel_SOM_MClust)$stability_resamples, 0L)
  expect_equal(standard_one$ModelInfo_SOM$som_codes,
    standard_two$ModelInfo_SOM$som_codes)
  expect_identical(standard_one$ModelInfo_SOM$training$mode, "batch")
  expect_identical(standard_one$ModelInfo_SOM$training$epochs, 100L)
  expect_equal(dim(standard_one$ModelInfo_SOM$training$initial_codebook),
    dim(standard_one$ModelInfo_SOM$som_codes))
  expect_error(do.call(CreateClusterModel_SOM_MClust, list(
    data = df_Training, variables = paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, consensus_resamples = 1L)),
    "unused argument")

  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    stability_resamples = 2L,
    stability_seed = 889L, stability_cores = 1L, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  stability <- model$Stability
  expect_true(all(stability$replicates$ParticipantsCompared == nrow(df_Training)))
  expect_true(all(c("ARI_NonBoundary", "NonBoundaryProportion") %in%
    names(stability$replicates)))
  expect_true(all(c("StabilityARI_NonBoundary", "StabilityARI_NonBoundary_P05",
    "StabilityJaccard_NonBoundaryMean") %in% names(stability$summary)))
  expect_true(all(c("overall", "reference_cluster") %in%
    stability$boundary_summary$Scope))
  expect_identical(stability$settings$non_boundary_mask,
    "reference_soft_margin_and_entropy")
  expect_identical(stability$settings$som_training$mode, "batch")
  expect_identical(stability$settings$som_training$epochs, 100L)
  expect_equal(stability$settings$assignment, "soft")
  expect_identical(model$ModelInfo_SOM$som_grid_info$GridDesignPolicy,
    "explicit_dimensions")
  expect_equal(model$ModelInfo_SOM$som_grid_info$ExpectedNodesPerProfile,
    16 / 3)

  parallel_model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1, som_xdim = 4, som_ydim = 4,
    assignment = "soft", stability_resamples = 2L,
    stability_seed = 889L, stability_cores = 2L, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  expect_identical(stability$replicates$ARI, parallel_model$Stability$replicates$ARI)
  expect_identical(stability$replicates$ARI_NonBoundary,
    parallel_model$Stability$replicates$ARI_NonBoundary)

  adaptive <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "finalize",
    final_k = 3, final_model = 1,
    stability_resamples = 0L, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  expected_side <- SOMMapSupport(nrow(df_Training), 6, 3, 1)$RecommendedSide[[1]]
  expect_equal(adaptive$ModelInfo_SOM$som_grid_info$som_xdim_used, expected_side)
  expect_equal(adaptive$ModelInfo_SOM$som_grid_info$som_ydim_used, expected_side)
  expect_identical(adaptive$ModelInfo_SOM$som_grid_info$GridDesignK, 3L)
  expect_identical(adaptive$ModelInfo_SOM$som_grid_info$TargetNodesPerProfile, 20L)
  expect_identical(adaptive$ModelInfo_SOM$som_grid_info$GridDesignPolicy,
    "k_informed_parameter_and_occupancy_aware")

  exploratory <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:6), method = "exploratory",
    k_range = 3:5, models = 1, stability_resamples = 0L,
    Relabel = FALSE, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30
  ))
  expect_identical(exploratory$ModelInfo_SOM$som_grid_info$GridDesignK, 5L)
  expected_exploratory_side <- SOMMapSupport(nrow(df_Training), 6, 3:5, 1)$RecommendedSide[[1]]
  expect_equal(exploratory$ModelInfo_SOM$som_grid_info$som_xdim_used,
    expected_exploratory_side)
  expect_equal(exploratory$ModelInfo_SOM$som_grid_info$som_ydim_used,
    expected_exploratory_side)
  expect_true(all(exploratory$ModelInfo_MClust$fit_table$Classes %in% 3:5))
})

test_that("SOM map support combines parameter and occupancy constraints", {
  expect_equal(SciDataReportR:::.SOMMclustParameterCount(2, 15, 5),
    mclust::nMclustParams("VVI", 15, 5))
  support <- SOMMapSupport(18196, 15, 5, 2)
  expect_equal(support$RecommendedSide[[1]], ceiling(sqrt(3 * 154)))
  expect_identical(support$BindingConstraint[[1]], "parameter_floor")
  capped <- SOMMapSupport(400, 15, 5, 2)
  expect_equal(capped$RecommendedSide[[1]], floor(sqrt(400 / 10)))
  expect_identical(capped$BindingConstraint[[1]], "participant_occupancy")
  expect_identical(capped$SupportStatus[[1]], "over-parameterized")

  set.seed(6031)
  df_ParameterTest <- as.data.frame(matrix(stats::rnorm(200 * 6), ncol = 6))
  names(df_ParameterTest) <- paste0("V", seq_len(6))

  df_ExplicitTest <- df_ParameterTest[, seq_len(6)]
  expect_warning(explicit <- CreateClusterModel_SOM_MClust(
    df_ExplicitTest, variables = names(df_ExplicitTest), method = "exploratory",
    k_range = 5, models = 1, som_xdim = 6, som_ydim = 6,
    stability_resamples = 0L, Relabel = FALSE, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = 30
  ), "Explicit SOM dimensions")
  expect_equal(explicit$ModelInfo_SOM$som_grid_info$som_xdim_used, 6L)
  expect_true(explicit$ModelInfo_MClust$fit_table$Overparameterized[[1]])
  expect_false(explicit$ModelInfo_MClust$fit_table$Eligible[[1]])
  expect_true(all(c("ParticipantsPerNode", "OccupancyCompliant",
    "SupportStatus") %in% names(explicit$ModelInfo_MClust$fit_table)))
})

test_that("map_from reuses map design but not training preprocessing", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  set.seed(6032)
  df_Explore <- data.frame(A = stats::rnorm(150), B = stats::rnorm(150))
  exploratory <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Explore, c("A", "B"), k_range = 2, models = 1,
    som_xdim = 4, som_ydim = 4, Relabel = FALSE,
    lpa_timeout_seconds = NULL))
  df_Final <- dplyr::mutate(df_Explore, A = .data$A + 100)
  finalized <- suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Final, c("A", "B"), method = "finalize", final_k = 2,
    final_model = 1, map_from = exploratory, Relabel = FALSE,
    lpa_timeout_seconds = NULL))
  expect_equal(finalized$ModelInfo_SOM$som_grid_info$som_xdim_used, 4L)
  expect_equal(finalized$ModelInfo_SOM$som_grid_info$som_ydim_used, 4L)
  expect_identical(finalized$ModelInfo_SOM$som_grid$topo,
    exploratory$ModelInfo_SOM$som_grid$topo)
  expect_identical(as.character(finalized$ModelInfo_SOM$som_grid$neighbourhood.fct),
    as.character(exploratory$ModelInfo_SOM$som_grid$neighbourhood.fct))
  expect_identical(finalized$ModelInfo_SOM$training$seed,
    exploratory$ModelInfo_SOM$training$seed)
  expect_identical(finalized$ModelInfo_SOM$som_grid_info$BindingConstraint,
    "map_from")
  expect_false(identical(finalized$ZScoreObject, exploratory$ZScoreObject))
})

test_that("fit review plots expose SOM mixture support status", {
  fit_table <- data.frame(Classes = c(2L, 3L, 4L), Model = 1L,
    BIC = c(-10, -12, -14),
    SupportStatus = c("supported", "under-supported", "over-parameterized"))
  plot <- PlotClusterFitReview(fit_table, metrics = "BIC")
  expect_identical(ggplot2::ggplot_build(plot)$plot$scales$
    get_scales("shape")$name, "SOM mixture support")
})

test_that("tidyLPA ICL is higher-is-better while its BIC is lower-is-better", {
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  X <- scale(datasets::iris[, 1:4])
  fit <- .FitTidyLPAMclust(X, 3L, 1L, nboot = 0L)
  fit_info <- tidyLPA::get_fit(fit)
  z_map <- apply(fit[[1]]$model$z, 1, max)
  # tidyLPA: BIC = -2 logLik + penalty, ICL = -BIC + 2 sum(log max z), i.e.
  # mclust's ICL on its higher-is-better scale.
  expect_equal(fit_info$BIC,
    -2 * fit[[1]]$model$loglik + fit[[1]]$model$df * log(nrow(X)))
  expect_equal(fit_info$ICL, -fit_info$BIC + 2 * sum(log(z_map + 1e-12)),
    tolerance = 1e-6)
  expect_equal(fit_info$ICL, mclust::icl(fit[[1]]$model), tolerance = 1e-6)
  expect_lt(fit_info$ICL, 0)
  expect_gt(fit_info$BIC, 0)

  som_table <- .SetClusterMetricDirection(
    data.frame(Classes = 2:3, Model = 1L, BIC = c(880, 876), ICL = c(-890, -884)),
    higher = "ICL", lower = "BIC")
  plot <- PlotClusterFitReview(som_table, metrics = c("BIC", "ICL"))
  labels <- plot$facet$params$labeller(data.frame(Metric = c("BIC", "ICL")))[[1]]
  expect_identical(labels, c("BIC (lower is better)", "ICL (higher is better)"))

  # mclust's own BIC is higher-is-better, and an override wins over the table.
  mclust_table <- .ClusterAHP(som_table, higher = c("BIC", "ICL"))$fit_table
  expect_identical(.ClusterMetricDirection(c("BIC", "ICL"), mclust_table),
    c("higher", "higher"))
  expect_identical(
    .ClusterMetricDirection("BIC", mclust_table, c(BIC = "lower")), "lower")
  # Without a recorded direction BIC is left unlabelled rather than guessed.
  expect_identical(.ClusterMetricDirectionLabel(c("BIC", "ICL")),
    c("BIC", "ICL (higher is better)"))
})
