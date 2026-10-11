test_that("SOM grid candidates include k-informed resolutions", {
  expect_identical(eval(formals(ScreenSOMGridStability)$screen_scope),
    c("map", "joint"))
  expect_identical(formals(ScreenSOMGridStability)$screen_resamples, 3L)
  grids <- SciDataReportR:::.SOMGridCandidates(1000, k_range = 3:6,
    policy = "k_informed")
  k_informed_sides <- ceiling(sqrt(max(3:6) * c(10, 20, 30)))

  expect_true(all(c(3L, 10L, k_informed_sides) %in%
    grids$som_xdim))
  expect_true(all(grids$som_xdim == grids$som_ydim))
  expect_true(all(c(8L, 11L, 14L) %in%
    SciDataReportR:::.SOMGridCandidates(18196, k_range = 3:6,
      policy = "k_informed")$som_xdim))
  explicit <- SciDataReportR:::.SOMGridCandidates(1000,
    grid_candidates = data.frame(som_xdim = c(4L, 4L), som_ydim = c(6L, 6L)),
    k_range = 3:6, policy = "k_informed")
  expect_identical(as.data.frame(explicit),
    data.frame(som_xdim = 4L, som_ydim = 6L))
})

test_that("default SOM map screen is k-agnostic and projects the full cohort", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  vars_SOM <- c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin")

  screen <- suppressWarnings(ScreenSOMGridStability(
    data = df_Labelled, variables = vars_SOM, k_range = 2, models = 1,
    grid_candidates = 3, screen_resamples = 1, stability_cores = 1,
    pair_sample = 50, neighborhood_sample = 20, lpa_timeout_seconds = NULL
  ))

  expect_s3_class(screen, "SOMGridStabilityScreen")
  expect_setequal(screen$decision_table$Topology, c("hexagonal", "rectangular"))
  expect_false("Classes" %in% names(screen$decision_table))
  expect_false("NodeWeighting" %in% names(screen$decision_table))
  expect_true(all(c("TopologyPreservation_Mean", "NeighborhoodOverlap_Mean",
    "ProjectedDistance_Mean", "ProjectedDistance_P95",
    "EmptyNodes", "SingletonNodes", "NodeUtilization", "SOMDistance_P95") %in%
    names(screen$decision_table)))
  expect_true("MapStabilityPlateau" %in% names(screen$decision_table))
  expect_identical(screen$settings$scoring_cohort, "all_reference_participants")
  expect_identical(screen$settings$screen_scope, "map")
  expect_false(screen$settings$mclust_fitted)
  expect_identical(screen$settings$grid_design_k, 2L)
  expect_true(all(c("ExpectedNodesPerProfile_K2",
    "ExpectedNodesPerMaximumProfile") %in% names(screen$decision_table)))
  expect_true(all(screen$map_replicates$ParticipantsCompared ==
    sum(stats::complete.cases(df_Labelled[, vars_SOM]))))
  expect_error(ConfirmSOMGridStability(screen, 1), "screen_scope = 'joint'")
})

test_that("joint SOM grid screen remains available for explicit confirmation", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  data("SampleVariableTypes", package = "SciDataReportR")
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  vars_SOM <- c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin")
  screen <- suppressWarnings(ScreenSOMGridStability(
    data = df_Labelled, variables = vars_SOM, k_range = 2, models = 1,
    grid_candidates = 3, screen_scope = "joint", screen_resamples = 1,
    stability_cores = 1, lpa_timeout_seconds = NULL
  ))

  expect_setequal(screen$decision_table$NodeWeighting, c("uniform", "occupancy"))
  expect_true(all(screen$decision_table$Classes == 2))
  expect_identical(screen$settings$screen_scope, "joint")

  confirmed <- suppressWarnings(ConfirmSOMGridStability(
    screen, screen$decision_table$.ScreenRow[[1]], confirm_resamples = 1
  ))
  expect_s3_class(confirmed, "SOMGridStabilityConfirmation")
  expect_equal(length(confirmed$models), 1)
  model <- confirmed$models[[1]]
  expect_true(all(model$Stability$replicates$ParticipantsCompared ==
    sum(stats::complete.cases(df_Labelled[, vars_SOM]))))
  expect_true("Flag_InfluenceCandidate" %in% names(model$Influence))
})

test_that("map topology metrics are invariant to coordinate orientation", {
  reference <- dplyr::tibble(Grid_X = c(0, 1, 0, 1), Grid_Y = c(0, 0, 1, 1))
  reflected <- dplyr::mutate(reference, Grid_X = - .data$Grid_X)
  comparison <- SciDataReportR:::.SOMMapComparison(reference, reflected,
    pairs = rbind(c(1L, 2L), c(1L, 3L), c(1L, 4L), c(2L, 3L),
      c(2L, 4L), c(3L, 4L)),
    neighborhood_rows = 1:4, neighborhood_k = 2)
  expect_equal(comparison$topology, 1)
  expect_equal(comparison$neighborhood, 1)
})

test_that("map stability plateau requires topology, distance, and singleton criteria", {
  decision_table <- dplyr::tibble(
    TopologyPreservation_Mean = c(.90, .89, .89, .89),
    ProjectedDistance_Mean = c(1, 1.20, 1.02, 1.02),
    SingletonNodeProportion = c(.04, .04, .06, .04),
    NeighborhoodOverlap_Mean = c(.8, .8, .8, .8))
  flagged <- SciDataReportR:::.SOMMapStabilityPlateau(decision_table)
  expect_identical(flagged$MapStabilityPlateau, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("occupancy weighting uses participant node counts without changing the default", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  vars_SOM <- c("age", "sex", "AXL", "Adiponectin")
  uniform <- suppressWarnings(CreateClusterModel_SOM_MClust(
    SampleData, variables = vars_SOM, method = "finalize", final_k = 2,
    final_model = 1, som_xdim = 3, som_ydim = 3,
    lpa_min_node_n = 1, lpa_timeout_seconds = NULL
  ))
  weighted <- suppressWarnings(CreateClusterModel_SOM_MClust(
    SampleData, variables = vars_SOM, method = "finalize", final_k = 2,
    final_model = 1, som_xdim = 3, som_ydim = 3,
    node_weighting = "occupancy", lpa_timeout_seconds = NULL
  ))

  expect_identical(uniform$ModelInfo_MClust$diagnostics$lpa_preprocess$node_weighting,
    "uniform")
  expect_identical(weighted$ModelInfo_MClust$diagnostics$lpa_preprocess$node_weighting,
    "occupancy")
  expect_equal(uniform$ModelInfo_MClust$diagnostics$lpa_preprocess$effective_node_weight,
    sum(uniform$ModelInfo_SOM$SOMFit$node_occupancy$n >= 1))
  expect_equal(weighted$ModelInfo_MClust$diagnostics$lpa_preprocess$effective_node_weight,
    sum(weighted$ModelInfo_SOM$SOMFit$node_occupancy$n))
  expect_true(all(c("WeightedLogLikelihood", "EffectiveSampleSize",
    "WeightedConverged") %in% names(weighted$ModelInfo_MClust$fit_table)))
  projected <- ProjectCluster(weighted, SampleData[seq_len(10), , drop = FALSE])
  expect_true(all(!is.na(projected$DataWithClusters$Cluster)))
  expect_true(all(projected$ProbFit$individual$Cluster ==
    projected$ProbFit$individual$NodeCluster))
  expect_error(CreateClusterModel_SOM_MClust(
    SampleData, variables = vars_SOM, method = "finalize", final_k = 2,
    final_model = 6, node_weighting = "occupancy", lpa_timeout_seconds = NULL
  ), "could not be estimated")
})

test_that("occupancy-weighted stability is deterministic in serial and parallel", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  skip_if(future::availableCores() < 2L,
    "This test environment permits only one parallel worker.")
  skip_if(isTRUE(pkgload::is_dev_package("SciDataReportR")),
    "Multisession workers load the installed package; this assertion runs after package installation.")

  data("SampleData", package = "SciDataReportR")
  vars_SOM <- c("age", "sex", "AXL", "Adiponectin")
  fit <- function(cores) suppressWarnings(CreateClusterModel_SOM_MClust(
    SampleData, variables = vars_SOM, method = "finalize", final_k = 2,
    final_model = 1, som_xdim = 3, som_ydim = 3,
    node_weighting = "occupancy", stability_resamples = 2,
    stability_seed = 901, stability_cores = cores, lpa_timeout_seconds = NULL
  ))
  serial <- fit(1L)
  parallel <- fit(2L)
  expect_identical(serial$Stability$replicates, parallel$Stability$replicates)
  expect_identical(serial$Stability$summary, parallel$Stability$summary)
})

test_that("SOM influence flags require high distance and low inclusion together", {
  model <- list(
    Specification = list(selected = list(model = 1L, k = 2L)),
    Stability = list(participant_inclusion = dplyr::tibble(
      .row_id = 1:3, Model = 1L, Classes = 2L,
      InclusionProbability = c(0.5, 0.9, 0.5), SuccessfulRefits = 10L)),
    ProbFit = list(individual = dplyr::tibble(
      .row_id = 1:3, SOM_Distance = 1:3,
      SOMDist_percentile_train = c(0.96, 0.96, 0.50),
      SOM_Node = 1:3, Cluster = c(1L, 1L, 2L),
      Top_Distance_Variable = "AXL", Top_Distance_Share = 0.5))
  )
  report <- SciDataReportR:::.SOMInfluenceReport(model)
  expect_identical(report$Flag_InfluenceCandidate, c(TRUE, FALSE, FALSE))
})

test_that("SOM sensitivity refits require explicit exclusions and preserve the source model", {
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))

  data("SampleData", package = "SciDataReportR")
  vars_SOM <- c("age", "sex", "AXL", "Adiponectin")
  model <- suppressWarnings(CreateClusterModel_SOM_MClust(
    SampleData, variables = vars_SOM, method = "finalize", final_k = 2,
    final_model = 1, som_xdim = 3, som_ydim = 3,
    stability_resamples = 1, lpa_timeout_seconds = NULL
  ))
  expect_error(RefitSOMSensitivity(model, SampleData,
    exclusion_reason = "reviewed"), "exclude_row_ids")
  sensitivity <- suppressWarnings(RefitSOMSensitivity(model, SampleData,
    exclude_row_ids = c(1, 2), exclusion_reason = "Reviewer-approved check",
    stability_resamples = 1, lpa_timeout_seconds = NULL
  ))

  expect_null(model$Sensitivity)
  expect_identical(sensitivity$Sensitivity$excluded_row_ids, c(1, 2))
  expect_identical(sensitivity$Sensitivity$exclusion_reason, "Reviewer-approved check")
  expect_equal(sensitivity$Sensitivity$n_excluded, 2)
  expect_identical(sensitivity$Specification$projection$assignment,
    model$Specification$projection$assignment)
  expect_identical(sensitivity$Specification$projection$soft_n_nodes,
    model$Specification$projection$soft_n_nodes)
  expect_identical(sensitivity$Specification$projection$boundary_margin_threshold,
    model$Specification$projection$boundary_margin_threshold)
  expect_identical(sensitivity$Specification$projection$boundary_entropy_threshold,
    model$Specification$projection$boundary_entropy_threshold)
})
