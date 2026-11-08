SuppressOrdinaryWarnings <- function(expr) {
  withCallingHandlers(expr, warning = function(condition) {
    if (!inherits(condition, "lifecycle_warning_deprecated")) invokeRestart("muffleWarning")
  })
}

SOMOutputFixture <- function(...) {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  set.seed(6032)
  df_Test <- EmbeddedAB()
  suppressMessages(suppressWarnings(CreateClusterModel_SOM_MClust(
    df_Test, c("A", "B"), k_range = 2:3, models = 1,
    som_xdim = 4, som_ydim = 4, Relabel = FALSE,
    mixture_bootstrap_n = 0, ...)))
}

test_that("candidate boundary shares use the final soft rule without changing default AHP", {
  model <- SOMOutputFixture()
  fits <- model$ModelInfo_MClust$fit_table
  selected <- model$ModelInfo_MClust$best_fit_name
  ids <- paste0("model_", fits$Model, "_class_", fits$Classes)
  expect_equal(fits$BoundaryPct[match(selected, ids)],
    mean(model$ProbFit$individual$Flag_BoundaryParticipant, na.rm = TRUE))
  expect_equal(fits$ahp_index, rowMeans(fits[c("AIC_scaled", "BIC_scaled", "Entropy_scaled")]))
  expect_true(all(c("equal", "only_AIC", "only_BIC", "only_Entropy",
    "only_Boundary", "with_Boundary", "without_Boundary") %in%
    model$ModelInfo_MClust$ahp_sensitivity$scheme))
  expect_true("BoundaryPct" %in% as.character(model$fit_plot$data$name))
  boundary <- SOMOutputFixture(criteria_weights = c(Boundary = 1))
  expect_equal(boundary$ModelInfo_MClust$fit_table$ahp_index,
    boundary$ModelInfo_MClust$fit_table$Boundary_scaled)
})

test_that("individual node sizes match training occupancy", {
  model <- SOMOutputFixture()
  individual <- model$ProbFit$individual
  counts <- tabulate(individual$SOM_Node, nrow(model$ProbFit$node))
  expect_equal(individual$NodeN, counts[individual$SOM_Node])
  expect_equal(individual$Flag_OnePersonNode, individual$NodeN == 1L)
  expect_equal(model$DataWithClusters$NodeN, individual$NodeN)
  expect_equal(model$DataWithClusters$Flag_OnePersonNode, individual$Flag_OnePersonNode)
  expect_equal(sum(counts), sum(!is.na(individual$NodeN)))
})

test_that("map diagnostics reuse the screen quality and preserve RNG", {
  set.seed(19)
  zmat <- scale(as.matrix(EmbeddedAB()))
  map <- .FitSOMGridScreen(zmat, 4L, 934521L, "hexagonal", "gaussian", 20L)
  before <- .Random.seed
  diagnostic <- .SOMMapDiagnostics(map, zmat, rep(1:2, each = 8), 934521L)
  expect_identical(.Random.seed, before)
  quality <- .SOMGridMapQuality(map, zmat, 3L)$metrics
  expect_equal(unlist(diagnostic$map_quality[names(quality)]), quality)
  expect_equal(nrow(diagnostic$u_matrix), 16L)
  model <- SOMOutputFixture(som_epochs = 20L)
  expect_equal(model$ModelInfo_SOM$training$epochs, 20L)
  z_model <- scale(as.matrix(model$DataWithClusters[, c("A", "B"), drop = FALSE]))
  colnames(z_model) <- colnames(model$ModelInfo_SOM$som_codes)
  expected <- .FitSOMGridScreen(z_model, 4L, 934521L, "hexagonal", "gaussian", 20L)
  expect_equal(model$ModelInfo_SOM$som_codes, expected$codes[[1]], ignore_attr = TRUE)
})

test_that("boundary summaries account for all complete participants", {
  model <- SOMOutputFixture()
  summary <- model$ProbFit$boundary_summary
  expect_equal(sum(summary$N), 320L)
  expect_equal(sum(summary$Share), 1)
  grid <- model$ProbFit$boundary_threshold_grid
  expect_equal(nrow(grid), 15L)
  expect_equal(grid$BoundaryPct[grid$margin == .20 & grid$entropy == .60],
    mean(model$ProbFit$individual$Flag_BoundaryParticipant))
})

test_that("explicit-grid defaults remain identical to ad035a60", {
  model <- SOMOutputFixture()
  expect_equal(as.integer(table(model$ProbFit$individual$Cluster)), c(160L, 160L))
  expect_equal(round(model$ModelInfo_MClust$fit_table$ahp_index, 6), c(.707107, -.707107))
  expect_identical(model$ModelInfo_MClust$best_fit_name, "model_1_class_2")
})

test_that("Mclust names resolve case-insensitively", {
  expect_identical(.ResolveMclustModels(c("eei", "VVI", "Eee", "vvv")), c(1L, 2L, 3L, 6L))
  expect_identical(.ResolveMclustModels("VVV"), .ResolveMclustModels(6))
  expect_error(.ResolveMclustModels("unknown"), "Mclust names")
})

test_that("deprecated SOM arguments warn and reproduce the new names", {
  rlang::local_options(lifecycle_verbosity = "warning")
  set.seed(6032)
  df_Test <- EmbeddedAB()
  arguments <- list(data = df_Test, variables = c("A", "B"), k_range = 2:3,
    models = 1, som_xdim = 4, som_ydim = 4, Relabel = FALSE,
    candidate_min_nodes_per_cluster = NULL)
  pairs <- c(seed_lpa = "seed_mixture", stability_seed = "seed_stability",
    lpa_em_itmax = "mixture_em_itmax", lpa_em_tol = "mixture_em_tol",
    lpa_timeout_seconds = "mixture_timeout_seconds", lpa_progress = "mixture_progress",
    lpa_bootstrap_n = "mixture_bootstrap_n", lpa_drop_zero_sd = "mixture_drop_zero_sd",
    lpa_zero_sd_tol = "mixture_zero_sd_tol", lpa_min_node_n = "mixture_min_node_people",
    min_nodes_per_cluster = "candidate_min_nodes_per_cluster",
    high_dist_quantile = "far_distance_quantile", low_prob_threshold = "clear_membership_threshold")
  reference <- suppressWarnings(suppressMessages(do.call(CreateClusterModel_SOM_MClust, arguments)))
  for (old in names(pairs)) {
    new <- pairs[[old]]
    value <- eval(formals(CreateClusterModel_SOM_MClust)[[new]])
    call_args <- arguments[setdiff(names(arguments), new)]
    call_args[old] <- list(value)
    # Retain the intentional grid/support warnings alongside the expected shim.
    expect_warning(result <- SuppressOrdinaryWarnings(suppressMessages(do.call(CreateClusterModel_SOM_MClust, call_args))),
      class = "lifecycle_warning_deprecated")
    expect_identical(result$ProbFit$individual$Cluster, reference$ProbFit$individual$Cluster)
    expect_identical(result$ModelInfo_MClust$fit_table$ahp_index,
      reference$ModelInfo_MClust$fit_table$ahp_index)
  }
})

test_that("projection thresholds and retired controls retain warning shims", {
  rlang::local_options(lifecycle_verbosity = "warning")
  model <- SOMOutputFixture()
  data <- EmbeddedAB("Projection")
  expected <- suppressWarnings(suppressMessages(ProjectCluster(model, data,
    far_distance_quantile = .90, clear_membership_threshold = .65)))
  expect_warning(actual <- SuppressOrdinaryWarnings(suppressMessages(ProjectCluster(model, data,
    high_dist_quantile = .90, clear_membership_threshold = .65))),
    class = "lifecycle_warning_deprecated")
  expect_equal(actual$ProbFit$individual, expected$ProbFit$individual)
  expect_warning(actual <- SuppressOrdinaryWarnings(suppressMessages(ProjectCluster(model, data,
    far_distance_quantile = .90, low_prob_threshold = .65))),
    class = "lifecycle_warning_deprecated")
  expect_equal(actual$ProbFit$individual, expected$ProbFit$individual)
  expect_warning(.ClusterMethod("explore"), class = "lifecycle_warning_deprecated")
  for (argument in c("min_nodes_per_parameter", "min_participants_per_node")) {
    expect_warning(expect_error(do.call(CreateClusterModel_SOM_MClust,
      c(list(data = data, variables = c("A", "B"), som_epochs = 0),
        stats::setNames(list(3), argument))), "som_epochs"), class = "lifecycle_warning_deprecated")
  }
})

test_that("seed stability is deterministic and combines types by mean ARI", {
  first <- SOMOutputFixture(stability_resamples = 2L, stability_type = "seed",
    stability_cores = 1L)
  second <- SOMOutputFixture(stability_resamples = 2L, stability_type = "seed",
    stability_cores = 1L)
  expect_equal(first$Stability$replicates, second$Stability$replicates)
  expect_true(all(first$Stability$replicates$RetainedParticipants == 320L))
  expect_true(all(first$Stability$replicates$StabilityType == "seed"))
  expect_length(unique(first$Stability$replicates$ModelSeed), 2L)
  expect_true(all(first$Stability$summary$StabilitySuccessRate == 1))
  expect_true("MeanJaccard" %in% names(first$Stability$cluster_recovery_summary))
  both <- SOMOutputFixture(stability_resamples = 2L,
    stability_type = c("subsample", "seed"), stability_cores = 1L)
  expect_setequal(both$Stability$summary$StabilityType, c("subsample", "seed"))
  summary <- both$Stability$summary %>% dplyr::group_by(Model, Classes) %>%
    dplyr::summarise(Expected = mean(StabilityARI_Mean), .groups = "drop")
  actual <- both$ModelInfo_MClust$fit_table %>%
    dplyr::left_join(summary, by = c("Model", "Classes"))
  expect_equal(actual$ReproducibilityScore, actual$Expected)
  expect_false(any(grepl("ARI.*P05", names(actual))))
})

test_that("Gaussian soft assignment is normalized and retained for projection", {
  inverse <- SOMOutputFixture(soft_kernel = "inverse")
  expect_equal(as.integer(table(inverse$ProbFit$individual$Cluster)), c(160L, 160L))
  gaussian <- SOMOutputFixture(soft_kernel = "gaussian")
  expect_gt(gaussian$Specification$projection$soft_bandwidth, 0)
  probabilities <- gaussian$ProbFit$individual %>% dplyr::select(dplyr::starts_with("SoftProb_"))
  expect_equal(rowSums(probabilities), rep(1, 320))
  fits <- gaussian$ModelInfo_MClust$fit_table
  selected <- paste0("model_", fits$Model, "_class_", fits$Classes) == gaussian$ModelInfo_MClust$best_fit_name
  expect_equal(fits$BoundaryPct[selected], mean(gaussian$ProbFit$individual$Flag_BoundaryParticipant))
  projection <- suppressWarnings(suppressMessages(ProjectCluster(gaussian,
    EmbeddedAB("Projection"))))
  z_projected <- ProjectZScore(EmbeddedAB("Projection"), c("A", "B"),
    parameters = gaussian$ZScoreObject, ParameterInputType = "ZScoreObj", names_prefix = "Z_")$ZScores
  expected <- .SOMSoftProjection(as.matrix(z_projected), gaussian$ModelInfo_SOM$som_codes,
    gaussian$ProbFit$node, soft_kernel = "gaussian",
    soft_bandwidth = gaussian$Specification$projection$soft_bandwidth)
  expect_equal(projection$ProbFit$individual$SoftMargin, expected$margin, tolerance = 1e-8)
})

test_that("grid refits are identical in serial and parallel and failures are counted", {
  skip_if_not_installed("kohonen")
  set.seed(56)
  zmat <- scale(as.matrix(EmbeddedAB()))
  settings <- list(zmat = zmat, sides = c(3L, 4L), refits = 2L,
    seed_som = 1234L, pair_sample = 40L, trust_sample = 30L, som_epochs = 10L)
  serial <- suppressWarnings(do.call(.ScreenSOMGridMatrix, c(settings, list(cores = 1L))))
  parallel <- suppressWarnings(do.call(.ScreenSOMGridMatrix, c(settings, list(cores = 2L))))
  expect_equal(serial$table, parallel$table)
  expect_equal(serial$refits, parallel$refits)
  expect_equal(serial$pairs, parallel$pairs)
  expect_true(all(serial$table$SuccessRate_seed == 1))
  original <- .FitSOMGridScreen
  local_mocked_bindings(.FitSOMGridScreen = function(zmat, side, seed, ...) {
    if (seed == 1235L) stop("planted refit failure")
    original(zmat, side, seed, ...)
  }, .package = "SciDataReportR")
  partial <- suppressWarnings(do.call(.ScreenSOMGridMatrix, settings))
  expect_true(all(partial$table$SuccessRate_seed == .5))
  expect_equal(sum(partial$refit_status$Status == "failed"), 2L)
})

test_that("finalize failures name the requested model and k", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  local_mocked_bindings(.FitTidyLPAMclust = function(...) stop("planted mixture failure"),
    .package = "SciDataReportR")
  set.seed(19)
  df_Test <- EmbeddedAB()
  expect_error(suppressWarnings(suppressMessages(CreateClusterModel_SOM_MClust(
    df_Test, c("A", "B"), method = "finalize", final_model = "VVV", final_k = 2,
    som_xdim = 4, som_ydim = 4, mixture_bootstrap_n = 0))), "VVV, k = 2")
})

test_that("all failed stability refits retain a zero success rate", {
  original <- .FitTidyLPAMclust
  fits <- 0L
  local_mocked_bindings(.FitTidyLPAMclust = function(...) {
    fits <<- fits + 1L
    if (fits > 2L) stop("planted stability failure")
    original(...)
  }, .package = "SciDataReportR")
  model <- SOMOutputFixture(stability_type = "seed", stability_resamples = 1L,
    stability_cores = 1L)
  expect_true(all(model$Stability$summary$StabilitySuccessRate == 0))
  expect_true(all(model$Stability$replicates$Status != "success"))
  expect_equal(nrow(model$Stability$cluster_recovery_summary), 0L)
})

test_that("sensitivity keeps schemes with no remaining positive weights", {
  df_Fits <- tibble::tibble(Model = 1L, Classes = 2L, BoundaryPct = .1,
    Boundary_scaled = 0, Eligible = TRUE)
  schemes <- .SOMAHPWeightSensitivity(df_Fits, c(Boundary = 1), "model_1_class_2")
  expect_true("without_Boundary" %in% schemes$scheme)
  expect_true(is.na(schemes$top_model[schemes$scheme == "without_Boundary"]))
})
