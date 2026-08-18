test_that("numeric clustering models fit and project with frozen preprocessing", {
  skip_if_not_installed("mclust")
  set.seed(101)
  df_Test <- data.frame(
    x = c(rnorm(20), rnorm(20, 4), NA_real_),
    y = c(rnorm(20), rnorm(20, 4), NA_real_)
  )

  mclust_model <- CreateClusterModel_MClust(
    data = df_Test, variables = c("x", "y"), method = "finalize",
    final_k = 2, final_model = 1
  )
  mclust_projection <- ProjectCluster(mclust_model, df_Test)
  expect_equal(nrow(mclust_projection$DataWithClusters), nrow(df_Test))
  expect_true(is.na(mclust_projection$DataWithClusters$Cluster[nrow(df_Test)]))
  expect_equal(mclust_model$DataWithClusters$Cluster[-nrow(df_Test)], mclust_projection$DataWithClusters$Cluster[-nrow(df_Test)])

  kmeans_model <- CreateClusterModel_KMeans(
    data = df_Test, variables = c("x", "y"), method = "finalize", final_k = 2,
    ZScoreType = "None"
  )
  expect_equal(kmeans_model$Preprocessing$Scaling, "None")
  expect_equal(
    kmeans_model$DataWithClusters$Cluster[-nrow(df_Test)],
    ProjectCluster(kmeans_model, df_Test)$DataWithClusters$Cluster[-nrow(df_Test)]
  )

  kmeans_projection <- ProjectCluster(kmeans_model, df_Test)
  # Figures live beside the object they describe, never in a flat list.
  expect_null(kmeans_model$plots)
  expect_null(kmeans_projection$plots)
  expect_s3_class(kmeans_model$fit_plot, "ggplot")
  expect_true(all(c("elbow", "silhouette", "profiles") %in%
    names(kmeans_model$ModelInfo$plots)))
  # Cluster occupancy was dropped from the pipeline figures.
  expect_false("occupancy" %in% names(kmeans_model$ModelInfo$plots))
  expect_true("distance_hist" %in% names(kmeans_model$ModelInfo$FitDiagnostics$plots))
  expect_true(all(c("projection_fit_class_bar", "separation_map") %in%
    names(kmeans_projection$ProjectionFit$plots)))
  expect_null(kmeans_projection$ModelInfo$plots$fit_plot)
})

test_that("ProjectCluster dispatches by fitted pipeline class", {
  skip_if_not_installed("mclust")
  set.seed(111)
  df_Test <- data.frame(x = c(rnorm(20), rnorm(20, 4)),
    y = c(rnorm(20), rnorm(20, 4)))
  model <- CreateClusterModel_MClust(
    df_Test, c("x", "y"), method = "finalize", final_k = 2,
    final_model = 1)
  expect_s3_class(model, "Pipeline_MClust")
  expect_equal(nrow(ProjectCluster(model, df_Test)$DataWithClusters), nrow(df_Test))
  expect_error(ProjectCluster(list(), df_Test), "not a supported finalized")
})

test_that("numeric Mclust model IDs and preprocessing aliases are standardized", {
  skip_if_not_installed("mclust")
  set.seed(112)
  df_Test <- data.frame(x = c(rnorm(20), rnorm(20, 4)),
    y = c(rnorm(20), rnorm(20, 4)))
  model <- CreateClusterModel_MClust(
    df_Test, c("x", "y"), method = "finalize", final_k = 2,
    final_model = 1, ZScoreType = "None")
  expect_identical(model$ModelInfo$final_model, 1L)
  expect_identical(model$ModelInfo$final_model_name, "EEI")
  expect_true(all(c("Model", "ModelName") %in% names(model$ModelInfo$fit_table)))
  expect_identical(model$ModelInfo$fit_table$ModelName[[1]], "EEI")
  expect_identical(model$Preprocessing$ZScoreType, "None")
  expect_error(
    CreateClusterModel_MClust(df_Test, c("x", "y"), method = "finalize",
      final_k = 2, final_model = "EEI"),
    "integer model IDs"
  )
  expect_error(
    CreateClusterModel_KMeans(df_Test, c("x", "y"), method = "finalize",
      final_k = 2, ZScoreType = "None", Scaling = "Center and Scale"),
    "different preprocessing"
  )
  expect_silent(CreateClusterModel_KMeans(
    df_Test, c("x", "y"), method = "finalize", final_k = 2,
    ZScoreType = "None", Scaling = "None"))
})

test_that("established SOM projection aliases delegate to ProjectCluster", {
  assign("ProjectCluster.Pipeline_AliasTest",
    function(object, new_df, ...) "projected", envir = .GlobalEnv)
  on.exit(rm(ProjectCluster.Pipeline_AliasTest, envir = .GlobalEnv), add = TRUE)
  object <- structure(list(), class = "Pipeline_AliasTest")
  expect_identical(Project_SOMClust(object, data.frame(x = 1)), "projected")
  expect_warning(
    expect_identical(ProjectSOMCluster(object, data.frame(x = 1)), "projected"),
    "deprecated"
  )
})

test_that("PlotClusterComposition supports both clinically useful facet layouts", {
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")

  by_variable <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster)
  by_cluster <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster,
    facet_by = "cluster")
  enrichment <- PlotClusterComposition(
    df_Test, c("CatVar1", "CatVar2"), df_Test$TruthCluster,
    style = "enrichment")

  expect_s3_class(by_variable, "ggplot")
  expect_s3_class(by_cluster, "ggplot")
  expect_s3_class(enrichment, "ggplot")
  expect_identical(
    rlang::as_string(rlang::get_expr(by_variable$facet$params$facets[[1]])),
    "Variable"
  )
  expect_identical(
    rlang::as_string(rlang::get_expr(by_cluster$facet$params$facets[[1]])),
    "Cluster"
  )
  expect_error(
    PlotClusterComposition(df_Test, "CatVar1", df_Test$TruthCluster,
      facet_by = "outcome"),
    "should be one of"
  )
})

test_that("clustering defaults review the full 2:10 candidate range", {
  DefaultRange <- function(fun, argument) eval(formals(fun)[[argument]])
  expect_identical(DefaultRange(CreateClusterModel_MClust, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_KMeans, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_PCA_MClust, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_PCA_KMeans, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_Gower_PAM, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_LatentClass, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_MCA_MClust, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_SOM_MClust, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_HDBSCAN, "minPts_range"), 2:10)
  expect_identical(DefaultRange(CreateClusterModel_SOM_HDBSCAN, "minPts_range"), 2:10)
})

test_that("categorical enrichment identifies a defining category", {
  df_Test <- data.frame(
    Sex = factor(c(rep("Male", 9), "Female", rep("Male", 2), rep("Female", 8))),
    Cluster = rep(c(1L, 2L), each = 10L))
  plot <- PlotClusterComposition(df_Test, "Sex", df_Test$Cluster,
    style = "enrichment")
  expect_true(any(plot$data$Enrichment > 0))
  expect_true(any(plot$data$Category == "Sex: Male" & plot$data$Enrichment > 0))
})

test_that("clustering figures are ggplots and are never aliased or NULL", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("cluster")
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
  df_New <- subset(SimulatedPhenotypeData, Cohort == "Projection")[1:60, ]
  model <- CreateClusterModel_KMeans(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    nstart = 5, stability_resamples = 1)
  projection <- ProjectCluster(model, df_New)

  blocks <- list(
    model$ModelInfo$plots, model$ModelInfo$FitDiagnostics$plots,
    model$ProbFit$plots, model$Stability$plots,
    projection$ProjectionFit$plots, projection$ProbFit$plots)
  for (block in blocks) {
    expect_gt(length(block), 0)
    expect_false(anyDuplicated(names(block)) > 0)
    expect_false(any(vapply(block, is.null, logical(1))))
    expect_true(all(vapply(block, inherits, logical(1), what = "ggplot")))
    expect_false(any(grepl("(MaxProb|Posterior).*Boxplot", names(block))))
    # No two entries in a block may be the same object under different names.
    for (i in seq_along(block)) for (j in seq_along(block)) {
      if (i < j) expect_false(identical(block[[i]], block[[j]]))
    }
  }
})

test_that("PCA, HDBSCAN, and Gower PAM pipelines retain projectable rows", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("dbscan")
  set.seed(102)
  df_Numeric <- data.frame(a = c(rnorm(20), rnorm(20, 3), NA), b = c(rnorm(20), rnorm(20, 3), NA), c = c(rnorm(20), rnorm(20, 3), NA))
  pca_model <- CreateClusterModel_PCA_KMeans(df_Numeric, c("a", "b", "c"), method = "finalize", final_k = 2)
  expect_equal(nrow(ProjectCluster(pca_model, df_Numeric)$DataWithClusters), nrow(df_Numeric))
  hdbscan_model <- CreateClusterModel_HDBSCAN(df_Numeric, c("a", "b", "c"),
    method = "finalize", final_minPts = 5, final_cluster_selection_epsilon = 0)
  expect_equal(nrow(ProjectCluster(hdbscan_model, df_Numeric)$DataWithClusters), nrow(df_Numeric))

  df_Mixed <- data.frame(age = c(20, 22, 65, 67, NA_real_), group = factor(c("a", "a", "b", "b", "a")))
  pam_model <- CreateClusterModel_Gower_PAM(df_Mixed, c("age", "group"),
    method = "finalize", final_k = 2)
  expect_equal(nrow(ProjectCluster(pam_model, df_Mixed)$DataWithClusters), nrow(df_Mixed))
})

test_that("categorical clustering models project into their frozen class spaces", {
  skip_if_not_installed("poLCA")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("mclust")
  set.seed(103)
  df_Categorical <- data.frame(
    symptom = factor(sample(c("none", "mild", "severe"), 60, replace = TRUE)),
    diagnosis = factor(sample(c("control", "case"), 60, replace = TRUE))
  )
  lca_model <- CreateClusterModel_LatentClass(df_Categorical, c("symptom", "diagnosis"), final_k = 2, nrep = 2)
  expect_equal(nrow(ProjectCluster(lca_model, df_Categorical)$DataWithClusters), nrow(df_Categorical))
  mca_model <- CreateClusterModel_MCA_MClust(df_Categorical, c("symptom", "diagnosis"), method = "finalize", final_k = 2, final_model = 1)
  expect_equal(nrow(ProjectCluster(mca_model, df_Categorical)$DataWithClusters), nrow(df_Categorical))
})

test_that("exploratory clustering models return AHP and deterministic stability", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("cluster")
  data("SimulatedPhenotypeData")
  df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:120, ]

  first <- CreateClusterModel_KMeans(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
  second <- CreateClusterModel_KMeans(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
  expect_equal(first$Stability$replicates, second$Stability$replicates)
  expect_equal(sum(first$ModelInfo$fit_table$Recommended), 1)
  expect_true(all(c("ReproducibilityScore", "ahp_index") %in% names(first$ModelInfo$fit_table)))
  expect_true(all(c("VI", "NMI", "FowlkesMallows") %in% names(first$Stability$replicates)))
  expect_true(all(c("participant_inclusion", "cluster_inclusion", "coassignment") %in% names(first$Stability)))
  expect_true(all(vapply(first$Stability$coassignment, function(x) {
    identical(x$status, "available")
  }, logical(1))))
  expect_false("bootstrap_agreement" %in% names(first$Stability$plots))
  expect_s3_class(first$Stability$plots$cluster_recovery, "ggplot")
  expect_s3_class(first$Stability$plots$cluster_recovery, "ggplot")

  pam <- CreateClusterModel_Gower_PAM(df_Test, c("Var1", "Var2", "CatVar1"),
    k_range = 2:3, stability_resamples = 2)
  expect_equal(sum(pam$ModelInfo$fit_table$Recommended), 1)
  expect_true(nrow(pam$Stability$cluster_recovery) > 0)
})

test_that("stability partition metrics and coassignment guard have expected behaviour", {
  metrics_same <- .ClusterPartitionMetrics(c(1, 1, 2, 2), c(2, 2, 1, 1))
  metrics_diff <- .ClusterPartitionMetrics(c(1, 1, 2, 2), c(1, 2, 1, 2))
  expect_equal(metrics_same[["VI"]], 0)
  expect_equal(metrics_same[["NMI"]], 1)
  expect_equal(metrics_same[["FowlkesMallows"]], 1)
  expect_gt(metrics_diff[["VI"]], metrics_same[["VI"]])
  reference <- rep(c(1L, 2L), each = 1001)
  diagnostics <- .ClusterStabilityDiagnostics(
    reference, list(reference), seq_along(reference), coassignment_limit = 2000L)
  expect_identical(diagnostics$coassignment$status, "skipped")
})

test_that("HDBSCAN stability reports noise recovery and handles candidate grids", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:160, ]
  model <- CreateClusterModel_HDBSCAN(df_Test, c("DensityX", "DensityY"),
    minPts_range = c(5, 8), cluster_selection_epsilon_range = c(0, 0.05),
    stability_resamples = 2)
  expect_equal(nrow(model$ModelInfo$fit_table), 4)
  expect_true(all(c("NoiseSensitivity", "NoiseSpecificity") %in%
    names(model$Stability$summary)))
  expect_equal(sum(model$ModelInfo$fit_table$Recommended), 1)
})

test_that("every non-SOM pipeline exposes full-pipeline bootstrap stability", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("cluster")
  skip_if_not_installed("dbscan")
  skip_if_not_installed("poLCA")
  skip_if_not_installed("FactoMineR")
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
  vars_Numeric <- paste0("Var", 1:12)
  vars_Categorical <- paste0("CatVar", 1:3)

  models <- list(
    mclust = CreateClusterModel_MClust(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      final_model = 1, stability_resamples = 1),
    kmeans = CreateClusterModel_KMeans(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      nstart = 5, stability_resamples = 1),
    pca_mclust = CreateClusterModel_PCA_MClust(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      final_model = 1, stability_resamples = 1),
    pca_kmeans = CreateClusterModel_PCA_KMeans(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      nstart = 5, stability_resamples = 1),
    hdbscan = CreateClusterModel_HDBSCAN(
      df_Test, c("DensityX", "DensityY"), method = "finalize",
      final_minPts = 8, final_cluster_selection_epsilon = 0,
      stability_resamples = 1),
    gower = CreateClusterModel_Gower_PAM(
      df_Test, c("Var1", "Var2", vars_Categorical), method = "finalize",
      final_k = 4, stability_resamples = 1),
    lca = CreateClusterModel_LatentClass(
      df_Test, vars_Categorical, method = "finalize", final_k = 4,
      nrep = 2, stability_resamples = 1),
    mca_mclust = CreateClusterModel_MCA_MClust(
      df_Test, vars_Categorical, method = "finalize", final_k = 4,
      final_model = 1, stability_resamples = 1)
  )

  # Every method carries the same structural blocks, populated with figures
  # appropriate to that method rather than a shared lowest common denominator.
  for (model in models) {
    expect_identical(model$Stability$settings$refit_scope, "full_pipeline")
    expect_equal(nrow(model$Stability$replicates), 1)
    expect_null(model$plots)
    expect_s3_class(model$fit_plot, "ggplot")
    expect_false("occupancy" %in% names(model$ModelInfo$plots))
    expect_gt(length(model$ModelInfo$plots), 3)
    expect_true("distance_hist" %in% names(model$ModelInfo$FitDiagnostics$plots))
    expect_gt(length(model$ProbFit$plots), 0)
    expect_false(any(grepl("(MaxProb|Posterior).*Boxplot", names(model$ProbFit$plots))))
    expect_false("bootstrap_agreement" %in% names(model$Stability$plots))
    expect_s3_class(model$Stability$plots$cluster_recovery, "ggplot")
  }
  # Method-specific figures a reader would expect from each approach.
  expect_true("silhouette" %in% names(models$kmeans$ModelInfo$plots))
  expect_true("silhouette" %in% names(models$gower$ModelInfo$plots))
  expect_true("bic" %in% names(models$mclust$ModelInfo$plots))
  expect_false("classification_uncertainty" %in% names(models$mclust$ModelInfo$plots))
  expect_true("cluster_persistence" %in% names(models$hdbscan$ModelInfo$plots))
  expect_true("response_probabilities" %in% names(models$lca$ModelInfo$plots))
  expect_true("scree" %in% names(models$pca_kmeans$ModelInfo$plots))
  expect_true("scree" %in% names(models$mca_mclust$ModelInfo$plots))
  expect_gte(ncol(models$pca_kmeans$PCA$Scores), 4)
})

test_that("SOM HDBSCAN freezes node-level assignments for projection", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  skip_if_not_installed("dbscan")
  df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training") %>%
    dplyr::slice(1:40)
  model <- CreateClusterModel_SOM_HDBSCAN(
    df_Training, variables = paste0("Var", 1:4), method = "finalize",
    final_minPts = 3L, final_cluster_selection_epsilon = 0,
    som_xdim = 3, som_ydim = 3,
    models = 1, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30,
    Relabel = FALSE)
  projected <- ProjectCluster(model, df_Training[1:5, ])
  expect_s3_class(model, "Pipeline_SOM_HDBSCAN")
  expect_true(all(c("ModelInfo_SOM", "ModelInfo_HDBSCAN", "Specification") %in%
    names(model)))
  expect_s3_class(projected, "Project_SOMHDBSCAN")
  expect_true(all(c("individual", "summary", "by_cluster", "out_of_support",
    "policy", "plots") %in% names(projected$ProjectionFit)))
})

test_that("SOM HDBSCAN exploratory models are not projectable", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("tidyLPA")
  skip_if_not_installed("dbscan")
  df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training") %>%
    dplyr::slice(1:40)
  review <- CreateClusterModel_SOM_HDBSCAN(
    df_Training, variables = paste0("Var", 1:4), minPts_range = 3L,
    cluster_selection_epsilon_range = 0, som_xdim = 3, som_ydim = 3,
    models = 1, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30,
    Relabel = FALSE)
  expect_s3_class(review, "Pipeline_SOM_HDBSCAN_Exploratory")
  expect_null(review$ProbFit)
  expect_error(ProjectCluster(review, df_Training[1:5, ]), "requires a finalized")
})

test_that("projections triage cases against the frozen training reference", {
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
  df_New <- subset(SimulatedPhenotypeData, Cohort == "Projection")[1:80, ]
  model <- CreateClusterModel_MClust(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    final_model = 1)
  projection <- ProjectCluster(model, df_New)

  expect_true(is.finite(
    model$ModelInfo$FitDiagnostics$distance_reference$high_distance_cutoff))
  expect_s3_class(projection$ProjectionFit$fit_class, "factor")
  expect_identical(
    levels(projection$ProjectionFit$fit_class),
    c("Good fit", "Uncertain membership", "Poor fit to training structure",
      "Potential novel phenotype"))
  expect_equal(
    nrow(projection$ProjectionFit$overall_fit_summary), 1)
  expect_true("Projection_Fit_Class" %in% names(projection$DataWithClusters))
  # The cutoff is the training one, not one recomputed from the new cohort.
  expect_equal(
    projection$ProjectionFit$overall_fit_summary$high_distance_cutoff,
    model$ModelInfo$FitDiagnostics$distance_reference$high_distance_cutoff)
})

test_that("SOM finalized models can report full-pipeline stability", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:100, ]
  model <- CreateClusterModel_SOM_MClust(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    final_model = 1, som_xdim = 4, som_ydim = 4,
    stability_resamples = 1, Relabel = FALSE,
    min_nodes_per_cluster = NULL, lpa_timeout_seconds = NULL)
  expect_equal(nrow(model$Stability$summary), 1)
  expect_identical(model$Stability$settings$refit_scope, "full_pipeline")
  expect_identical(model$method, "finalize")
  # The SOM pipeline keeps its own figures where they were; it gains stability
  # figures, and never a flat plots list.
  expect_null(model$plots)
  expect_false("bootstrap_agreement" %in% names(model$Stability$plots))
  expect_s3_class(model$Stability$plots$cluster_recovery, "ggplot")
  expect_s3_class(model$ModelInfo_SOM$plots$Circular, "htmlwidget")
  expect_false("occupancy" %in% names(model$ModelInfo_SOM$plots))
})

test_that("stability validation and failed resamples are explicit", {
  expect_error(
    CreateClusterModel_KMeans(
      data.frame(x = 1:5, y = 6:10), c("x", "y"),
      stability_resamples = -1),
    "non-negative integer")
  failed <- .ClusterBootstrapStability(
    data.frame(x = 1:6), reference = rep(1:2, each = 3),
    fit_project = function(boot, original) stop("expected failure"),
    resamples = 2, seed = 1)
  expect_equal(nrow(failed$failures), 2)
  expect_true(all(failed$replicates$Status == "failed"))
  expect_equal(failed$summary$StabilitySuccessRate, 0)
})

test_that("all-noise HDBSCAN candidates remain projectable", {
  skip_if_not_installed("dbscan")
  set.seed(290)
  df_Test <- data.frame(x = rnorm(20), y = rnorm(20))
  model <- CreateClusterModel_HDBSCAN(
    df_Test, c("x", "y"), method = "finalize", final_minPts = 20,
    final_cluster_selection_epsilon = 0, stability_resamples = 1)
  projected <- ProjectCluster(model, df_Test)
  expect_true(all(model$ProbFit$individual$Cluster == 0))
  expect_true(all(projected$ProbFit$individual$Cluster == 0))
  expect_true(all(!projected$ProbFit$individual$InTrainingSupport))
})

test_that("final clustering outputs use clean canonical fields", {
  skip_if_not_installed("mclust")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:80, ]
  df_Test$Var1[[3]] <- NA_real_
  model <- CreateClusterModel_MClust(
    df_Test, paste0("Var", 1:4), method = "finalize", final_k = 2,
    final_model = 1, ClusterVariableName = "Phenotype")
  expect_identical(nrow(model$DataWithClusters), nrow(df_Test))
  expect_identical(names(model$DataWithClusters), c(names(df_Test), "Phenotype"))
  expect_true(is.na(model$DataWithClusters$Phenotype[[3]]))
  expect_false(any(c("df_with_clusters", "ClusterName", "complete_rows",
    "CandidateAudit", "ModelInfo_Mclust") %in% names(model)))
  expect_identical(model$ClusterVariableName, "Phenotype")
  expect_true("ModelInfo_MClust" %in% names(model))
  expect_false(".scidr_rowid" %in% names(model$ProbFit$individual))
})
