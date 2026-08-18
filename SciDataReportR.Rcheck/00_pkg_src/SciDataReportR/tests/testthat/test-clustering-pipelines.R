test_that("numeric clustering models fit and project with frozen preprocessing", {
  skip_if_not_installed("mclust")
  set.seed(101)
  df_Test <- data.frame(
    x = c(rnorm(20), rnorm(20, 4), NA_real_),
    y = c(rnorm(20), rnorm(20, 4), NA_real_)
  )

  mclust_model <- CreateMclustClusterModel(
    data = df_Test, variables = c("x", "y"), method = "finalize",
    final_k = 2, final_model = "EEI"
  )
  mclust_projection <- ProjectMclustCluster(mclust_model, df_Test)
  expect_equal(nrow(mclust_projection$df_with_clusters), nrow(df_Test))
  expect_true(is.na(mclust_projection$df_with_clusters$Cluster[nrow(df_Test)]))
  expect_equal(mclust_model$df_with_clusters$Cluster[-nrow(df_Test)], mclust_projection$df_with_clusters$Cluster[-nrow(df_Test)])

  kmeans_model <- CreateKMeansClusterModel(
    data = df_Test, variables = c("x", "y"), method = "finalize", final_k = 2,
    Scaling = "None"
  )
  expect_equal(kmeans_model$Preprocessing$Scaling, "None")
  expect_equal(
    kmeans_model$df_with_clusters$Cluster[-nrow(df_Test)],
    ProjectKMeansCluster(kmeans_model, df_Test)$df_with_clusters$Cluster[-nrow(df_Test)]
  )

  kmeans_projection <- ProjectKMeansCluster(kmeans_model, df_Test)
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
  expect_identical(DefaultRange(CreateMclustClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateKMeansClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreatePCAMclustClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreatePCAKMeansClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateGowerPAMClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateLatentClassClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateMCAMclustClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateSOMMclustClusterModel, "k_range"), 2:10)
  expect_identical(DefaultRange(CreateHDBSCANClusterModel, "minPts_range"), 2:10)
  expect_identical(DefaultRange(CreateSOMHDBSCANClusterModel, "minPts_range"), 2:10)
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
  model <- CreateKMeansClusterModel(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    nstart = 5, stability_resamples = 1)
  projection <- ProjectKMeansCluster(model, df_New)

  blocks <- list(
    model$ModelInfo$plots, model$ModelInfo$FitDiagnostics$plots,
    model$ProbFit$plots, model$Stability$plots,
    projection$ProjectionFit$plots, projection$ProbFit$plots)
  for (block in blocks) {
    expect_gt(length(block), 0)
    expect_false(anyDuplicated(names(block)) > 0)
    expect_false(any(vapply(block, is.null, logical(1))))
    expect_true(all(vapply(block, inherits, logical(1), what = "ggplot")))
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
  pca_model <- CreatePCAKMeansClusterModel(df_Numeric, c("a", "b", "c"), method = "finalize", final_k = 2)
  expect_equal(nrow(ProjectPCAKMeansCluster(pca_model, df_Numeric)$df_with_clusters), nrow(df_Numeric))
  hdbscan_model <- CreateHDBSCANClusterModel(df_Numeric, c("a", "b", "c"), minPts = 5)
  expect_equal(nrow(ProjectHDBSCANCluster(hdbscan_model, df_Numeric)$df_with_clusters), nrow(df_Numeric))

  df_Mixed <- data.frame(age = c(20, 22, 65, 67, NA_real_), group = factor(c("a", "a", "b", "b", "a")))
  pam_model <- CreateGowerPAMClusterModel(df_Mixed, c("age", "group"), k = 2)
  expect_equal(nrow(ProjectGowerPAMCluster(pam_model, df_Mixed)$df_with_clusters), nrow(df_Mixed))
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
  lca_model <- CreateLatentClassClusterModel(df_Categorical, c("symptom", "diagnosis"), final_k = 2, nrep = 2)
  expect_equal(nrow(ProjectLatentClassCluster(lca_model, df_Categorical)$df_with_clusters), nrow(df_Categorical))
  mca_model <- CreateMCAMclustClusterModel(df_Categorical, c("symptom", "diagnosis"), method = "finalize", final_k = 2, final_model = "EEI")
  expect_equal(nrow(ProjectMCAMclustCluster(mca_model, df_Categorical)$df_with_clusters), nrow(df_Categorical))
})

test_that("exploratory clustering models return AHP and deterministic stability", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("cluster")
  data("SimulatedPhenotypeData")
  df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:120, ]

  first <- CreateKMeansClusterModel(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
  second <- CreateKMeansClusterModel(df_Test, paste0("Var", 1:6),
    k_range = 2:3, nstart = 5, stability_resamples = 2, stability_seed = 44)
  expect_equal(first$Stability$replicates, second$Stability$replicates)
  expect_equal(sum(first$ModelInfo$fit_table$Recommended), 1)
  expect_true(all(c("ReproducibilityScore", "ahp_index") %in% names(first$ModelInfo$fit_table)))
  expect_false("bootstrap_agreement" %in% names(first$Stability$plots))
  expect_s3_class(first$Stability$plots$cluster_recovery, "ggplot")
  expect_s3_class(first$Stability$plots$cluster_recovery, "ggplot")

  pam <- CreateGowerPAMClusterModel(df_Test, c("Var1", "Var2", "CatVar1"),
    k_range = 2:3, stability_resamples = 2)
  expect_equal(sum(pam$ModelInfo$fit_table$Recommended), 1)
  expect_true(nrow(pam$Stability$cluster_recovery) > 0)
})

test_that("HDBSCAN stability reports noise recovery and handles candidate grids", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Test <- SimulatedPhenotypeData[SimulatedPhenotypeData$Cohort == "Training", ][1:160, ]
  model <- CreateHDBSCANClusterModel(df_Test, c("DensityX", "DensityY"),
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
    mclust = CreateMclustClusterModel(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      final_model = "EEI", stability_resamples = 1),
    kmeans = CreateKMeansClusterModel(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      nstart = 5, stability_resamples = 1),
    pca_mclust = CreatePCAMclustClusterModel(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      final_model = "EEI", stability_resamples = 1),
    pca_kmeans = CreatePCAKMeansClusterModel(
      df_Test, vars_Numeric, method = "finalize", final_k = 4,
      nstart = 5, stability_resamples = 1),
    hdbscan = CreateHDBSCANClusterModel(
      df_Test, c("DensityX", "DensityY"), method = "finalize",
      final_minPts = 8, final_cluster_selection_epsilon = 0,
      stability_resamples = 1),
    gower = CreateGowerPAMClusterModel(
      df_Test, c("Var1", "Var2", vars_Categorical), method = "finalize",
      final_k = 4, stability_resamples = 1),
    lca = CreateLatentClassClusterModel(
      df_Test, vars_Categorical, method = "finalize", final_k = 4,
      nrep = 2, stability_resamples = 1),
    mca_mclust = CreateMCAMclustClusterModel(
      df_Test, vars_Categorical, method = "finalize", final_k = 4,
      final_model = "EEI", stability_resamples = 1)
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
  model <- CreateSOMHDBSCANClusterModel(
    df_Training, variables = paste0("Var", 1:4), minPts_range = 3L,
    cluster_selection_epsilon_range = 0, som_xdim = 3, som_ydim = 3,
    models = 1, min_nodes_per_cluster = NULL, lpa_timeout_seconds = 30,
    Relabel = FALSE)
  projected <- ProjectSOMHDBSCANCluster(model, df_Training[1:5, ])
  expect_s3_class(model, "Pipeline_SOMHDBSCAN")
  expect_true(all(c("ModelInfo_SOM", "ModelInfo_HDBSCAN", "Specification") %in%
    names(model)))
  expect_s3_class(projected, "Project_SOMHDBSCAN")
  expect_true(all(c("individual", "summary", "by_cluster", "out_of_support",
    "policy", "plots") %in% names(projected$ProjectionFit)))
})

test_that("projections triage cases against the frozen training reference", {
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:120, ]
  df_New <- subset(SimulatedPhenotypeData, Cohort == "Projection")[1:80, ]
  model <- CreateMclustClusterModel(
    df_Test, paste0("Var", 1:6), method = "finalize", final_k = 3,
    final_model = "EEI")
  projection <- ProjectMclustCluster(model, df_New)

  expect_true(is.finite(
    model$ModelInfo$FitDiagnostics$distance_reference$high_distance_cutoff))
  expect_s3_class(projection$ProjectionFit$fit_class, "factor")
  expect_identical(
    levels(projection$ProjectionFit$fit_class),
    c("Good fit", "Uncertain membership", "Poor fit to training structure",
      "Potential novel phenotype"))
  expect_equal(
    nrow(projection$ProjectionFit$overall_fit_summary), 1)
  expect_true("Projection_Fit_Class" %in% names(projection$df_with_clusters))
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
  model <- CreateSOMMclustClusterModel(
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
    CreateKMeansClusterModel(
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
  model <- CreateHDBSCANClusterModel(
    df_Test, c("x", "y"), method = "finalize", final_minPts = 20,
    final_cluster_selection_epsilon = 0, stability_resamples = 1)
  projected <- ProjectHDBSCANCluster(model, df_Test)
  expect_true(all(model$ProbFit$individual$Cluster == 0))
  expect_true(all(projected$ProbFit$individual$Cluster == 0))
  expect_true(all(!projected$ProbFit$individual$InTrainingSupport))
})
