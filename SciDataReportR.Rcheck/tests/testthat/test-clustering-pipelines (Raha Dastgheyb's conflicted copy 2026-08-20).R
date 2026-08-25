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
  lca_model <- CreateClusterModel_LatentClass(df_Categorical, c("symptom", "diagnosis"),
    method = "finalize", final_k = 2, nrep = 2)
  expect_equal(nrow(ProjectCluster(lca_model, df_Categorical)$DataWithClusters), nrow(df_Categorical))
  mca_model <- CreateClusterModel_MCA_MClust(df_Categorical, c("symptom", "diagnosis"), method = "finalize", final_k = 2, final_model = 1)
  expect_equal(nrow(ProjectCluster(mca_model, df_Categorical)$DataWithClusters), nrow(df_Categorical))
})

test_that("categorical constructors do not infer the lifecycle from final settings", {
  df_Categorical <- data.frame(
    symptom = factor(rep(c("none", "mild"), each = 20)),
    diagnosis = factor(rep(c("control", "case"), 20))
  )
  skip_if_not_installed("poLCA")
  skip_if_not_installed("FactoMineR")
  expect_error(CreateClusterModel_LatentClass(
    df_Categorical, c("symptom", "diagnosis"), method = "exploratory",
    final_k = 2), "final_k")
  expect_error(CreateClusterModel_MCA_MClust(
    df_Categorical, c("symptom", "diagnosis"), method = "finalize",
    k_range = 2:3, final_k = 2, final_model = 1), "k_range")
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

test_that("SOM HDBSCAN rejects density settings from the wrong lifecycle", {
  skip_if_not_installed("dbscan")
  df_Test <- data.frame(x = rnorm(30), y = rnorm(30))
  expect_error(CreateClusterModel_SOM_HDBSCAN(
    df_Test, c("x", "y"), method = "finalize", minPts_range = 2:3,
    final_minPts = 3, final_cluster_selection_epsilon = 0), "minPts_range")
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
  expect_identical(model$Stability$settings$resample_type,
    "subsample_without_replacement")
  expect_equal(model$Stability$settings$resample_fraction, 0.80)
  expect_identical(model$method, "finalize")
  # The SOM pipeline keeps its own figures where they were; it gains stability
  # figures, and never a flat plots list.
  expect_null(model$plots)
  expect_false("bootstrap_agreement" %in% names(model$Stability$plots))
  expect_s3_class(model$Stability$plots$cluster_recovery, "ggplot")
  expect_s3_class(model$ModelInfo_SOM$plots$Circular, "htmlwidget")
  expect_false("occupancy" %in% names(model$ModelInfo_SOM$plots))
})

test_that("SOM projection separates assignment certainty from frozen-map fit", {
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  skip_if_not_installed("mclust")
  data("SimulatedPhenotypeData")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  df_Projection <- subset(SimulatedPhenotypeData, Cohort == "Projection")
  model <- CreateClusterModel_SOM_MClust(
    df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4,
    final_model = 3, som_xdim = 5, som_ydim = 5,
    lpa_timeout_seconds = NULL, min_nodes_per_cluster = NULL)
  warning_messages <- character()
  projected <- withCallingHandlers(
    ProjectCluster(model, df_Projection),
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("ProjectionFit\\$out_of_support", warning_messages)))
  expect_true(any(grepl("ProjectionFit\\$individual.*ProjectionFit\\$by_cluster",
    warning_messages)))
  expect_true(all(c("individual", "summary", "by_cluster", "out_of_support",
    "policy", "plots") %in% names(projected$ProjectionFit)))
  expect_false(any(c("SOMProj", "out_of_range_summary") %in% names(projected)))
  expect_false("table" %in% names(projected$ProjectionFit))
  expect_identical(nrow(projected$DataWithClusters), nrow(df_Projection))
  expect_identical(nrow(projected$ProbFit$individual), nrow(df_Projection))
  expect_identical(nrow(projected$ProjectionFit$individual), nrow(df_Projection))
  expect_false("SOM_Distance" %in% names(projected$ProbFit$individual))
  expect_true(all(c("SOM_Distance", "Flag_SOMDist_overallHigh",
    "Flag_OutsideTrainingRange", "Projection_Fit_Class") %in%
    names(projected$ProjectionFit$individual)))
  expect_true(nrow(projected$ProjectionFit$out_of_support) > 0)
  expect_true("node_occupancy_js_divergence" %in%
    projected$ProjectionFit$summary$metric)
  expect_true(all(c("training_vs_projected_distance_density",
    "training_vs_projected_distance_qq") %in%
    names(projected$ProjectionFit$plots)))
  expect_false(any(c("training_vs_projected_distance_ecdf",
    "projected_posterior_by_cluster") %in% names(projected$ProjectionFit$plots)))
  expect_s3_class(projected$ProbFit$plots$individual_ProbAssignedDensity, "ggplot")
})

test_that("node-occupancy Jensen-Shannon divergence has the expected direction", {
  expect_equal(.NodeOccupancyJSD(c(5, 5, 10), c(5, 5, 10)), 0)
  expect_gt(.NodeOccupancyJSD(c(5, 5, 10), c(0, 0, 20)), 0)
})

test_that("SOM projected posterior plot has a clear near-constant fallback", {
  fallback <- .SOMAssignedProbabilityPlot(dplyr::tibble(
    Cluster = factor(c(1, 1, 2)), prob_assigned = c(1, 1, 1)))
  expect_s3_class(fallback, "ggplot")
  expect_match(fallback$labels$subtitle, "Near-constant probabilities")
  expect_equal(fallback$coordinates$limits$x, c(0, 1))
})

test_that("SOM plus HDBSCAN uses the same unique-participant stability contract", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:80, ]
  model <- CreateClusterModel_SOM_HDBSCAN(
    df_Test, paste0("Var", 1:6), method = "finalize",
    final_minPts = 2, final_cluster_selection_epsilon = 0,
    som_xdim = 4, som_ydim = 4, min_nodes_per_cluster = NULL,
    lpa_timeout_seconds = NULL, stability_resamples = 1)
  expect_identical(model$Stability$settings$resample_type,
    "subsample_without_replacement")
  expect_equal(model$Stability$settings$resample_fraction, 0.80)
  expect_true(all(c("MinPts", "Epsilon") %in% names(model$Stability$summary)))
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
  expect_identical(failed$settings$resample_type,
    "subsample_without_replacement")
  expect_equal(failed$settings$resample_fraction, 0.80)
})

test_that("primary stability refits use unique 80 percent participant subsamples", {
  observed_ids <- integer()
  result <- .ClusterBootstrapStability(
    data.frame(id = 1:10), reference = rep(1:2, each = 5),
    fit_project = function(boot, original) {
      observed_ids <<- boot$id
      rep(1:2, each = 5)
    }, resamples = 1, seed = 11)
  expect_equal(length(observed_ids), 8)
  expect_equal(length(unique(observed_ids)), 8)
  expect_equal(result$settings$resample_fraction, 0.80)
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

test_that("constructors reject explicitly supplied settings from the other lifecycle", {
  df_Test <- data.frame(x = rnorm(40), y = rnorm(40))
  expect_error(CreateClusterModel_MClust(df_Test, c("x", "y"),
    method = "finalize", k_range = 2:3, final_k = 2, final_model = 1), "k_range")
  expect_error(CreateClusterModel_KMeans(df_Test, c("x", "y"),
    method = "exploratory", final_k = 2), "final_k")
  expect_error(CreateClusterModel_PCA_MClust(df_Test, c("x", "y"),
    method = "finalize", models = 1, final_k = 2, final_model = 1), "models")
  expect_error(CreateClusterModel_HDBSCAN(df_Test, c("x", "y"),
    method = "finalize", minPts_range = 2:3, final_minPts = 3,
    final_cluster_selection_epsilon = 0), "minPts_range")
  expect_error(CreateClusterModel_Gower_PAM(df_Test, c("x", "y"),
    method = "exploratory", final_k = 2), "final_k")
  expect_error(CreateClusterModel_HDBSCAN(df_Test, c("x", "y"),
    method = "exploratory", final_minPts = 3,
    final_cluster_selection_epsilon = 0), "final_minPts")
  expect_true(length(.ClusterSpecification("test", "x")$rng_kind) == 3L)
})

test_that("Gower PAM fits and assigns with one metric", {
  skip_if_not_installed("cluster")
  set.seed(11)
  n <- 60
  df_Mixed <- data.frame(
    Num = stats::rnorm(n),
    Ord = factor(sample(c("none", "mild", "severe"), n, TRUE),
      levels = c("none", "mild", "severe"), ordered = TRUE),
    Nom = factor(sample(c("a", "b", "c"), n, TRUE)),
    Flag = sample(c(TRUE, FALSE), n, TRUE))
  variables <- names(df_Mixed)

  # The frozen specification has to reproduce the coefficient PAM was fitted
  # on, including daisy's ordinal and asymmetric-binary handling.
  specification <- .GowerSpecification(df_Mixed, variables)
  expect_equal(
    .GowerToMedoids(df_Mixed, df_Mixed, variables, specification),
    as.matrix(suppressWarnings(cluster::daisy(df_Mixed, metric = "gower"))),
    ignore_attr = TRUE)

  # Projecting the training cohort through its own model must return the
  # training labels; a second, different metric used to move about 8% of them.
  model <- suppressWarnings(CreateClusterModel_Gower_PAM(
    df_Mixed, variables, method = "finalize", final_k = 3))
  projected <- suppressWarnings(ProjectCluster(model, df_Mixed))
  expect_identical(
    projected$ProbFit$individual$Cluster, model$ProbFit$individual$Cluster)
})

test_that("confidence figures only pin a probability to a 0-1 axis", {
  data("SimulatedPhenotypeData")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  variables <- paste0("Var", 1:6)

  # An assignment margin is an unbounded distance, so the density has to show
  # the values above 1 rather than silently dropping them.
  km <- CreateClusterModel_KMeans(
    df_Training, variables, method = "finalize", final_k = 3)
  expect_gt(max(km$ProbFit$individual$AssignmentMargin, na.rm = TRUE), 1)
  built <- ggplot2::ggplot_build(km$ProbFit$plots$confidence_density)
  expect_gt(max(built$data[[1]]$x), 1)

  # A posterior probability still gets the fixed 0-1 axis.
  mc <- CreateClusterModel_MClust(
    df_Training, variables, method = "finalize", final_k = 3, final_model = 1)
  limits <- lapply(mc$ProbFit$plots$confidence_density$scales$scales,
    function(scale) scale$limits)
  expect_true(any(vapply(limits, identical, logical(1), c(0, 1))))
})

test_that("HDBSCAN triages projections against a matching training reference", {
  skip_if_not_installed("dbscan")
  data("SimulatedPhenotypeData")
  df_Training <- subset(SimulatedPhenotypeData, Cohort == "Training")
  model <- CreateClusterModel_HDBSCAN(
    df_Training, c("DensityX", "DensityY"), method = "finalize",
    final_minPts = 10, final_cluster_selection_epsilon = 0)

  # The frozen cutoff has to come from the same quantity a projected case is
  # scored on, not from the training core distances.
  expect_true("NearestCoreDistance" %in% names(model$ProbFit$individual))
  reference <- model$ModelInfo$FitDiagnostics$distance_reference
  expect_equal(
    reference$high_distance_cutoff,
    unname(stats::quantile(
      model$ProbFit$individual$NearestCoreDistance, 0.95, na.rm = TRUE)))

  projected <- ProjectCluster(model, df_Training)
  expect_equal(
    mean(projected$ProbFit$individual$Projection_Fit_Class ==
      "Good fit", na.rm = TRUE),
    0.95, tolerance = 0.05)
})

test_that("SOM plus HDBSCAN clusters the codebook without fitting the LPA grid", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("kohonen")
  skip_if_not_installed("aweSOM")
  suppressPackageStartupMessages(skip_if_not_installed("tidyLPA"))
  data("SimulatedPhenotypeData")
  df_Test <- subset(SimulatedPhenotypeData, Cohort == "Training")[1:80, ]
  model <- suppressWarnings(CreateClusterModel_SOM_HDBSCAN(
    df_Test, paste0("Var", 1:6), method = "finalize", final_minPts = 2,
    final_cluster_selection_epsilon = 0, som_xdim = 4, som_ydim = 4,
    min_nodes_per_cluster = NULL, Relabel = FALSE))
  expect_null(model$ModelInfo_MClust$lpa_models)
  expect_equal(nrow(model$ModelInfo_MClust$diagnostics$lpa_fit_diagnostics), 0)

  # The SOM training diagnostics have to describe the clusters the model
  # reports, not a discarded latent-profile solution.
  expect_identical(
    model$ModelInfo_SOM$SOMFit$table$Cluster,
    model$ProbFit$individual$Cluster)
  expect_true(inherits(model$fit_plot, "ggplot"))
})
