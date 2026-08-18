# Extracted from test-clustering-pipelines.R:242

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
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
