## -----------------------------------------------------------------------------
#| label: LoadPackages
library(SciDataReportR)
library(dplyr)
library(ggplot2)


## -----------------------------------------------------------------------------
#| label: LoadData
data("SimulatedPhenotypeData", package = "SciDataReportR")
data("SimulatedPhenotypeVariableTypes", package = "SciDataReportR")

df_Raw <- SimulatedPhenotypeData
VariableTypes <- SimulatedPhenotypeVariableTypes
df_Revalued <- RevalueData(df_Raw, VariableTypes)$RevaluedData

vars_Continuous <- paste0("Var", 1:12)
vars_Categorical <- paste0("CatVar", 1:3)
vars_Mixed <- c(vars_Continuous, vars_Categorical)
vars_Density <- c("DensityX", "DensityY")

df_Training <- df_Revalued %>%
  dplyr::filter(.data$Cohort == "Training")
df_Projected <- df_Revalued %>%
  dplyr::filter(.data$Cohort == "Projection")


## -----------------------------------------------------------------------------
#| label: BenchmarkHelpers
PrintPlots <- function(plots, plot_names = NULL) {
  if (is.null(plot_names)) plot_names <- names(plots)
  for (plot_name in intersect(plot_names, names(plots))) {
    if (inherits(plots[[plot_name]], "ggplot")) print(plots[[plot_name]])
  }
  invisible(NULL)
}

PrintFitMetrics <- function(object, metrics) {
  object$ModelInfo$fit_table %>%
    dplyr::select(dplyr::any_of(c(
      "Model", "Classes", "MinPts", "Epsilon", metrics,
      "ReproducibilityScore", "ahp_index", "Recommended")))
}

CalculateTruthARI <- function(truth, assignment) {
  valid <- !is.na(truth) & !is.na(assignment)
  if (sum(valid) < 2L) return(NA_real_)
  mclust::adjustedRandIndex(truth[valid], assignment[valid])
}

CalculateNoisePerformance <- function(truth, assignment) {
  valid <- !is.na(truth) & !is.na(assignment)
  truth_noise <- as.character(truth[valid]) == "Noise"
  assigned_noise <- assignment[valid] == 0
  dplyr::tibble(
    NoiseSensitivity = mean(assigned_noise[truth_noise]),
    NoiseSpecificity = mean(!assigned_noise[!truth_noise]))
}


## -----------------------------------------------------------------------------
#| label: MclustReview
df_MclustReview <- CreateClusterModel_MClust(
  df_Training, vars_Continuous, k_range = 2:10, model_names = "EEI",
  stability_resamples = 2
)
PrintFitMetrics(
  df_MclustReview,
  c("BIC", "ICL", "Entropy", "MinClusterN", "SizeBalance"))
df_MclustReview$ModelInfo$AHP$recommendation
print(df_MclustReview$fit_plot)
PrintPlots(df_MclustReview$ModelInfo$plots, c("bic", "icl", "entropy"))
PrintPlots(df_MclustReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: MclustFinalize
df_Mclust <- CreateClusterModel_MClust(
  df_Training, vars_Continuous, method = "finalize", final_k = 4,
  final_model = "EEI", ClusterName = "MclustCluster",
  stability_resamples = 2
)
df_MclustProjected <- ProjectCluster(df_Mclust, df_Projected)
PrintPlots(
  df_Mclust$ModelInfo$plots,
  c("separation_map", "classification_uncertainty", "centre_heatmap",
    "profiles"))
PrintPlots(df_Mclust$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_Mclust$ProbFit$plots)
PrintPlots(df_MclustProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: KMeansReview
df_KMeansReview <- CreateClusterModel_KMeans(
  df_Training, vars_Continuous, k_range = 2:10, nstart = 20,
  stability_resamples = 2
)
PrintFitMetrics(
  df_KMeansReview,
  c("WSS", "Silhouette", "CalinskiHarabasz", "MinClusterN"))
df_KMeansReview$ModelInfo$AHP$recommendation
print(df_KMeansReview$fit_plot)
PrintPlots(
  df_KMeansReview$ModelInfo$plots,
  c("elbow", "silhouette_by_k", "calinski_harabasz"))
PrintPlots(df_KMeansReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: KMeansFinalize
df_KMeans <- CreateClusterModel_KMeans(
  df_Training, vars_Continuous, method = "finalize", final_k = 4,
  nstart = 20, ClusterName = "KMeansCluster", stability_resamples = 2
)
df_KMeansProjected <- ProjectCluster(df_KMeans, df_Projected)
PrintPlots(
  df_KMeans$ModelInfo$plots,
  c("silhouette", "separation_map", "centre_heatmap", "profiles"))
PrintPlots(df_KMeans$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_KMeans$ProbFit$plots)
PrintPlots(df_KMeansProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: PCAMclustReview
df_PCAMclustReview <- CreateClusterModel_PCA_MClust(
  df_Training, vars_Continuous, k_range = 2:10, model_names = "EEI",
  stability_resamples = 2
)
PrintFitMetrics(
  df_PCAMclustReview,
  c("BIC", "ICL", "Entropy", "MinClusterN"))
df_PCAMclustReview$ModelInfo$AHP$recommendation
ncol(df_PCAMclustReview$PCA$Scores)
print(df_PCAMclustReview$fit_plot)
PrintPlots(df_PCAMclustReview$ModelInfo$plots, c("scree", "bic", "icl"))
PrintPlots(df_PCAMclustReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: PCAMclustFinalize
df_PCAMclust <- CreateClusterModel_PCA_MClust(
  df_Training, vars_Continuous, method = "finalize", final_k = 4,
  final_model = "EEI", ClusterName = "PCAMclustCluster",
  stability_resamples = 2
)
df_PCAMclustProjected <- ProjectCluster(df_PCAMclust, df_Projected)
PrintPlots(
  df_PCAMclust$ModelInfo$plots,
  c("scree", "separation_map", "profiles"))
PrintPlots(df_PCAMclust$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_PCAMclust$ProbFit$plots)
PrintPlots(df_PCAMclustProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: PCAKMeansReview
df_PCAKMeansReview <- CreateClusterModel_PCA_KMeans(
  df_Training, vars_Continuous, k_range = 2:10, nstart = 20,
  stability_resamples = 2
)
PrintFitMetrics(
  df_PCAKMeansReview,
  c("WSS", "Silhouette", "CalinskiHarabasz", "MinClusterN"))
df_PCAKMeansReview$ModelInfo$AHP$recommendation
ncol(df_PCAKMeansReview$PCA$Scores)
print(df_PCAKMeansReview$fit_plot)
PrintPlots(
  df_PCAKMeansReview$ModelInfo$plots,
  c("scree", "elbow", "silhouette_by_k"))
PrintPlots(df_PCAKMeansReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: PCAKMeansFinalize
df_PCAKMeans <- CreateClusterModel_PCA_KMeans(
  df_Training, vars_Continuous, method = "finalize", final_k = 4,
  nstart = 20, ClusterName = "PCAKMeansCluster",
  stability_resamples = 2
)
df_PCAKMeansProjected <- ProjectCluster(df_PCAKMeans, df_Projected)
PrintPlots(
  df_PCAKMeans$ModelInfo$plots,
  c("scree", "silhouette", "separation_map", "profiles"))
PrintPlots(df_PCAKMeans$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_PCAKMeans$ProbFit$plots)
PrintPlots(df_PCAKMeansProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: HDBSCANReview
df_HDBSCANReview <- CreateClusterModel_HDBSCAN(
  df_Training, vars_Density, minPts_range = c(6, 10, 14),
  cluster_selection_epsilon_range = c(0, 0.05), stability_resamples = 2
)
PrintFitMetrics(
  df_HDBSCANReview,
  c("Clusters", "Persistence", "MeanMembershipProbability",
    "NoiseProportion", "MinClusterN"))
df_HDBSCANReview$ModelInfo$AHP$recommendation
print(df_HDBSCANReview$fit_plot)
PrintPlots(df_HDBSCANReview$ModelInfo$plots, c("density_review", "persistence"))
PrintPlots(df_HDBSCANReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: HDBSCANFinalize
df_HDBSCAN <- CreateClusterModel_HDBSCAN(
  df_Training, vars_Density, method = "finalize", final_minPts = 10,
  final_cluster_selection_epsilon = 0, ClusterName = "HDBSCANCluster",
  stability_resamples = 2
)
df_HDBSCANProjected <- ProjectCluster(df_HDBSCAN, df_Projected)
PrintPlots(
  df_HDBSCAN$ModelInfo$plots,
  c("separation_map", "cluster_persistence", "profiles"))
PrintPlots(df_HDBSCAN$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_HDBSCAN$ProbFit$plots)
PrintPlots(df_HDBSCANProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: SOMReview
df_SOMReview <- CreateClusterModel_SOM_MClust(
  df_Training, vars_Continuous, method = "exploratory", k_range = 2:10,
  models = 1, som_xdim = 5, som_ydim = 5, stability_resamples = 2,
  Relabel = FALSE, min_nodes_per_cluster = NULL
)
df_SOMReview$ModelInfo_MClust$fit_table
df_SOMReview$ModelInfo_MClust$AHP$recommendation
df_SOMReview$Stability$summary
print(df_SOMReview$fit_plot)
df_SOMReview$ModelInfo_SOM$plots$Circular
df_SOMReview$ModelInfo_SOM$plots$Line
df_SOMReview$ModelInfo_SOM$plots$Cloud


## -----------------------------------------------------------------------------
#| label: SOMFinalize
df_SOM <- CreateClusterModel_SOM_MClust(
  df_Training, vars_Continuous, method = "finalize", final_k = 4,
  final_model = 1, som_xdim = 5, som_ydim = 5,
  ClusterName = "SOMCluster", stability_resamples = 2,
  Relabel = FALSE, min_nodes_per_cluster = NULL
)
df_SOMProjected <- ProjectCluster(df_SOM, df_Projected)
PrintPlots(df_SOM$ModelInfo_SOM$plots)
PrintPlots(df_SOM$ModelInfo_SOM$SOMFit$plots)
PrintPlots(df_SOM$ModelInfo_SOM$PhenotypeReference$plots)
PrintPlots(df_SOM$ProbFit$plots)
PrintPlots(df_SOM$Stability$plots)
PrintPlots(
  df_SOMProjected$ProjectionFit$plots,
  c("distance_by_cluster_box", "training_vs_projected_cluster_occupancy",
    "projection_fit_class_bar", "som_grid_map"))


## -----------------------------------------------------------------------------
#| label: GowerPAMReview
df_GowerReview <- CreateClusterModel_Gower_PAM(
  df_Training, vars_Mixed, k_range = 2:10, stability_resamples = 2
)
PrintFitMetrics(
  df_GowerReview,
  c("Silhouette", "MeanMedoidDistance", "MinClusterN", "SizeBalance"))
df_GowerReview$ModelInfo$AHP$recommendation
print(df_GowerReview$fit_plot)
PrintPlots(
  df_GowerReview$ModelInfo$plots,
  c("silhouette_by_k", "medoid_distance_by_k"))
PrintPlots(df_GowerReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: GowerPAMFinalize
df_Gower <- CreateClusterModel_Gower_PAM(
  df_Training, vars_Mixed, method = "finalize", final_k = 4,
  ClusterName = "GowerCluster", stability_resamples = 2
)
df_GowerProjected <- ProjectCluster(df_Gower, df_Projected)
PrintPlots(
  df_Gower$ModelInfo$plots,
  c("silhouette", "gower_map", "profiles",
    "categorical_composition"))
PrintPlots(df_Gower$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_Gower$ProbFit$plots)
PrintPlots(df_GowerProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: LatentClassReview
df_LCAReview <- CreateClusterModel_LatentClass(
  df_Training, vars_Categorical, k_range = 2:10, nrep = 5,
  stability_resamples = 2
)
PrintFitMetrics(
  df_LCAReview,
  c("BIC", "AIC", "Entropy", "MinClassN", "SizeBalance"))
df_LCAReview$ModelInfo$AHP$recommendation
print(df_LCAReview$fit_plot)
PrintPlots(df_LCAReview$ModelInfo$plots, c("bic", "entropy"))
PrintPlots(df_LCAReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: LatentClassFinalize
df_LCA <- CreateClusterModel_LatentClass(
  df_Training, vars_Categorical, method = "finalize", final_k = 4,
  nrep = 5, ClusterName = "LCACluster", stability_resamples = 2
)
df_LCAProjected <- ProjectCluster(df_LCA, df_Projected)
PrintPlots(
  df_LCA$ModelInfo$plots,
  c("response_probabilities", "item_profiles", "posterior_map",
    "categorical_composition"))
PrintPlots(df_LCA$ProbFit$plots)
PrintPlots(df_LCAProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: MCAMclustReview
df_MCAMclustReview <- CreateClusterModel_MCA_MClust(
  df_Training, vars_Categorical, k_range = 2:10, model_names = "EEI",
  stability_resamples = 2
)
PrintFitMetrics(
  df_MCAMclustReview,
  c("BIC", "ICL", "Entropy", "MinClusterN", "SizeBalance"))
df_MCAMclustReview$ModelInfo$AHP$recommendation
print(df_MCAMclustReview$fit_plot)
PrintPlots(df_MCAMclustReview$ModelInfo$plots, c("scree", "bic", "icl"))
PrintPlots(df_MCAMclustReview$Stability$plots)


## -----------------------------------------------------------------------------
#| label: MCAMclustFinalize
df_MCAMclust <- CreateClusterModel_MCA_MClust(
  df_Training, vars_Categorical, method = "finalize", final_k = 4,
  final_model = "EEI", ClusterName = "MCAMclustCluster",
  stability_resamples = 2
)
df_MCAMclustProjected <- ProjectCluster(df_MCAMclust, df_Projected)
PrintPlots(
  df_MCAMclust$ModelInfo$plots,
  c("scree", "separation_map", "categorical_composition"))
PrintPlots(df_MCAMclust$ModelInfo$FitDiagnostics$plots)
PrintPlots(df_MCAMclust$ProbFit$plots)
PrintPlots(df_MCAMclustProjected$ProjectionFit$plots)


## -----------------------------------------------------------------------------
#| label: TruthAgreement
df_TruthAgreement <- dplyr::tibble(
  Method = c(
    "Mclust", "K-means", "PCA + Mclust", "PCA + K-means",
    "Gower + PAM", "LCA", "MCA + Mclust", "SOM + Mclust"),
  TrainingARI = c(
    CalculateTruthARI(df_Training$TruthCluster, df_Mclust$df_with_clusters$MclustCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_KMeans$df_with_clusters$KMeansCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_PCAMclust$df_with_clusters$PCAMclustCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_PCAKMeans$df_with_clusters$PCAKMeansCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_Gower$df_with_clusters$GowerCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_LCA$df_with_clusters$LCACluster),
    CalculateTruthARI(df_Training$TruthCluster, df_MCAMclust$df_with_clusters$MCAMclustCluster),
    CalculateTruthARI(df_Training$TruthCluster, df_SOM$df_with_clusters$SOMCluster)),
  ProjectionARI = c(
    CalculateTruthARI(df_Projected$TruthCluster, df_MclustProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_KMeansProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_PCAMclustProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_PCAKMeansProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_GowerProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_LCAProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_MCAMclustProjected$ProbFit$individual$Cluster),
    CalculateTruthARI(df_Projected$TruthCluster, df_SOMProjected$ProbFit$individual$Cluster))
)
df_TruthAgreement


## -----------------------------------------------------------------------------
#| label: DensityTruthAgreement
df_DensityTruth <- dplyr::tibble(
  Cohort = c("Training", "Projection"),
  ARI = c(
    CalculateTruthARI(
      df_Training$TruthDensityGroup,
      df_HDBSCAN$df_with_clusters$HDBSCANCluster),
    CalculateTruthARI(
      df_Projected$TruthDensityGroup,
      df_HDBSCANProjected$ProbFit$individual$Cluster)),
  NoiseProportion = c(
    mean(df_HDBSCAN$df_with_clusters$HDBSCANCluster == 0, na.rm = TRUE),
    mean(df_HDBSCANProjected$ProbFit$individual$Cluster == 0, na.rm = TRUE))) %>%
  dplyr::bind_cols(dplyr::bind_rows(
    CalculateNoisePerformance(
      df_Training$TruthDensityGroup,
      df_HDBSCAN$df_with_clusters$HDBSCANCluster),
    CalculateNoisePerformance(
      df_Projected$TruthDensityGroup,
      df_HDBSCANProjected$ProbFit$individual$Cluster)))
df_DensityTruth


## -----------------------------------------------------------------------------
#| label: CompareClusterOccupancy
df_Comparison <- df_Mclust$df_with_clusters %>%
  dplyr::select(.data$.scidr_rowid, .data$MclustCluster) %>%
  dplyr::left_join(
    df_KMeans$df_with_clusters %>%
      dplyr::select(.data$.scidr_rowid, .data$KMeansCluster),
    by = ".scidr_rowid")

df_Comparison %>%
  dplyr::count(.data$MclustCluster, .data$KMeansCluster)


## -----------------------------------------------------------------------------
#| label: PlotClusterProfiles
PlotClusterBoxplot(
  data = df_Mclust$df_with_clusters,
  ClusterVar = "MclustCluster",
  variables = vars_Continuous,
  codebook = VariableTypes,
  Scale = TRUE
)


## -----------------------------------------------------------------------------
#| label: Reproducibility
# save.image()
print(sessionInfo())

