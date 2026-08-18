#' Project cases through a fitted clustering model
#'
#' @description Projects `new_df` through the frozen preprocessing, reduction,
#' and clustering layers stored in `object`. Dispatch is determined by the
#' fitted model's `Pipeline_*` class; callers do not select a method-specific
#' projector.
#'
#' @param object A finalized object returned by `CreateClusterModel_*()`.
#' @param new_df New data to project into the frozen cluster structure.
#' @param ... Method-specific projection options. SOM + Mclust accepts its
#'   frozen-distance and posterior-probability threshold overrides.
#' @return A method-specific projection result with the common `ProjectionFit`
#'   contract and `ProbFit` assignment table.
#'
#' @section SOM + Mclust projection results:
#' SOM + Mclust projections keep assignment certainty and training-space fit
#' separate. `ProbFit$individual` is full length and contains the assigned
#' `Cluster`, `SOM_Node`, available node-derived `prob_*` columns,
#' `prob_assigned`, `max_prob`, and `uncertainty`. It answers: *which
#' phenotype was assigned, and how certain is that node-derived assignment?*
#'
#' `ProjectionFit$individual` is also full length. It contains the frozen
#' best-matching-unit `SOM_Node`, `SOM_Distance`,
#' `SOMDist_percentile_train`, `Flag_SOMDist_overallHigh`,
#' `Flag_SOMDist_clusterHigh`, `Flag_OutsideTrainingRange`, and
#' `Projection_Fit_Class`. It answers: *how well does this participant fit the
#' training phenotype space?* A high-distance flag compares the participant's
#' BMU distance with the frozen training 95th-percentile reference; it does
#' not mean the cluster assignment is necessarily wrong.
#'
#' `ProjectionFit$summary` is a metric/value table for cohort size,
#' completeness, training and projected distance summaries, high-distance
#' burden, `phenotype_drift_index`, and `node_occupancy_js_divergence`.
#' The latter is Jensen--Shannon divergence between the training and projected
#' proportions assigned to each frozen SOM node. It is 0 when the two node
#' occupancy distributions are identical and increases as their cohort-level
#' map use differs; it is neither an individual fit score nor proof that the
#' model failed.
#'
#' `ProjectionFit$by_cluster` summarizes projected sample size, BMU distance,
#' high-distance burden, and assignment probabilities by assigned cluster.
#' `ProjectionFit$out_of_support` records variables with projected values below
#' the observed training minimum or above the observed training maximum,
#' including counts and the training/projected limits. These violations produce
#' a warning but do not block complete cases from being mapped.
#'
#' `ProjectionFit$policy` records the frozen distance and probability
#' thresholds. `ProjectionFit$plots` contains projected-distance, QQ,
#' frozen-grid, and fit-class views. `ProbFit$plots$individual_ProbAssignedDensity`
#' uses a 0--1 density when assignment probabilities vary; when they are
#' essentially constant, it uses an explicit spike/rug display instead of an
#' uninformative empty density.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_New <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:6), method = "finalize", final_k = 3
#' )
#' projected <- ProjectCluster(model, df_New)
#' projected$ProjectionFit$summary
#' }
#' @export
ProjectCluster <- function(object, new_df, ...) {
  UseMethod("ProjectCluster")
}

#' @export
ProjectCluster.default <- function(object, new_df, ...) {
  stop(
    "object is not a supported finalized clustering model. Expected a Pipeline_* object returned by CreateClusterModel_*().",
    call. = FALSE
  )
}
