#' Project cases through a fitted clustering model
#'
#' @description Projects `data` through the frozen preprocessing, reduction,
#' and clustering layers stored in `object`. Dispatch is determined by the
#' fitted model's `Pipeline_*` class; callers do not select a method-specific
#' projector.
#'
#' @param object A finalized object returned by `CreateClusterModel_*()`.
#' @param new_df Deprecated alias for data, accepted by every method.
#' @param data New data to project into the frozen cluster structure.
#' @param distance_reference SOM distance reference: `loo` (default) or `training`.
#'   LOO removes each training person from their own node mean, allowing another
#'   node instead, without refits. Older models recover training scores or
#'   message a fallback. A shared quantile q_star calibrates the combined overall OR cluster LOO
#'   rule to 1-q within the attainable empirical step. The stored training
#'   rate is used by chance_flag_rate and novelty tests. Report the held-out
#'   flag rate; remaining excess comes from the map moving without the person
#'   and would need refits to remove. No further correction factor is applied.
#' @param novel_group_min_n Minimum far people in a significant merged SOM
#'   neighbourhood group, default 3. Each neighbourhood uses its training LOO
#'   far rate, supplemented to ten observations at the calibrated training rate if training n < 10.
#'   One-sided binomial tests use Holm adjustment. Overall excess is a separate
#'   one-sided binomial test of all far people against the calibrated training rate.
#' @param ... Method-specific projection options. SOM + Mclust accepts its
#'   frozen-distance and posterior-probability threshold overrides.
#' @return A method-specific projection result with the common `ProjectionFit`
#'   contract and `ProbFit` assignment table.
#'
#' @section SOM + Mclust projection results:
#' `NoveltyGroups` contains `groups` (Group, Nodes, N, Members row IDs,
#' ShareOfFar, Clusters, NearestCluster, MeanZ, TrainingClusterMeanZ, DifferenceZ,
#' Candidate, Label), `neighbourhood_tests` (NodeID, ProjectedN, TrainingN, TrainingFarN,
#' TrainingFarRate, FarN, ExpectedN, P, HolmP), `overall_excess` (N, FarN,
#' ExpectedN, P, Excess) and `verdict`. Candidate groups meet novel_group_min_n;
#' neighbourhood groups below that size remain descriptive. Group evidence and
#' overall far-from-map excess are separate from each person's Flag_FarFromMap.
#' Per-person `FitClass` has levels Good fit, Uncertain membership, Far from map,
#' clear membership, and Far from map, unclear membership. Projection_Fit_Class
#' mirrors FitClass. FitClass_Legacy is deprecated and retained for one release.
#' Soft models judge clarity using max(SoftProb_*); hard models use max_prob.
#' NodeN and Flag_OnePersonNode appear in individual outputs and DataWithClusters.
#' Specification$projection records distance_reference, far_distance_quantile,
#' cutoff_quantile (the shared q_star), chance_flag_rate (the attained rate),
#' clear_membership_threshold, novel_group_min_n and membership_probability.
#' MembershipProbability in ProjectionFit$individual contains the chosen
#' per-person clarity probability.
#' Summary metrics add one_person_node_share, rough_expected_rate_unfamiliar_profiles
#' (training share on one-person nodes; descriptive only, no inference),
#' chance_flag_rate and far_from_map_share.
#'

#' SOM + Mclust projections keep assignment certainty and training-space fit
#' separate. `ProbFit$individual` is full length and contains the assigned
#' `Cluster`, `SOM_Node`, available node-derived `prob_*` columns,
#' `prob_assigned`, `max_prob`, and `uncertainty`. It answers: *which
#' phenotype was assigned, and how certain is that node-derived assignment?*
#' Models fitted with the hard/soft projection interface also retain
#' `Cluster_Hard`, `Cluster_Soft`, `SoftProb_*`, `SoftEntropy`, `SoftMargin`,
#' `Flag_BoundaryParticipant`, and `BoundaryReason`. Node-level `prob_*`
#' values remain the Mclust posterior for a map node; `SoftProb_*` is a
#' participant-level proximity-weighted combination of nearby node posteriors.
#' Boundary flags are confidence diagnostics and never suppress projection.
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
ProjectCluster <- function(object, data, ...) {
  UseMethod("ProjectCluster")
}

#' @export
ProjectCluster.default <- function(object, data, new_df = lifecycle::deprecated(), ...) {
  if (lifecycle::is_present(new_df)) {
    lifecycle::deprecate_warn("21.18.0", "ProjectCluster(new_df)", "ProjectCluster(data)")
    data <- new_df
  }

  stop(
    "object is not a supported finalized clustering model. Expected a Pipeline_* object returned by CreateClusterModel_*().",
    call. = FALSE
  )
}
