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
#' @param distance_reference SOM distance reference for the novelty
#'   p-values: `loo` (default) or `training`. LOO removes each training person
#'   from their own node mean, allowing another node instead, without refits,
#'   so training people are scored as a new person would be. `training` uses
#'   in-sample distances, which are optimistic. Older models recover training
#'   scores or message a fallback to `training`.
#' @param novel_group_min_n Minimum unusual people in a group for it to be a
#'   candidate new phenotype (default 3).
#' @param ... Method-specific projection options. SOM + Mclust accepts its
#'   frozen-distance and posterior-probability threshold overrides.
#' @return A method-specific projection result with the common `ProjectionFit`
#'   contract and `ProbFit` assignment table.
#'
#' @section SOM + Mclust novelty: three questions, answered separately:
#' \describe{
#'   \item{Is this person outside the training data?}{Each projected person's
#'     distance to their best-matching node is compared with the training
#'     people's distances (leave-one-out by default). `NoveltyP` is the
#'     conformal p-value (1 + training people at least as far) / (n + 1),
#'     which is valid without distributional assumptions. `NoveltyQ` is its
#'     Benjamini-Hochberg adjustment across the projected cohort, and
#'     `Flag_Outlier` is `NoveltyQ <= novelty_fdr`. A cohort that resembles
#'     training therefore gives few or no outliers. `Top_Distance_Variable`
#'     and `Top_Distance_Share` name the variable that contributes most to
#'     the distance.}
#'   \item{Does the cohort as a whole differ from training?}{
#'     `ProjectionFit$cohort_shift` has one row per test (Test, Question,
#'     Statistic, StatisticName, Expected, P): a one-sided rank test that
#'     projected distances are larger (statistic = probability that a
#'     projected person is farther than a training person; 0.5 = no shift), a
#'     Kolmogorov-Smirnov test of the distance distribution, a binomial test
#'     of the share with `NoveltyP <= group_screen_alpha`, and Monte Carlo
#'     chi-square tests of node occupancy and cluster prevalence against
#'     training proportions.}
#'   \item{Do unusual people form a new group?}{People with
#'     `NoveltyP <= group_screen_alpha` are `Flag_Unusual`. For each map node
#'     and its neighbours, a one-sided binomial test compares the number of
#'     unusual projected people with that neighbourhood's own training rate
#'     (shrunk toward chance with ten pseudo-people, so no neighbourhood has
#'     an expected rate of zero);
#'     Holm-significant neighbourhoods are merged into groups. A group's
#'     `Spread` (mean pairwise distance between members) is compared with
#'     same-size random sets of training's unusual people (`CoherenceP`). A
#'     group is a `Candidate` new phenotype when it has at least
#'     `novel_group_min_n` members and `CoherenceP < 0.05`.}
#' }
#' `NoveltyGroups` contains `groups` (Group, Nodes, N, Members row IDs,
#' ShareOfUnusual, Clusters, NearestCluster, Spread, ReferenceSpread,
#' CoherenceP, MeanZ, TrainingClusterMeanZ, DifferenceZ, Candidate, Label),
#' `neighbourhood_tests`, `overall_excess` and a one-line `verdict`.
#'
#' Per-person `FitClass` has levels Good fit, Uncertain membership (membership
#' probability below `clear_membership_threshold`: between known clusters),
#' Outlier, and Member of candidate new group, in increasing precedence.
#' Membership probability is not used to judge novelty: a Gaussian mixture
#' gives a confident label even to someone far from every cluster. Soft models
#' judge clarity using max(SoftProb_*); hard models use max_prob. Hard-label
#' SOM pipelines (K-means, hierarchical, HDBSCAN) have no membership
#' probability and use Good fit, Outlier and Member of candidate new group.
#'
#' `ProbFit$individual` is full length and contains the assigned `Cluster`,
#' `SOM_Node`, node-derived `prob_*`, `prob_assigned`, `max_prob`,
#' `uncertainty`, the hard/soft assignment columns (`Cluster_Hard`,
#' `Cluster_Soft`, `SoftProb_*`, `SoftEntropy`, `SoftMargin`,
#' `Flag_BoundaryParticipant`, `BoundaryReason`) and the novelty columns.
#' `ProjectionFit$individual` adds `SOM_Distance`, `SOMDist_percentile_train`,
#' `Flag_OutsideTrainingRange` and `Flag_Unusual`.
#'
#' `ProjectionFit$summary` is a metric/value table for cohort size,
#' completeness, training and projected distance summaries,
#' `phenotype_drift_index`, `outlier_n`, `outlier_share`, `unusual_share`,
#' `expected_unusual_share`, `cluster_occupancy_js_divergence` and
#' `node_occupancy_js_divergence` (Jensen--Shannon divergences; 0 = identical), `one_person_node_share` and
#' `rough_expected_rate_unfamiliar_profiles` (training share on one-person
#' nodes; descriptive only). `ProjectionFit$by_cluster` summarizes distance,
#' probability, outlier, unusual, new-group and uncertain shares by cluster.
#' `ProjectionFit$out_of_support` records variables with projected values
#' beyond the observed training range; these produce a warning but complete
#' cases are still mapped.
#'
#' `ProjectionFit$policy` records the novelty and probability settings.
#' `ProjectionFit$plots` contains projected-distance, QQ, novelty p-value
#' histogram (flat when the cohort resembles training), unusual share by
#' cluster, frozen-grid and fit-class views. `ProbFit$plots$individual_ProbAssignedDensity`
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
