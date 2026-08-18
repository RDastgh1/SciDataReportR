#' Fit HDBSCAN clusters on a frozen self-organizing map
#'
#' @description Trains the same frozen SOM used by [CreateSOMMclustClusterModel()],
#' then applies HDBSCAN to its node codebook. Participants inherit the cluster
#' (or noise) label of their best-matching node. This is a node-based phenotype
#' model: projected participants are mapped to the original nodes and are not
#' refit with HDBSCAN.
#' @param data Data frame used to train the SOM and node-level HDBSCAN model.
#' @param variables Numeric variables used for SOM training.
#' @param minPts_range Candidate HDBSCAN minimum-point settings for SOM nodes.
#' @param cluster_selection_epsilon_range Candidate HDBSCAN extraction epsilon settings.
#' @param final_minPts Optional finalized HDBSCAN minimum-point setting.
#' @param final_cluster_selection_epsilon Optional finalized extraction epsilon.
#' @param ClusterName Name of the appended cluster column.
#' @param seed_som Seed used for SOM training.
#' @param seed_hdbscan Seed retained in the model specification.
#' @param ... Additional arguments passed to [CreateSOMMclustClusterModel()].
#' @return A `Pipeline_SOMHDBSCAN` object containing frozen SOM and HDBSCAN
#'   node models, assignments, diagnostics, and a projection specification.
#' @export
CreateSOMHDBSCANClusterModel <- function(data, variables = NULL,
    minPts_range = 2:10,
    cluster_selection_epsilon_range = c(0, 0.05, 0.10),
    final_minPts = NULL, final_cluster_selection_epsilon = NULL,
    ClusterName = "Cluster", seed_som = 934521L, seed_hdbscan = 93421L,
    ...) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required.")
  }
  # Reuse the established, tested SOM fitting and frozen scaling workflow.
  som <- CreateSOMMclustClusterModel(
    data = data, variables = variables, method = "exploratory", ...,
    ClusterName = ClusterName, seed_som = seed_som, seed_lpa = seed_hdbscan)
  codes <- som$ModelInfo_SOM$som_codes
  grid <- if (!is.null(final_minPts)) {
    expand.grid(MinPts = final_minPts,
      Epsilon = if (is.null(final_cluster_selection_epsilon)) 0 else final_cluster_selection_epsilon)
  } else expand.grid(MinPts = minPts_range,
    Epsilon = cluster_selection_epsilon_range)
  fits <- lapply(seq_len(nrow(grid)), function(i) {
    dbscan::hdbscan(codes, minPts = grid$MinPts[[i]],
      cluster_selection_epsilon = grid$Epsilon[[i]])
  })
  fit_table <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    sizes <- table(fit$cluster[fit$cluster > 0])
    dplyr::tibble(MinPts = grid$MinPts[[i]], Epsilon = grid$Epsilon[[i]],
      Classes = length(sizes), Persistence = mean(fit$cluster_scores, na.rm = TRUE),
      NoiseProportion = mean(fit$cluster == 0),
      MinClusterN = if (length(sizes)) min(sizes) else 0L)
  }))
  review <- .ClusterAHP(fit_table,
    higher = c("Persistence", "MinClusterN"), lower = "NoiseProportion",
    setting = "SOM-node HDBSCAN setting")
  best_i <- which(review$fit_table$Recommended)[1]
  node_fit <- fits[[best_i]]
  node_cluster <- node_fit$cluster
  individual <- som$ProbFit$individual
  individual$Cluster <- node_cluster[individual$SOM_Node]
  individual$Cluster[is.na(individual$SOM_Node)] <- NA_integer_
  som$df_with_clusters[[ClusterName]] <- individual$Cluster
  som$ProbFit$individual <- individual
  ModelInfo_HDBSCAN <- list(hdbscan_model = node_fit, node_cluster = node_cluster,
    fit_table = review$fit_table, AHP = review$AHP,
    final_minPts = grid$MinPts[[best_i]],
    final_cluster_selection_epsilon = grid$Epsilon[[best_i]])
  som$ModelInfo_HDBSCAN <- ModelInfo_HDBSCAN
  som$ModelInfo <- ModelInfo_HDBSCAN
  som$fit_plot <- PlotClusterFitReview(review$fit_table, x = "MinPts")
  som$CandidateAudit <- dplyr::mutate(review$fit_table, Status = "success",
    RuntimeSeconds = NA_real_, Warnings = NA_character_, Error = NA_character_)
  som$Specification <- .ClusterSpecification("SOM + HDBSCAN", som$vars_used,
    list(seed_som = seed_som, seed_hdbscan = seed_hdbscan),
    list(ZScoreObject = som$ZScoreObject, SOM = som$ModelInfo_SOM$som_grid_info),
    grid, list(minPts = grid$MinPts[[best_i]], Epsilon = grid$Epsilon[[best_i]]),
    som$complete_rows, list(projection = "frozen SOM node to HDBSCAN label"))
  class(som) <- c("Pipeline_SOMHDBSCAN", class(som))
  som
}

#' Project cases onto a SOM + HDBSCAN model
#'
#' @description Maps cases to the frozen SOM, then assigns the HDBSCAN label
#' of each best-matching node. This is not a refitted or native HDBSCAN
#' prediction; its distance diagnostic is distance to the frozen SOM BMU.
#' @param object A `CreateSOMHDBSCANClusterModel()` result.
#' @param new_df New cases to map to the frozen SOM nodes.
#' @param ClusterName Name of the projected cluster column.
#' @param ... Additional arguments passed to [ProjectSOMMclustCluster()].
#' @return A `Project_SOMHDBSCAN` result. Assignment is inherited from the
#'   projected case's frozen best-matching SOM node; it is not a native
#'   HDBSCAN prediction.
#' @export
ProjectSOMHDBSCANCluster <- function(object, new_df,
    ClusterName = object$ClusterName, ...) {
  if (!inherits(object, "Pipeline_SOMHDBSCAN")) {
    stop("object must be a Pipeline_SOMHDBSCAN result.")
  }
  # The canonical SOM + Mclust projector accepts the inherited class.
  projected <- ProjectSOMMclustCluster(object, new_df, ClusterName = ClusterName, ...)
  node_cluster <- object$ModelInfo_HDBSCAN$node_cluster
  projected$ProbFit$individual$Cluster <- node_cluster[projected$ProbFit$individual$SOM_Node]
  projected$df_with_clusters[[ClusterName]] <- projected$ProbFit$individual$Cluster
  projected$ProbFit$individual$Projection_Fit_Class <- dplyr::case_when(
    is.na(projected$ProbFit$individual$SOM_Distance) ~ NA_character_,
    projected$ProbFit$individual$Cluster == 0 ~ "Potential novel phenotype",
    TRUE ~ "Good fit")
  projected$df_with_clusters$Projection_Fit_Class <-
    projected$ProbFit$individual$Projection_Fit_Class
  projected$ProjectionFit$individual <- projected$ProbFit$individual
  projected$ProjectionFit$policy <- list(
    assignment = "frozen SOM best-matching node HDBSCAN label",
    distance_metric = "distance to frozen SOM best-matching unit")
  projected$ModelInfo_HDBSCAN <- object$ModelInfo_HDBSCAN
  projected$Specification <- object$Specification
  class(projected) <- c("Project_SOMHDBSCAN", class(projected))
  projected
}
