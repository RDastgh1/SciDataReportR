#' Fit HDBSCAN clusters on a frozen self-organizing map
#'
#' @description Trains the same frozen SOM used by
#'   [CreateClusterModel_SOM_MClust()],
#' then applies HDBSCAN to its node codebook. Participants inherit the cluster
#' (or noise) label of their best-matching node. This is a node-based phenotype
#' model: projected participants are mapped to the original nodes and are not
#' refit with HDBSCAN.
#' @param data Data frame used to train the SOM and node-level HDBSCAN model.
#' @param variables Numeric variables used for SOM training.
#' @param method Either `"exploratory"` or `"finalize"`.
#' @param minPts_range Candidate HDBSCAN minimum-point settings for SOM nodes.
#' @param cluster_selection_epsilon_range Candidate HDBSCAN extraction epsilon settings.
#' @param final_minPts Optional finalized HDBSCAN minimum-point setting.
#' @param final_cluster_selection_epsilon Optional finalized extraction epsilon.
#' @param ClusterVariableName Name of the appended cluster column.
#' @param seed_som Seed used for SOM training.
#' @param seed_hdbscan Seed retained in the model specification.
#' @param stability_resamples Number of 80% participant subsample refits used
#'   for final-model stability. Subsamples are drawn without replacement.
#' @param stability_seed Seed controlling participant subsampling.
#' @param stability_progress Whether to print subsample progress messages.
#' @param ... Additional arguments passed to [CreateClusterModel_SOM_MClust()].
#' @inheritSection cluster-stability-output Stability output
#' @return A `Pipeline_SOM_HDBSCAN` object containing frozen SOM and HDBSCAN
#'   node models, assignments, diagnostics, and a projection specification.
#'   Persistence and minimum node-cluster size are higher-is-better, noise
#'   proportion is lower-is-better, and extracted class count is data-derived.
#' @export
CreateClusterModel_SOM_HDBSCAN <- function(data, variables = NULL,
    method = c("exploratory", "finalize"),
    minPts_range = 2:10,
    cluster_selection_epsilon_range = c(0, 0.05, 0.10),
    final_minPts = NULL, final_cluster_selection_epsilon = NULL,
    ClusterVariableName = "Cluster", seed_som = 934521L, seed_hdbscan = 93421L,
    stability_resamples = 0L, stability_seed = seed_hdbscan + 1L,
    stability_progress = FALSE,
    ...) {
  supplied <- list(minPts_range = !missing(minPts_range),
    cluster_selection_epsilon_range = !missing(cluster_selection_epsilon_range),
    final_minPts = !missing(final_minPts),
    final_cluster_selection_epsilon = !missing(final_cluster_selection_epsilon))
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required.")
  }
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  method <- .ClusterMethod(method)
  .ValidateClusterLifecycle(method,
    supplied[c("minPts_range", "cluster_selection_epsilon_range")],
    supplied[c("final_minPts", "final_cluster_selection_epsilon")])
  if (method == "finalize" &&
      (is.null(final_minPts) || is.null(final_cluster_selection_epsilon))) {
    stop("final_minPts and final_cluster_selection_epsilon are required for method = 'finalize'.")
  }
  grid <- if (method == "finalize") {
    expand.grid(MinPts = final_minPts,
      Epsilon = if (is.null(final_cluster_selection_epsilon)) 0 else final_cluster_selection_epsilon)
  } else expand.grid(MinPts = minPts_range,
    Epsilon = cluster_selection_epsilon_range)

  # This pipeline clusters the SOM codebook with HDBSCAN, so the SOM's own
  # latent-profile grid is not just unused, it is the most expensive step in
  # the run. Handing the node clustering to CreateClusterModel_SOM_MClust()
  # keeps its frozen scaling, SOM fit, node-to-individual mapping, and
  # projection contract while skipping the profile fits entirely.
  ModelInfo_HDBSCAN <- NULL
  ClusterSOMNodes <- function(codes) {
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
    best_i <- if (method == "finalize") 1L else {
      recommended <- which(review$fit_table$Recommended)
      if (length(recommended)) recommended[[1]] else 1L
    }
    node_fit <- fits[[best_i]]
    ModelInfo_HDBSCAN <<- list(
      hdbscan_model = node_fit, node_cluster = node_fit$cluster,
      fit_table = review$fit_table, AHP = review$AHP,
      final_minPts = grid$MinPts[[best_i]],
      final_cluster_selection_epsilon = grid$Epsilon[[best_i]],
      model_status = method)
    list(
      node_cluster = node_fit$cluster,
      fit_table = review$fit_table,
      ahp_best_row = review$AHP$ahp_best_row,
      recommendation = review$AHP$recommendation,
      best_fit_name = paste0(
        "hdbscan_minPts_", grid$MinPts[[best_i]],
        "_epsilon_", grid$Epsilon[[best_i]]),
      fit_plot = PlotClusterFitReview(review$fit_table, x = "MinPts"))
  }

  # Reuse the established, tested SOM fitting and frozen scaling workflow.
  som <- CreateClusterModel_SOM_MClust(
    data = data, variables = variables, method = "exploratory", ...,
    ClusterVariableName = ClusterVariableName, seed_som = seed_som,
    seed_lpa = seed_hdbscan, stability_resamples = 0L,
    .NodeClusterFn = ClusterSOMNodes)
  som$ModelInfo_HDBSCAN <- ModelInfo_HDBSCAN
  som$ModelInfo <- ModelInfo_HDBSCAN
  selected_minPts <- ModelInfo_HDBSCAN$final_minPts
  selected_epsilon <- ModelInfo_HDBSCAN$final_cluster_selection_epsilon
  som$Specification <- .ClusterSpecification("SOM + HDBSCAN", som$vars_used,
    list(seed_som = seed_som, seed_hdbscan = seed_hdbscan),
    list(ZScoreObject = som$ZScoreObject, SOM = som$ModelInfo_SOM$som_grid_info),
    grid, list(minPts = selected_minPts, Epsilon = selected_epsilon),
    !is.na(som$ModelInfo_SOM$SOMFit$table$SOM_Distance),
    list(projection = "frozen SOM node to HDBSCAN label"))
  som$method <- method
  if (method == "exploratory") {
    som$DataWithClusters <- NULL
    som$ProbFit <- NULL
    som$ProjectionFit <- NULL
    som$Stability <- NULL
    class(som) <- c("Pipeline_SOM_HDBSCAN_Exploratory", class(som))
    return(som)
  }
  # The node labels reached the individuals through the shared SOM assembly,
  # so ProbFit$individual and DataWithClusters already carry them.
  individual <- som$ProbFit$individual
  if (stability_resamples > 0L) {
    reference <- individual$Cluster
    Stability <- .ClusterBootstrapStability(data, reference, function(boot, original) {
      fitted <- CreateClusterModel_SOM_HDBSCAN(
        data = boot, variables = variables, method = "finalize",
        final_minPts = selected_minPts,
        final_cluster_selection_epsilon = selected_epsilon,
        ClusterVariableName = ".scidr_stability_cluster",
        seed_som = seed_som, seed_hdbscan = seed_hdbscan,
        stability_resamples = 0L, ...)
      ProjectCluster(fitted, original,
        ClusterVariableName = ".scidr_stability_cluster")$ProbFit$individual$Cluster
    }, resamples = stability_resamples, seed = stability_seed,
    candidate = list(MinPts = selected_minPts, Epsilon = selected_epsilon),
    noise_label = 0L, progress = stability_progress)
    Stability$plots <- .ClusterStabilityPlots(Stability)
    som$Stability <- Stability
    som$ModelInfo_HDBSCAN$fit_table <- dplyr::left_join(
      som$ModelInfo_HDBSCAN$fit_table, Stability$summary,
      by = c("MinPts", "Epsilon"))
    som$ModelInfo <- som$ModelInfo_HDBSCAN
  } else {
    som$Stability <- NULL
  }
  class(som) <- c("Pipeline_SOM_HDBSCAN", class(som))
  som
}

#' Project cases onto a SOM + HDBSCAN model
#'
#' @description Maps cases to the frozen SOM, then assigns the HDBSCAN label
#' of each best-matching node. This is not a refitted or native HDBSCAN
#' prediction; its distance diagnostic is distance to the frozen SOM BMU.
#' @param object A `CreateClusterModel_SOM_HDBSCAN()` result.
#' @param new_df New cases to map to the frozen SOM nodes.
#' @param ClusterVariableName Name of the projected cluster column.
#' @param ... Additional arguments passed to [ProjectCluster()].
#' @return A `Project_SOMHDBSCAN` result. Assignment is inherited from the
#'   projected case's frozen best-matching SOM node; it is not a native
#'   HDBSCAN prediction.
#' @noRd
#' @export
ProjectCluster.Pipeline_SOM_HDBSCAN_Exploratory <- function(object, new_df, ...) {
  stop("ProjectCluster() requires a finalized SOM + HDBSCAN model.", call. = FALSE)
}

#' @noRd
#' @export
ProjectCluster.Pipeline_SOM_HDBSCAN <- function(object, new_df,
    ClusterVariableName = object$ClusterVariableName, ...) {
  # The canonical SOM + Mclust projector accepts the inherited class.
  projected <- ProjectCluster.Pipeline_SOM_MClust(
    object, new_df, ClusterVariableName = ClusterVariableName, ...
  )
  node_cluster <- object$ModelInfo_HDBSCAN$node_cluster
  projected$ProjectionFit$individual$Cluster <- node_cluster[
    projected$ProjectionFit$individual$SOM_Node]
  projected$ProbFit$individual$Cluster <- projected$ProjectionFit$individual$Cluster
  projected$ProjectionFit$individual$Projection_Fit_Class <- dplyr::case_when(
    is.na(projected$ProjectionFit$individual$SOM_Distance) ~ NA_character_,
    projected$ProjectionFit$individual$Cluster == 0 ~ "Potential novel phenotype",
    TRUE ~ "Good fit")
  projected$DataWithClusters[[ClusterVariableName]] <- projected$ProbFit$individual$Cluster
  projected$DataWithClusters$Projection_Fit_Class <-
    projected$ProjectionFit$individual$Projection_Fit_Class
  projected$ProjectionFit$by_cluster <- projected$ProjectionFit$individual %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data$SOM_Distance, na.rm = TRUE),
      median_distance = stats::median(.data$SOM_Distance, na.rm = TRUE),
      mean_distance_percentile = mean(.data$SOMDist_percentile_train, na.rm = TRUE),
      median_distance_percentile = stats::median(.data$SOMDist_percentile_train, na.rm = TRUE),
      prop_high_distance = mean(.data$Flag_SOMDist_overallHigh |
        .data$Flag_SOMDist_clusterHigh, na.rm = TRUE),
      prop_potential_novel = mean(.data$Projection_Fit_Class ==
        "Potential novel phenotype", na.rm = TRUE),
      .groups = "drop")
  probability_plot_data <- projected$ProbFit$individual %>%
    dplyr::filter(!is.na(.data$Cluster), is.finite(.data$prob_assigned)) %>%
    dplyr::mutate(Cluster = factor(.data$Cluster))
  projected$ProbFit$plots$individual_ProbAssignedDensity <-
    .SOMAssignedProbabilityPlot(probability_plot_data)
  projected$ProjectionFit$plots$poor_fit_by_cluster <- ggplot2::ggplot(
    projected$ProjectionFit$by_cluster,
    ggplot2::aes(x = factor(.data$Cluster), y = .data$prop_high_distance)
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "High-distance projected cases by HDBSCAN cluster",
      x = "HDBSCAN cluster", y = "Proportion"
    )
  projected$ProjectionFit$policy$assignment <-
    "frozen SOM best-matching node HDBSCAN label"
  projected$ModelInfo_HDBSCAN <- object$ModelInfo_HDBSCAN
  projected$Specification <- object$Specification
  class(projected) <- c("Project_SOMHDBSCAN", class(projected))
  projected
}
