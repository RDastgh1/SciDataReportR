# Freeze one shared metric for all candidate cuts of a map. Shrinkage is
# spherical in the SOM's existing feature units, not a second scaling pass.
.FitSOMNodeWhitening <- function(codes, covariance_shrinkage) {
  if (!is.numeric(covariance_shrinkage) || length(covariance_shrinkage) != 1L ||
      !is.finite(covariance_shrinkage) || covariance_shrinkage < 0 || covariance_shrinkage > 1) {
    stop("covariance_shrinkage must be a single number between 0 and 1.")
  }
  if (nrow(codes) < 3L || any(!is.finite(codes))) {
    stop("Whitening requires at least three finite SOM node profiles.")
  }
  covariance <- stats::cov(codes)
  average_variance <- mean(diag(covariance))
  if (!is.finite(average_variance) || average_variance <= 0) {
    stop("SOM node codebook has no variation.")
  }
  regularized <- (1 - covariance_shrinkage) * covariance +
    covariance_shrinkage * diag(average_variance, ncol(codes))
  decomposition <- eigen(regularized, symmetric = TRUE)
  tolerance <- max(decomposition$values) * ncol(codes) * .Machine$double.eps
  if (min(decomposition$values) <= tolerance) {
    stop("Node covariance is not positive definite; use positive covariance_shrinkage.")
  }
  transform <- decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values), ncol(codes))
  inverse <- diag(sqrt(decomposition$values), ncol(codes)) %*%
    t(decomposition$vectors)
  list(center = colMeans(codes), raw_covariance = covariance,
    covariance = regularized, transform = transform, inverse_transform = inverse,
    covariance_shrinkage = covariance_shrinkage,
    variables = colnames(codes), fit_scope = "all_training_SOM_nodes_uniform_weight")
}

.ApplySOMNodeWhitening <- function(X, whitening) {
  sweep(X, 2L, whitening$center, "-") %*% whitening$transform
}

.SOMKMeansFeatures <- function(object, data) {
  if (identical(object$ZScoreType, "PreZScored")) {
    X <- as.matrix(data[, object$ZScoreVars, drop = FALSE])
  } else {
    projected <- ProjectZScore(data = data, variables = object$vars_used,
      parameters = object$ZScoreObject, ParameterInputType = "ZScoreObj",
      names_prefix = "Z_")
    X <- as.matrix(projected$ZScores[, object$ZScoreVars, drop = FALSE])
  }
  if (any(is.infinite(X))) stop("SOM K-means feature values must be finite or missing.")
  X
}

# Distances always refer to the inherited cluster, even when the closest
# centroid belongs to another cluster. Negative margins retain that evidence.
.SOMKMeansAssignmentFit <- function(X, centers, assigned) {
  distances <- matrix(NA_real_, nrow(X), nrow(centers))
  for (cluster_index in seq_len(nrow(centers))) {
    distances[, cluster_index] <- sqrt(rowSums(sweep(X, 2L, centers[cluster_index, ], "-")^2))
  }
  valid <- stats::complete.cases(X) & !is.na(assigned)
  nearest <- rep(NA_integer_, nrow(X))
  own <- alternative <- rep(NA_real_, nrow(X))
  if (any(valid)) {
    nearest[valid] <- max.col(-distances[valid, , drop = FALSE], ties.method = "first")
    own[valid] <- distances[cbind(which(valid), assigned[valid])]
    distances[cbind(which(valid), assigned[valid])] <- Inf
    alternative[valid] <- apply(distances[valid, , drop = FALSE], 1L, min)
  }
  dplyr::tibble(DistanceToCentroid = own, AssignmentMargin = alternative - own,
    NearestCentroidCluster = nearest, Flag_ProjectionDisagreement = nearest != assigned)
}

.SOMKMeansAddDiagnostics <- function(individual, X, ModelInfo) {
  whitened <- .ApplySOMNodeWhitening(X, ModelInfo$Whitening)
  fit <- .SOMKMeansAssignmentFit(whitened, ModelInfo$kmeans_model$centers, individual$Cluster)
  for (column in names(fit)) individual[[column]] <- fit[[column]]
  for (column in c("DistanceToCentroid", "AssignmentMargin", "NearestCentroidCluster",
      "Flag_ProjectionDisagreement")) {
    individual[[paste0("Node", column)]] <- ModelInfo$node_fit[[column]][individual$SOM_Node]
  }
  individual
}

#' Fit covariance-aware K-means clusters on a frozen SOM
#'
#' @description Trains the established SOM, then fits K-means to equally
#'   weighted, whitened node profiles. The frozen regularized node covariance
#'   defines a shared Mahalanobis metric for node clustering. SOM training and
#'   best-matching-node assignment retain their original geometry.
#' @inheritParams CreateClusterModel_SOM_Hierarchical
#' @param nstart Number of K-means random starts; defaults to 50.
#' @param covariance_shrinkage Number from zero to one, default 0.05.
#'   Covariance is shrunk toward a spherical matrix with the same average
#'   variance. Zero requires positive-definite node covariance.
#' @param seed_kmeans Seed for node K-means fits; candidate k uses seed_kmeans + k.
#' @param ... SOM preprocessing and map controls passed to
#'   [CreateClusterModel_SOM_MClust()], including `ZScoreType`, `som_xdim`,
#'   `som_ydim`, `som_epochs`, and `Relabel`. Assignment is hard and node
#'   weighting is uniform; mixture and stability controls in `...` are rejected.
#' @inheritSection cluster-stability-output Stability output
#' @return A `Pipeline_SOM_KMeans` object with `ModelInfo`, `ModelInfo_KMeans`,
#'   `ModelInfo_SOM`, `fit_plot`, `Specification`, and `Stability`. Finalized
#'   models include `DataWithClusters` and `ProbFit`. `ModelInfo$Whitening`
#'   freezes the node mean, raw and regularized covariance, and transforms.
#'   `ModelInfo$fit_table` reports whitened-node WSS, silhouette,
#'   Calinski-Harabasz, minimum node and participant counts, and size balance.
#'   `ModelInfo$centers` and profile figures use original SOM feature units;
#'   `kmeans_model$centers` uses whitened coordinates. `cluster_sizes` reports
#'   node and participant counts, including empty nodes.
#' @details Participants inherit their frozen best-matching node's cluster.
#'   `DistanceToCentroid` measures participant distance to that cluster's
#'   centroid in the frozen Mahalanobis metric. `AssignmentMargin` is distance
#'   to the closest alternative centroid minus assigned-centroid distance;
#'   it can be negative. `Flag_ProjectionDisagreement` compares inherited
#'   assignment to nearest-centroid assignment without changing the label.
#'   `Node*` fields describe the corresponding node profile. These measures
#'   are not membership probabilities. SOM distance diagnostics remain in
#'   `ModelInfo_SOM$SOMFit`; centroid diagnostics are in
#'   `ModelInfo$FitDiagnostics`. Stability retrains scaling, the SOM, whitening,
#'   and K-means on participant subsamples with reference grid dimensions.
#'   It accompanies the initial recommendation without changing node labels
#'   after map assembly. External preprocessing supplied by the caller remains
#'   externally fixed during refits.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' model <- CreateClusterModel_SOM_KMeans(
#'   SimulatedPhenotypeData, paste0("Var", 1:6), method = "finalize",
#'   final_k = 4, som_xdim = 6, som_ydim = 6)
#' projected <- ProjectCluster(model, SimulatedPhenotypeData)
#' projected$ProjectionFit$summary
#' projected$ProjectionFit$centroid_fit$summary
#' }
#' @export
CreateClusterModel_SOM_KMeans <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10, final_k = NULL,
    nstart = 50L, covariance_shrinkage = 0.05,
    ClusterVariableName = "Cluster", seed_som = 934521L, seed_kmeans = 93421L,
    stability_resamples = 0L, seed_stability = seed_som + 1L,
    stability_progress = FALSE, stability_cores = NULL, ...) {
  supplied <- list(k_range = !missing(k_range), final_k = !missing(final_k))
  method <- .ClusterMethod(method)
  .ValidateClusterLifecycle(method, supplied["k_range"], supplied["final_k"])
  if (method == "finalize" && length(final_k) != 1L) stop("final_k must be a single integer.")
  if (!is.numeric(nstart) || length(nstart) != 1L || !is.finite(nstart) ||
      nstart < 1 || nstart != floor(nstart)) stop("nstart must be a positive integer.")
  for (seed in list(seed_som, seed_kmeans)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max - 10000L) stop("Seeds must be finite nonnegative integers below the integer limit.")
  }
  stability_resamples <- .ValidateClusterStability(stability_resamples,
    seed_stability, stability_progress, stability_cores)
  som_args <- list(...)
  protected <- c("assignment", "node_weighting", ".NodeClusterFn", ".CandidatePairs",
    "seed_mixture", "seed_lpa", "models", "final_model", "stability_seed", "stability_type",
    "stability_fractions", "stability_candidates", "stability_top_n")
  if (any(names(som_args) %in% protected) || any(grepl("^(mixture_|lpa_)", names(som_args)))) {
    stop("Supply only SOM preprocessing/map controls in ...; node assignment, mixture, and stability settings are managed by this workflow.")
  }
  whitening <- fits <- review <- selected <- whitened <- NULL
  ClusterSOMNodes <- function(codes) {
    ks <- if (method == "finalize") final_k else k_range
    if (!is.numeric(ks) || !length(ks) || anyNA(ks) || any(!is.finite(ks)) ||
        any(ks != floor(ks)) || any(ks < 2L | ks >= nrow(codes))) {
      stop("Cluster counts must be integers between 2 and the number of SOM nodes minus 1.")
    }
    ks <- unique(as.integer(ks))
    whitening <<- .FitSOMNodeWhitening(codes, covariance_shrinkage)
    whitened <<- .ApplySOMNodeWhitening(codes, whitening)
    fits <<- lapply(ks, function(k) {
      set.seed(seed_kmeans + k)
      stats::kmeans(whitened, centers = k, nstart = nstart)
    })
    rows <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
      fit <- fits[[i]]
      dplyr::tibble(Classes = ks[[i]], WSS = fit$tot.withinss,
        BetweenSS = fit$betweenss,
        CalinskiHarabasz = if (fit$tot.withinss > 0) {
          (fit$betweenss / (length(fit$size) - 1L)) /
            (fit$tot.withinss / (nrow(codes) - length(fit$size)))
        } else NA_real_,
        Silhouette = mean(cluster::silhouette(fit$cluster, stats::dist(whitened))[, "sil_width"]),
        MinClusterN = min(fit$size), SizeBalance = min(fit$size) / max(fit$size))
    }))
    review <<- .ClusterAHP(rows,
      higher = c("Silhouette", "CalinskiHarabasz", "MinClusterN", "SizeBalance"),
      setting = "SOM node K-means cluster count")
    selected <<- which(review$fit_table$Recommended)[[1]]
    list(node_cluster = fits[[selected]]$cluster,
      fit_table = review$fit_table, ahp_best_row = review$AHP$ahp_best_row,
      recommendation = review$AHP$recommendation,
      best_fit_name = paste0("som_kmeans_k_", ks[[selected]]),
      fit_plot = PlotClusterFitReview(review$fit_table))
  }
  som <- do.call(CreateClusterModel_SOM_MClust, c(list(data = data,
    variables = variables, method = "exploratory",
    k_range = if (method == "finalize") final_k else k_range,
    ClusterVariableName = ClusterVariableName, seed_som = seed_som,
    assignment = "hard", node_weighting = "uniform", stability_resamples = 0L,
    .NodeClusterFn = ClusterSOMNodes), som_args))
  best <- fits[[selected]]
  codes <- som$ModelInfo_SOM$som_codes
  occupancy <- tabulate(som$ProbFit$individual$SOM_Node, nrow(codes))
  som$ProbFit$node$NodeN <- occupancy
  rows <- review$fit_table
  rows$MinParticipantN <- vapply(fits, function(fit) {
    min(vapply(split(occupancy, fit$cluster), sum, numeric(1)))
  }, numeric(1))
  Stability <- NULL
  if (stability_resamples > 0L) {
    refit_args <- som_args
    refit_args$map_from <- NULL
    refit_args$som_xdim <- som$ModelInfo_SOM$som_grid_info$som_xdim_used
    refit_args$som_ydim <- som$ModelInfo_SOM$som_grid_info$som_ydim_used
    stabilities <- lapply(seq_along(fits), function(i) {
      k <- length(fits[[i]]$size)
      reference <- fits[[i]]$cluster[som$ProbFit$individual$SOM_Node]
      .ClusterSubsampleStability(data, reference, function(boot, model_seed) {
        refitted <- do.call(CreateClusterModel_SOM_KMeans, c(list(
          data = boot, variables = variables, method = "finalize", final_k = k,
          nstart = nstart, covariance_shrinkage = covariance_shrinkage,
          seed_som = model_seed, seed_kmeans = model_seed, stability_resamples = 0L), refit_args))
        dplyr::select(refitted$ProbFit$individual, dplyr::all_of(c(".row_id", "Cluster")))
      }, stability_resamples, seed_stability + k, list(Classes = k),
      progress = stability_progress, stability_cores = stability_cores,
      memory_method = "som_kmeans", som_nodes = nrow(codes), memory_variables = som$vars_used)
    })
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, seed_stability)
    Stability$settings$som_grid <- list(xdim = refit_args$som_xdim,
      ydim = refit_args$som_ydim, policy = "fixed_to_reference_fit")
    Stability$settings$whitening <- "refit_from_subsample_SOM_nodes"
    Stability$plots <- .ClusterStabilityPlots(Stability)
    rows <- dplyr::left_join(rows, Stability$summary, by = "Classes")
  }
  centers <- sweep(best$centers %*% whitening$inverse_transform, 2L, whitening$center, "+")
  colnames(centers) <- colnames(codes)
  labels <- vapply(som$vars_used, function(v) .ClusterVariableLabel(data, v), character(1))
  silhouette <- cluster::silhouette(best$cluster, stats::dist(whitened))
  cluster_sizes <- dplyr::tibble(Cluster = best$cluster, ParticipantN = occupancy) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(NodeN = dplyr::n(), ParticipantN = sum(.data$ParticipantN), .groups = "drop")
  ModelInfo <- list(kmeans_model = best, Whitening = whitening,
    node_cluster = best$cluster, node_weighting = "uniform", centers = centers,
    final_k = length(best$size), fit_table = rows, AHP = review$AHP,
    silhouette = silhouette, cluster_sizes = cluster_sizes,
    node_fit = .SOMKMeansAssignmentFit(whitened, best$centers, best$cluster),
    plots = .DropNullPlots(list(
      elbow = .ClusterCandidateCurve(rows, "WSS", "", "Whitened-node within-cluster sum of squares"),
      silhouette = PlotClusterSilhouette(silhouette, title = NULL, observation_unit = "SOM node"),
      silhouette_by_k = .ClusterCandidateCurve(rows, "Silhouette", "", "Average node silhouette width"),
      calinski_harabasz = .ClusterCandidateCurve(rows, "CalinskiHarabasz", ""),
      centre_heatmap = PlotClusterCentreHeatmap(centers, labels, "", "SOM codebook units"),
      centre_profile = PlotClusterCentreProfile(centers, labels, "", "SOM codebook units"),
      profiles = PlotClusterProfiles(som$DataWithClusters, ClusterVariableName, som$vars_used))))
  individual <- .SOMKMeansAddDiagnostics(.SOMHardClusterFitClass(som$ProbFit$individual),
    .SOMKMeansFeatures(som, data), ModelInfo)
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(individual,
    "DistanceToCentroid", "Mahalanobis distance to inherited cluster centroid")
  ModelInfo$FitDiagnostics$projection_agreement <- mean(!individual$Flag_ProjectionDisagreement, na.rm = TRUE)
  som$ProbFit$individual <- individual
  som$ProbFit$plots <- .ClusterConfidencePlots(individual, "AssignmentMargin",
    "Mahalanobis distance margin to alternative centroid", bounded = FALSE)
  som$ModelInfo_SOM$SOMFit$table <- individual
  som$ModelInfo <- som$ModelInfo_KMeans <- ModelInfo
  som$ModelInfo_MClust <- NULL
  som$method <- method
  som$fit_plot <- PlotClusterFitReview(rows)
  som$Stability <- Stability
  som$Preprocessing <- som$Specification$preprocessing
  som$Preprocessing$Whitening <- whitening
  som$Specification$method <- "SOM + KMeans"
  som$Specification$candidate_grid <- dplyr::select(rows, "Classes")
  som$Specification$selected <- list(k = ModelInfo$final_k, nstart = nstart,
    covariance_shrinkage = covariance_shrinkage)
  som$Specification$preprocessing$Whitening <- whitening
  som$Specification$seeds$seed_mixture <- NULL
  som$Specification$seeds$seed_lpa <- NULL
  som$Specification$seeds$seed_kmeans <- seed_kmeans
  som$Specification$seeds$seed_stability <- seed_stability
  som$Specification$seeds$stability_seed <- seed_stability
  som$Specification$stability <- if (is.null(Stability)) NULL else Stability$settings
  som$Specification$projection$membership_probability <- "unavailable"
  som$Specification$support$node_metric <- "frozen regularized Mahalanobis metric; equal node weights"
  som$Specification$selection <- list(
    recommendation = "initial node-fit metrics; stability reported alongside",
    covariance_fit_scope = whitening$fit_scope, node_weighting = "uniform")
  if (method == "exploratory") {
    som$DataWithClusters <- NULL
    som$ProbFit <- NULL
  }
  class(som) <- c("Pipeline_SOM_KMeans", class(som))
  som
}

#' @noRd
#' @export
ProjectCluster.Pipeline_SOM_KMeans <- function(object, data,
    ClusterVariableName = object$ClusterVariableName, ...) {
  if (object$method != "finalize") stop("ProjectCluster() requires a finalized SOM K-means model.")
  projected <- ProjectCluster.Pipeline_SOM_MClust(object, data,
    ClusterVariableName = ClusterVariableName, ...)
  individual <- .SOMKMeansAddDiagnostics(
    .SOMHardClusterFitClass(projected$ProjectionFit$individual),
    .SOMKMeansFeatures(object, data), object$ModelInfo)
  projected$ProjectionFit$individual <- individual
  projected$ProjectionFit$fit_class <- individual$Projection_Fit_Class
  projected$ProbFit$individual <- individual
  projected$ProbFit$plots <- .ClusterConfidencePlots(individual, "AssignmentMargin",
    "Mahalanobis distance margin to alternative centroid", prefix = "Projected", bounded = FALSE)
  projected$ProjectionFit$centroid_fit <- .ClusterProjectionFit(individual,
    object$ModelInfo$FitDiagnostics$distance_reference,
    "DistanceToCentroid", "Mahalanobis distance to inherited cluster centroid",
    confidence_var = "AssignmentMargin", confidence_threshold = 0)
  projected$ProjectionFit$centroid_fit$policy$assignment <- "inherited frozen SOM node label; centroid distance does not reassign"
  projected$ProjectionFit$centroid_fit$policy$training_projection_agreement <-
    object$ModelInfo$FitDiagnostics$projection_agreement
  projected$ProjectionFit$plots$centroid_distance_hist <-
    projected$ProjectionFit$centroid_fit$plots$distance_hist
  projected$ProjectionFit$plots$projection_fit_class_bar <- .ClusterFitClassPlot(individual$Projection_Fit_Class)
  projected$ProjectionFit$by_cluster$prop_poor_fit <- projected$ProjectionFit$by_cluster$prop_high_distance
  projected$ProjectionFit$by_cluster$prop_potential_novel <- NA_real_
  projected$ProjectionFit$policy$assignment <- "frozen SOM best-matching node K-means label"
  projected$ProjectionFit$policy$membership_probability <- "unavailable; margins are not posterior probabilities"
  projected$ProjectionFit$policy$low_probability_threshold <- NULL
  for (column in c("Projection_Fit_Class", "FitClass", "FitClass_Legacy")) {
    projected$DataWithClusters[[column]] <- individual[[column]]
  }
  projected$ModelInfo <- projected$ModelInfo_KMeans <- object$ModelInfo
  projected$ModelInfo_MClust <- NULL
  projection_specification <- projected$Specification$projection
  projected$Specification <- object$Specification
  projected$Specification$projection <- utils::modifyList(
    projected$Specification$projection, projection_specification)
  projected$Specification$projection$membership_probability <- "unavailable"
  class(projected) <- c("Project_SOM_KMeans", class(projected))
  projected
}
