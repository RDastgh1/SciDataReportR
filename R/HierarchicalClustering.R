.ValidateHierarchicalSettings <- function(linkage, distance, ks, n) {
  linkage <- match.arg(linkage, c("ward.D2", "average", "complete", "single"))
  distance <- match.arg(distance, c("euclidean", "manhattan"))
  if (linkage == "ward.D2" && distance != "euclidean") {
    stop("Ward linkage requires Euclidean distance.", call. = FALSE)
  }
  if (!is.numeric(ks) || !length(ks) || anyNA(ks) ||
      any(!is.finite(ks)) || any(ks != floor(ks)) ||
      any(ks < 2L | ks >= n)) {
    stop("Cluster counts must be integers between 2 and the number of complete rows/nodes minus 1.",
      call. = FALSE)
  }
  list(linkage = linkage, distance = distance, ks = unique(as.integer(ks)))
}

# A single distance object and tree serve every candidate cut. Medoids are
# actual training profiles; they define an explicit frozen projection rule.
.FitHierarchicalCuts <- function(X, ks, linkage, distance) {
  settings <- .ValidateHierarchicalSettings(linkage, distance, ks, nrow(X))
  if (any(!is.finite(X))) stop("Clustering inputs must be finite after preprocessing.")
  distances <- stats::dist(X, method = settings$distance)
  tree <- stats::hclust(distances, method = settings$linkage)
  distance_matrix <- as.matrix(distances)
  fits <- lapply(settings$ks, function(k) {
    cluster <- unname(stats::cutree(tree, k = k))
    groups <- split(seq_len(nrow(X)), cluster)
    medoid_rows <- vapply(groups, function(rows) {
      rows[[which.min(rowSums(distance_matrix[rows, rows, drop = FALSE]))]]
    }, integer(1))
    centers <- do.call(rbind, lapply(groups, function(rows) colMeans(X[rows, , drop = FALSE])))
    within_ss <- sum(vapply(groups, function(rows) {
      sum(sweep(X[rows, , drop = FALSE], 2L, colMeans(X[rows, , drop = FALSE]), "-")^2)
    }, numeric(1)))
    total_ss <- sum(sweep(X, 2L, colMeans(X), "-")^2)
    silhouette <- cluster::silhouette(cluster, distances)
    list(cluster = cluster, medoid_rows = medoid_rows,
      medoids = X[medoid_rows, , drop = FALSE], centers = centers,
      silhouette = silhouette,
      row = dplyr::tibble(Classes = k, Silhouette = mean(silhouette[, "sil_width"]),
        WSS = if (distance == "euclidean") within_ss else NA_real_,
        CalinskiHarabasz = if (distance == "euclidean" && within_ss > 0) {
          ((total_ss - within_ss) / (k - 1L)) / (within_ss / (nrow(X) - k))
        } else NA_real_,
        MinClusterN = min(lengths(groups)),
        SizeBalance = min(lengths(groups)) / max(lengths(groups))))
  })
  list(tree = tree, fits = fits, rows = dplyr::bind_rows(lapply(fits, `[[`, "row")),
    settings = settings)
}

.HierarchicalMedoidDistances <- function(X, medoids, distance) {
  result <- matrix(NA_real_, nrow(X), nrow(medoids))
  for (cluster_index in seq_len(nrow(medoids))) {
    differences <- sweep(X, 2L, medoids[cluster_index, ], "-")
    result[, cluster_index] <- if (distance == "euclidean") {
      sqrt(rowSums(differences^2))
    } else rowSums(abs(differences))
  }
  result
}

.HierarchicalAssignmentFit <- function(X, medoids, distance, assigned = NULL) {
  distances <- .HierarchicalMedoidDistances(X, medoids, distance)
  nearest <- max.col(-distances, ties.method = "first")
  if (is.null(assigned)) assigned <- nearest
  own <- distances[cbind(seq_len(nrow(X)), assigned)]
  other <- distances
  other[cbind(seq_len(nrow(X)), assigned)] <- Inf
  alternative <- apply(other, 1L, min)
  # A negative training margin is meaningful: the dendrogram assigned this
  # participant to a different group than nearest-medoid projection would.
  dplyr::tibble(DistanceToMedoid = own,
    AssignmentMargin = alternative - own,
    NearestMedoidCluster = nearest,
    Flag_ProjectionDisagreement = nearest != assigned)
}

.HierarchicalDendrogram <- function(tree, cluster) {
  dendrogram <- stats::as.dendrogram(tree)
  segments <- list()
  leaf <- 0L
  Traverse <- function(branch) {
    if (stats::is.leaf(branch)) {
      leaf <<- leaf + 1L
      return(c(x = leaf, y = 0))
    }
    children <- lapply(branch, Traverse)
    parent <- c(x = mean(vapply(children, `[[`, numeric(1), "x")),
      y = attr(branch, "height"))
    for (child in children) {
      segments[[length(segments) + 1L]] <<- data.frame(
        x = c(child[["x"]], child[["x"]]),
        y = c(child[["y"]], parent[["y"]]),
        xend = c(child[["x"]], parent[["x"]]),
        yend = c(parent[["y"]], parent[["y"]]))
    }
    parent
  }
  Traverse(dendrogram)
  df_Leaves <- data.frame(x = seq_along(tree$order), y = 0,
    Cluster = factor(cluster[tree$order]))
  ggplot2::ggplot(dplyr::bind_rows(segments)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y,
      xend = .data$xend, yend = .data$yend), linewidth = 0.3) +
    ggplot2::geom_point(data = df_Leaves,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$Cluster), size = 1.5) +
    ggplot2::theme_bw() + ggplot2::labs(x = "Dendrogram order", y = "Merge height")
}

.HierarchicalReview <- function(rows) {
  .ClusterAHP(rows,
    higher = c("Silhouette", "CalinskiHarabasz", "MinClusterN", "SizeBalance",
      "ReproducibilityScore"), setting = "Hierarchical cluster count")
}

.HierarchicalSelected <- function(rows) {
  recommended <- which(rows$Recommended)
  if (length(recommended)) recommended[[1]] else which.max(rows$Silhouette)
}

.HierarchicalPlots <- function(fitted, selected, rows, variable_labels, scale_label) {
  best <- fitted$fits[[selected]]
  .DropNullPlots(list(
    dendrogram = .HierarchicalDendrogram(fitted$tree, best$cluster),
    silhouette = PlotClusterSilhouette(best$silhouette, title = NULL),
    silhouette_by_k = .ClusterCandidateCurve(rows, "Silhouette", "", "Average silhouette width"),
    calinski_harabasz = .ClusterCandidateCurve(rows, "CalinskiHarabasz", ""),
    centre_heatmap = PlotClusterCentreHeatmap(best$centers, variable_labels, "", scale_label),
    centre_profile = PlotClusterCentreProfile(best$centers, variable_labels, "", scale_label)))
}

#' Fit a hierarchical clustering model with frozen medoid projection
#'
#' @description Clusters numeric participant profiles with an agglomerative
#'   dendrogram. Exploratory mode compares fixed-k cuts of one tree; finalized
#'   models project new participants to the nearest frozen training medoid.
#'   Projection is an approximation to the training partition, not a native
#'   prediction from hierarchical clustering.
#' @param data Data frame containing numeric clustering variables.
#' @param variables Character vector of clustering variables.
#' @param method Either `"exploratory"` or `"finalize"`.
#' @param k_range Candidate cluster counts in exploratory mode.
#' @param final_k Explicit cluster count in finalize mode.
#' @param linkage Agglomeration linkage: `"ward.D2"` (default), `"average"`,
#'   `"complete"`, or `"single"`.
#' @param distance `"euclidean"` (default) or `"manhattan"`. Ward linkage
#'   requires Euclidean distance.
#' @param ZScoreType Frozen numeric preprocessing; defaults to center and scale.
#' @param ClusterVariableName Name of the appended cluster column.
#' @param stability_resamples Number of full-pipeline 90% participant subsample refits.
#' @param seed_stability Seed for participant subsampling.
#' @param stability_progress Whether to print stability progress.
#' @param stability_cores Worker count, or `NULL` for the shared memory-aware policy.
#' @inheritSection cluster-stability-output Stability output
#' @return A `Pipeline_Hierarchical` object with `ModelInfo`,
#'   `ModelInfo_Hierarchical`, `fit_plot`, `Preprocessing`, `Specification`,
#'   and `Stability`. `ModelInfo` contains the tree, fit table, medoids,
#'   silhouette, profile plots, dendrogram, and frozen distance diagnostics.
#'   Finalized models also contain `DataWithClusters` and `ProbFit`.
#'   `ProbFit` holds distances and margins, not posterior probabilities.
#' @details Training labels always come from the dendrogram cut. Training
#'   `AssignmentMargin` compares the assigned medoid against the closest other
#'   medoid and can be negative. `Flag_ProjectionDisagreement` identifies cases
#'   whose nearest medoid differs from their dendrogram label. Projection uses
#'   frozen scaling and medoids; fit flags use the frozen training distance
#'   reference. These are descriptive triage measures, not calibrated tests of
#'   cluster membership or novelty. Pairwise training distances require
#'   quadratic memory in the number of complete participants.
#' @examples
#' df_Example <- data.frame(x = c(1:6, 21:26), y = c(2:7, 22:27))
#' model <- CreateClusterModel_Hierarchical(
#'   df_Example, c("x", "y"), method = "finalize", final_k = 2)
#' projected <- ProjectCluster(model, df_Example)
#' projected$ProjectionFit$summary
#' @export
CreateClusterModel_Hierarchical <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10, final_k = NULL,
    linkage = "ward.D2", distance = "euclidean", ZScoreType = NULL,
    ClusterVariableName = "Cluster", stability_resamples = 0L,
    seed_stability = 93422L, stability_progress = FALSE, stability_cores = NULL) {
  supplied <- list(k_range = !missing(k_range), final_k = !missing(final_k))
  method <- .ClusterMethod(method)
  .ValidateClusterLifecycle(method, supplied["k_range"], supplied["final_k"])
  if (method == "finalize" && is.null(final_k)) stop("final_k is required for method = 'finalize'.")
  if (method == "finalize" && length(final_k) != 1L) stop("final_k must be a single integer.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, seed_stability, stability_progress, stability_cores)
  scaling <- .ResolveNumericClusterScaling(ZScoreType)
  prep <- .PrepareClusterNumeric(data, variables, scaling)
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  ks <- if (method == "finalize") final_k else k_range
  fitted <- .FitHierarchicalCuts(X, ks, linkage, distance)
  rows <- fitted$rows
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(seq_along(fitted$fits), function(i) {
      k <- fitted$settings$ks[[i]]
      reference <- rep(NA_integer_, nrow(prep$data))
      reference[prep$complete_rows] <- fitted$fits[[i]]$cluster
      .ClusterSubsampleStability(prep$data, reference, function(boot, model_seed) {
        refitted <- CreateClusterModel_Hierarchical(boot, prep$variables,
          method = "finalize", final_k = k, linkage = linkage, distance = distance,
          ZScoreType = scaling, stability_resamples = 0L)
        dplyr::select(refitted$ProbFit$individual, dplyr::all_of(c(".row_id", "Cluster")))
      }, stability_resamples, seed_stability + k, list(Classes = k),
      progress = stability_progress, stability_cores = stability_cores,
      memory_method = "hierarchical", memory_variables = prep$variables)
    })
    rows <- dplyr::left_join(rows,
      dplyr::bind_rows(lapply(stabilities, `[[`, "summary")), by = "Classes")
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, seed_stability)
    Stability$plots <- .ClusterStabilityPlots(Stability)
  }
  review <- .HierarchicalReview(rows)
  rows <- review$fit_table
  selected <- .HierarchicalSelected(rows)
  best <- fitted$fits[[selected]]
  df_Fit <- .HierarchicalAssignmentFit(X, best$medoids, distance, best$cluster)
  individual <- df_Fit[rep(NA_integer_, nrow(prep$data)), , drop = FALSE]
  individual[prep$complete_rows, ] <- df_Fit
  labels <- vapply(prep$variables, function(v) .ClusterVariableLabel(prep$data, v), character(1))
  ModelInfo <- list(hclust_model = fitted$tree, fit_table = rows, AHP = review$AHP,
    final_k = length(best$medoid_rows), linkage = linkage, distance = distance,
    medoids = best$medoids, medoid_row_ids = prep$data$.row_id[prep$complete_rows][best$medoid_rows],
    centers = best$centers, silhouette = best$silhouette,
    ReviewSpace = .ClusterReviewSpace(X))
  base <- .ClusterOutput(prep, best$cluster, ClusterVariableName, individual, ModelInfo, Stability)
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    base$ProbFit$individual, "DistanceToMedoid", "Distance to assigned medoid")
  ModelInfo$FitDiagnostics$projection_agreement <- mean(!df_Fit$Flag_ProjectionDisagreement)
  ModelInfo$FitDiagnostics$training_distances <- df_Fit$DistanceToMedoid
  ModelInfo$plots <- .HierarchicalPlots(fitted, selected, rows, labels, .ClusterCentreScaleLabel(scaling))
  ModelInfo$plots$profiles <- PlotClusterProfiles(base$DataWithClusters, ClusterVariableName, prep$variables)
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(base$ProbFit$individual, "AssignmentMargin",
    "Distance margin to the closest other medoid", bounded = FALSE)
  out <- c(list(method = method, vars_used = prep$variables,
    ClusterVariableName = ClusterVariableName,
    Preprocessing = list(ZScoreType = scaling, Scaling = scaling, ZScoreObject = prep$ZScoreObject),
    fit_plot = PlotClusterFitReview(rows)), base)
  out$ModelInfo_Hierarchical <- ModelInfo
  out$Specification <- .ClusterSpecification("Hierarchical", prep$variables,
    list(seed_stability = seed_stability), out$Preprocessing,
    dplyr::select(rows, "Classes"),
    list(k = ModelInfo$final_k, linkage = linkage, distance = distance), prep$complete_rows,
    list(projection = "nearest frozen training medoid; approximate assignment"))
  out$Specification$stability <- if (is.null(Stability)) NULL else Stability$settings
  out$Specification$support$training_ranges <- dplyr::tibble(
    Variable = prep$variables,
    Minimum = vapply(prep$data[prep$variables], function(x) {
      if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
    }, numeric(1)),
    Maximum = vapply(prep$data[prep$variables], function(x) {
      if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
    }, numeric(1)))
  if (method == "exploratory") {
    out$DataWithClusters <- NULL
    out$ProbFit <- NULL
  }
  class(out) <- c("Pipeline_Hierarchical", class(out))
  out
}

#' @noRd
#' @export
ProjectCluster.Pipeline_Hierarchical <- function(object, data,
    ClusterVariableName = object$ClusterVariableName, ...) {
  if (object$method != "finalize") stop("ProjectCluster() requires a finalized hierarchical model.")
  prep <- .ProjectClusterNumeric(object, data)
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  if (any(!is.finite(X))) stop("Projection inputs must be finite after preprocessing.")
  individual <- dplyr::tibble(DistanceToMedoid = rep(NA_real_, nrow(prep$data)),
    AssignmentMargin = rep(NA_real_, nrow(prep$data)),
    NearestMedoidCluster = rep(NA_integer_, nrow(prep$data)),
    Flag_ProjectionDisagreement = rep(NA, nrow(prep$data)))
  if (nrow(X)) individual[prep$complete_rows, ] <- .HierarchicalAssignmentFit(
    X, object$ModelInfo$medoids, object$ModelInfo$distance)
  # There is no dendrogram label for a new case against which to compare.
  individual$Flag_ProjectionDisagreement <- NA
  cluster <- individual$NearestMedoidCluster
  base <- .ClusterOutput(prep, cluster[prep$complete_rows], ClusterVariableName,
    individual, object$ModelInfo)
  fit <- .ClusterProjectionFit(base$ProbFit$individual,
    object$ModelInfo$FitDiagnostics$distance_reference,
    "DistanceToMedoid", "Distance to assigned medoid",
    confidence_var = "AssignmentMargin", confidence_threshold = 0)
  fit$policy$assignment <- "nearest frozen training medoid; approximate assignment"
  fit$policy$margin <- "distance to nearest other medoid minus assigned-medoid distance; not a probability"
  fit$policy$training_projection_agreement <- object$ModelInfo$FitDiagnostics$projection_agreement
  training_distances <- object$ModelInfo$FitDiagnostics$training_distances
  fit$individual$DistancePercentileTrain <- vapply(individual$DistanceToMedoid, function(x) {
    if (is.finite(x)) mean(training_distances <= x) else NA_real_
  }, numeric(1))
  fit$individual$Flag_HighDistance <- ifelse(is.na(individual$DistanceToMedoid), NA,
    individual$DistanceToMedoid > fit$policy$high_distance_cutoff)
  ranges <- object$Specification$support$training_ranges
  fit$out_of_support <- dplyr::bind_rows(lapply(seq_len(nrow(ranges)), function(i) {
    values <- prep$data[[ranges$Variable[[i]]]]
    dplyr::tibble(Variable = ranges$Variable[[i]],
      Label = .ClusterVariableLabel(prep$data, ranges$Variable[[i]]),
      TrainingMinimum = ranges$Minimum[[i]], TrainingMaximum = ranges$Maximum[[i]],
      BelowTrainingRange = sum(values < ranges$Minimum[[i]], na.rm = TRUE),
      AboveTrainingRange = sum(values > ranges$Maximum[[i]], na.rm = TRUE))
  })) %>% dplyr::filter(.data$BelowTrainingRange > 0 | .data$AboveTrainingRange > 0)
  fit$plots$separation_map <- .ClusterSeparationMap(object$ModelInfo$ReviewSpace,
    X, cluster[prep$complete_rows], object$ModelInfo$medoids, prefix = "Projected")
  base$ProbFit <- .ClusterProbFit(fit$individual, "AssignmentMargin",
    "Distance margin to the closest other medoid", prefix = "Projected", bounded = FALSE)
  base$DataWithClusters$Projection_Fit_Class <- fit$fit_class
  out <- c(list(vars_used = object$vars_used, ClusterVariableName = ClusterVariableName,
    Preprocessing = object$Preprocessing, Specification = object$Specification,
    ProjectionFit = fit), base)
  out$ModelInfo_Hierarchical <- object$ModelInfo
  class(out) <- c("Project_Hierarchical", class(out))
  out
}

#' Fit hierarchical clusters on a frozen self-organizing map
#'
#' @description Trains the shared SOM workflow and clusters all node codebook
#'   profiles with equal weight. Participants inherit their best-matching node's
#'   dendrogram label. New participants map to the frozen SOM without refitting.
#' @inheritParams CreateClusterModel_Hierarchical
#' @param seed_som Seed for SOM training.
#' @param ... Additional SOM controls passed to [CreateClusterModel_SOM_MClust()],
#'   including `som_xdim`, `som_ydim`, `som_epochs`, `ZScoreType`, and `Relabel`.
#'   Assignment is always hard and node weighting is always uniform.
#' @inheritSection cluster-stability-output Stability output
#' @return A `Pipeline_SOM_Hierarchical` object with the shared SOM outputs,
#'   `ModelInfo_Hierarchical`, and `ModelInfo` containing the node dendrogram,
#'   candidate fit table, medoids, profiles, silhouette, and cluster sizes.
#'   Candidate `MinClusterN` counts nodes; `MinParticipantN` counts people.
#'   `cluster_sizes` reports both for the selected cut, including empty nodes.
#'   Finalized objects retain `DataWithClusters`, `ProbFit`, and SOM distance
#'   diagnostics. Probability fields inherited from the SOM infrastructure are
#'   missing because hierarchical clustering has no posterior probabilities.
#' @details Node occupancy is reported, not used to weight the dendrogram.
#'   Stability refits scaling, the SOM, and the dendrogram on 90% participant
#'   subsamples using the reference grid dimensions. It compares sampled
#'   participant assignments. Stability accompanies the initial candidate
#'   recommendation without changing the selected node labels after map assembly.
#'   Projection retains frozen-map distance percentiles,
#'   high-distance flags, range checks, residuals, and cohort summaries. Node
#'   medoid margins describe node separation, not participant probabilities.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' model <- CreateClusterModel_SOM_Hierarchical(
#'   SimulatedPhenotypeData, paste0("Var", 1:6),
#'   method = "finalize", final_k = 3, som_xdim = 4, som_ydim = 4)
#' projected <- ProjectCluster(model, SimulatedPhenotypeData)
#' projected$ProjectionFit$summary
#' }
#' @export
CreateClusterModel_SOM_Hierarchical <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10, final_k = NULL,
    linkage = "ward.D2", distance = "euclidean",
    ClusterVariableName = "Cluster", seed_som = 934521L,
    stability_resamples = 0L, seed_stability = seed_som + 1L,
    stability_progress = FALSE, stability_cores = NULL, ...) {
  supplied <- list(k_range = !missing(k_range), final_k = !missing(final_k))
  method <- .ClusterMethod(method)
  .ValidateClusterLifecycle(method, supplied["k_range"], supplied["final_k"])
  if (method == "finalize" && is.null(final_k)) stop("final_k is required for method = 'finalize'.")
  if (method == "finalize" && length(final_k) != 1L) stop("final_k must be a single integer.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, seed_stability, stability_progress, stability_cores)
  som_args <- list(...)
  if (any(c("assignment", "node_weighting", ".NodeClusterFn", ".CandidatePairs",
      "seed_mixture", "models", "final_model", "stability_seed", "stability_type",
      "stability_fractions") %in% names(som_args))) {
    stop("Supply only SOM preprocessing/map controls in ...; assignment, node weighting, mixture, and stability controls are managed by this workflow.")
  }
  fitted <- NULL
  review <- NULL
  selected <- NULL
  ClusterSOMNodes <- function(codes) {
    fitted <<- .FitHierarchicalCuts(codes,
      if (method == "finalize") final_k else k_range, linkage, distance)
    review <<- .HierarchicalReview(fitted$rows)
    selected <<- .HierarchicalSelected(review$fit_table)
    list(node_cluster = fitted$fits[[selected]]$cluster,
      fit_table = review$fit_table, ahp_best_row = review$AHP$ahp_best_row,
      recommendation = review$AHP$recommendation,
      best_fit_name = paste0("hierarchical_k_", fitted$settings$ks[[selected]]),
      fit_plot = PlotClusterFitReview(review$fit_table))
  }
  som <- do.call(CreateClusterModel_SOM_MClust, c(list(
    data = data, variables = variables, method = "exploratory",
    k_range = if (method == "finalize") final_k else k_range,
    ClusterVariableName = ClusterVariableName, seed_som = seed_som,
    assignment = "hard", node_weighting = "uniform", stability_resamples = 0L,
    .NodeClusterFn = ClusterSOMNodes), som_args))
  best <- fitted$fits[[selected]]
  occupancy <- tabulate(som$ProbFit$individual$SOM_Node, nrow(som$ModelInfo_SOM$som_codes))
  som$ProbFit$node$NodeN <- occupancy
  rows <- review$fit_table
  rows$MinParticipantN <- vapply(fitted$fits, function(fit) {
    min(vapply(split(occupancy, fit$cluster), sum, numeric(1)))
  }, numeric(1))
  Stability <- NULL
  if (stability_resamples > 0L) {
    refit_args <- som_args
    refit_args$map_from <- NULL
    refit_args$som_xdim <- som$ModelInfo_SOM$som_grid_info$som_xdim_used
    refit_args$som_ydim <- som$ModelInfo_SOM$som_grid_info$som_ydim_used
    stabilities <- lapply(seq_along(fitted$fits), function(i) {
      k <- fitted$settings$ks[[i]]
      reference <- fitted$fits[[i]]$cluster[som$ProbFit$individual$SOM_Node]
      .ClusterSubsampleStability(data, reference, function(boot, model_seed) {
        refitted <- do.call(CreateClusterModel_SOM_Hierarchical, c(list(
          data = boot, variables = variables, method = "finalize", final_k = k,
          linkage = linkage, distance = distance, seed_som = model_seed,
          stability_resamples = 0L), refit_args))
        dplyr::select(refitted$ProbFit$individual, dplyr::all_of(c(".row_id", "Cluster")))
      }, stability_resamples, seed_stability + k, list(Classes = k),
      progress = stability_progress, stability_cores = stability_cores,
      memory_method = "som_hierarchical", som_nodes = length(occupancy),
      memory_variables = som$vars_used)
    })
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, seed_stability)
    Stability$settings$som_grid <- list(xdim = refit_args$som_xdim,
      ydim = refit_args$som_ydim, policy = "fixed_to_reference_fit")
    Stability$plots <- .ClusterStabilityPlots(Stability)
    rows <- dplyr::left_join(rows, Stability$summary, by = "Classes")
    # The map and assignments were assembled for the original candidate.
    # Keep that candidate fixed; stability supplies review evidence rather
    # than silently changing the fitted node labels after map assembly.
  }
  node_fit <- .HierarchicalAssignmentFit(som$ModelInfo_SOM$som_codes,
    best$medoids, distance, best$cluster)
  node_fit$NodeID <- seq_along(best$cluster)
  cluster_sizes <- dplyr::tibble(Cluster = best$cluster, ParticipantN = occupancy) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(NodeN = dplyr::n(), ParticipantN = sum(.data$ParticipantN), .groups = "drop")
  labels <- vapply(som$vars_used, function(v) .ClusterVariableLabel(data, v), character(1))
  ModelInfo <- list(hclust_model = fitted$tree, fit_table = rows, AHP = review$AHP,
    final_k = length(best$medoid_rows), linkage = linkage, distance = distance,
    node_cluster = best$cluster, medoids = best$medoids,
    medoid_nodes = best$medoid_rows, centers = best$centers,
    silhouette = best$silhouette, cluster_sizes = cluster_sizes,
    node_fit = node_fit, node_weighting = "uniform",
    FitDiagnostics = som$ModelInfo_SOM$SOMFit,
    plots = .HierarchicalPlots(fitted, selected, rows, labels, "SOM codebook units"))
  ModelInfo$plots$silhouette <- PlotClusterSilhouette(best$silhouette,
    title = NULL, observation_unit = "SOM node")
  ModelInfo$plots$profiles <- PlotClusterProfiles(som$DataWithClusters,
    ClusterVariableName, som$vars_used)
  som$ModelInfo <- ModelInfo
  som$ModelInfo_Hierarchical <- ModelInfo
  som$ModelInfo_MClust <- NULL
  som$fit_plot <- PlotClusterFitReview(rows)
  som$Stability <- Stability
  som$method <- method
  som$Specification$method <- "SOM + Hierarchical"
  som$Specification$candidate_grid <- dplyr::select(rows, "Classes")
  som$Specification$selected <- list(k = ModelInfo$final_k, linkage = linkage, distance = distance)
  som$Specification$stability <- if (is.null(Stability)) NULL else Stability$settings
  som$Specification$projection$assignment <- "hard"
  som$Specification$support$node_clustering <- "equal-weight SOM nodes; frozen node dendrogram labels"
  som$Specification$seeds$seed_stability <- seed_stability
  som$ProbFit$individual <- .SOMHardClusterFitClass(som$ProbFit$individual)
  som$ModelInfo_SOM$SOMFit$table <- som$ProbFit$individual
  som$ModelInfo$FitDiagnostics$table <- som$ProbFit$individual
  som$ModelInfo_Hierarchical$FitDiagnostics$table <- som$ProbFit$individual
  for (column in c("NodeDistanceToMedoid", "NodeAssignmentMargin", "NodeFlag_ProjectionDisagreement")) {
    source <- sub("^Node", "", column)
    som$ProbFit$individual[[column]] <- node_fit[[source]][som$ProbFit$individual$SOM_Node]
  }
  som$ProbFit$plots <- .ClusterConfidencePlots(som$ProbFit$individual,
    "NodeAssignmentMargin", "SOM node margin to the closest other cluster medoid",
    prefix = "Training", bounded = FALSE)
  if (method == "exploratory") {
    som$DataWithClusters <- NULL
    som$ProbFit <- NULL
    som$ProjectionFit <- NULL
  }
  class(som) <- c("Pipeline_SOM_Hierarchical", class(som))
  som
}

# Hard node clustering provides no membership posterior. Classify map fit
# from distance alone rather than treating absent probabilities as uncertainty.
.SOMHardClusterFitClass <- function(individual) {
  far <- individual$Flag_SOMDist_overallHigh | individual$Flag_SOMDist_clusterHigh
  individual$Projection_Fit_Class <- ifelse(is.na(individual$SOM_Distance), NA_character_,
    ifelse(far, "Poor fit to training structure", "Good fit"))
  individual$FitClass <- factor(individual$Projection_Fit_Class,
    levels = c("Good fit", "Poor fit to training structure"))
  if ("FitClass_Legacy" %in% names(individual)) {
    individual$FitClass_Legacy <- individual$Projection_Fit_Class
  }
  individual
}

#' @noRd
#' @export
ProjectCluster.Pipeline_SOM_Hierarchical <- function(object, data,
    ClusterVariableName = object$ClusterVariableName, ...) {
  if (object$method != "finalize") stop("ProjectCluster() requires a finalized SOM hierarchical model.")
  projected <- ProjectCluster.Pipeline_SOM_MClust(object, data,
    ClusterVariableName = ClusterVariableName, ...)
  # Preserve all frozen-map fit flags, range checks, residuals, and summaries.
  # Missing mixture probabilities must not be interpreted as low certainty.
  projected$ProjectionFit$policy$assignment <- "frozen SOM best-matching node hierarchical label"
  projected$ProjectionFit$policy$membership_probability <- "unavailable; hierarchical labels have no posterior probabilities"
  projected$ProjectionFit$policy$low_probability_threshold <- NULL
  individual <- .SOMHardClusterFitClass(projected$ProjectionFit$individual)
  projected$ProjectionFit$individual <- individual
  projected$ProjectionFit$fit_class <- individual$Projection_Fit_Class
  for (column in c("Projection_Fit_Class", "FitClass", "FitClass_Legacy")) {
    projected$ProbFit$individual[[column]] <- individual[[column]]
    projected$DataWithClusters[[column]] <- individual[[column]]
  }
  projected$ProjectionFit$plots$projection_fit_class_bar <-
    .ClusterFitClassPlot(individual$Projection_Fit_Class)
  projected$ProjectionFit$by_cluster$prop_poor_fit <-
    projected$ProjectionFit$by_cluster$prop_high_distance
  projected$ProjectionFit$by_cluster$prop_potential_novel <- NA_real_
  projected$ModelInfo_MClust <- NULL
  node_fit <- object$ModelInfo_Hierarchical$node_fit
  node <- projected$ProjectionFit$individual$SOM_Node
  for (column in c("DistanceToMedoid", "AssignmentMargin", "Flag_ProjectionDisagreement")) {
    name <- paste0("Node", column)
    projected$ProjectionFit$individual[[name]] <- node_fit[[column]][node]
    projected$ProbFit$individual[[name]] <- node_fit[[column]][node]
  }
  projected$ProbFit$plots <- .ClusterConfidencePlots(projected$ProbFit$individual,
    "NodeAssignmentMargin", "SOM node margin to the closest other cluster medoid",
    prefix = "Projected", bounded = FALSE)
  projected$ModelInfo <- object$ModelInfo_Hierarchical
  projected$ModelInfo_Hierarchical <- object$ModelInfo_Hierarchical
  projection_specification <- projected$Specification$projection
  projected$Specification <- object$Specification
  projected$Specification$projection <- utils::modifyList(
    projected$Specification$projection, projection_specification)
  projected$Specification$projection$membership_probability <- "unavailable"
  class(projected) <- c("Project_SOM_Hierarchical", class(projected))
  projected
}
