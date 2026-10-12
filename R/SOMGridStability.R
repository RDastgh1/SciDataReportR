#' Screen fixed SOM map designs before selecting cluster count
#'
#' @description In its default `"map"` scope, this function refits only SOM
#'   preprocessing and map training. It selects a representation for later,
#'   fixed-map Mclust comparisons without fitting a mixture or varying `k`.
#'   The legacy joint SOM + Mclust screen remains available with
#'   `screen_scope = "joint"`.
#'
#' @param screen_scope `"map"` (default) for k-informed, Mclust-free SOM screening or
#'   `"joint"` for the earlier grid-by-k Mclust workflow.
#' @param data A data frame supplied to [CreateClusterModel_SOM_MClust()].
#' @param variables Numeric clustering variables.
#' @param k_range Requested profile counts. In `"map"` scope, the largest
#'   requested value defines the automatic candidate-resolution ladder but
#'   Mclust is not fitted and no k is selected.
#' @param models Mclust covariance-model families. In `"map"` scope they are
#'   not fitted, but their largest requested parameter count supplies a
#'   candidate-map support floor.
#' @param grid_candidates Optional square side lengths or a data frame with
#'   `som_xdim` and `som_ydim`. The map-scope default includes approximately
#'   10, 20, and 30 nodes per maximum requested profile, plus 3x3 and 10x10
#'   anchors and a parameter-support floor for the largest requested Mclust
#'   family. User-supplied candidates are retained unchanged.
#' @param topologies SOM topologies to screen.
#' @param node_weighting Node-mixture definitions. Ignored in `"map"` scope.
#' @param min_participants_per_node Target complete participants per SOM node
#'   for automatic candidates. Default `10`; anchors and user grids remain in
#'   the review table even when they do not meet this map-health target.
#' @param screen_resamples Map refits per candidate. The map default is 3;
#'   the explicit joint workflow retains its prior default of 10.
#' @param seed_som Fixed SOM initialization seed used for every refit of a
#'   candidate design.
#' @param stability_seed Seed for participant subsampling and metric sampling.
#' @param stability_cores Optional parallel workers for joint-screen refits;
#'   map-screen refits are intentionally lightweight and serial.
#' @param pair_sample Number of complete-cohort participant pairs used for
#'   orientation-invariant grid-distance preservation.
#' @param neighborhood_sample Number of participants used for local-neighbor
#'   overlap; sampling is deterministic from `stability_seed`.
#' @param neighborhood_k Number of local grid neighbors per sampled participant.
#' @param ... Additional arguments passed to [CreateClusterModel_SOM_MClust()].
#'
#' @return A `SOMGridStabilityScreen` object. Map scope returns a decision
#'   table and per-refit `map_replicates`; joint scope returns the prior
#'   grid-by-k Mclust decision table for [ConfirmSOMGridStability()].
#' @export
ScreenSOMGridStability <- function(data,
    variables = NULL,
    k_range = 2:10,
    models = c(1, 2, 3),
    grid_candidates = NULL,
    topologies = c("hexagonal", "rectangular"),
    node_weighting = c("uniform", "occupancy"),
    min_participants_per_node = 10,
    screen_scope = c("map", "joint"),
    screen_resamples = 3L,
    seed_som = 934521L,
    stability_seed = 934522L,
    stability_cores = NULL,
    pair_sample = 10000L,
    neighborhood_sample = 500L,
    neighborhood_k = 10L,
    ...) {
  resamples_was_missing <- missing(screen_resamples)
  screen_scope <- match.arg(screen_scope)
  if (!is.numeric(min_participants_per_node) ||
      length(min_participants_per_node) != 1L ||
      !is.finite(min_participants_per_node) || min_participants_per_node <= 0) {
    stop("min_participants_per_node must be a single positive finite number.")
  }
  if (screen_scope == "joint") {
    if (resamples_was_missing) screen_resamples <- 10L
    return(.ScreenSOMGridJointStability(
      data = data, variables = variables, k_range = k_range, models = models,
      grid_candidates = grid_candidates, topologies = topologies,
      node_weighting = node_weighting, screen_resamples = screen_resamples,
      seed_som = seed_som, stability_seed = stability_seed,
      stability_cores = stability_cores,
      min_participants_per_node = min_participants_per_node, ...))
  }
  .ScreenSOMMapStability(
    data = data, variables = variables, k_range = k_range,
    models = models,
    grid_candidates = grid_candidates,
    topologies = topologies, screen_resamples = screen_resamples,
    seed_som = seed_som, stability_seed = stability_seed,
    pair_sample = pair_sample, neighborhood_sample = neighborhood_sample,
    neighborhood_k = neighborhood_k,
    min_participants_per_node = min_participants_per_node, ...)
}

#' Screen SOM map designs using participant-level cluster reproducibility
#'
#' @keywords internal
#'
#' @description Fits a separate SOM + Mclust workflow for each requested map
#'   resolution and topology.  For each requested number of classes, the
#'   best-fitting eligible Mclust family is carried into a small, fixed-design
#'   subsampling experiment.  Stability is the adjusted Rand index (ARI) of
#'   the final participant partition after every refit has projected the full
#'   complete reference cohort.  The result is a review table, not an automatic
#'   choice of map resolution or profile count.
#'
#' @param data A data frame supplied to [CreateClusterModel_SOM_MClust()].
#' @param variables Numeric clustering variables.  Defaults to the normal SOM
#'   workflow behaviour.
#' @param k_range Requested profile counts.
#' @param models Mclust covariance-model families.
#' @param grid_candidates Optional square side lengths or a data frame with
#'   `som_xdim` and `som_ydim`. The joint-screen default retains the legacy
#'   occupancy ladder; map-scope screening uses the k-informed 10/20/30
#'   nodes-per-maximum-profile ladder described above.
#' @param topologies SOM topologies to screen.
#' @param node_weighting Node-mixture definitions to screen. `"uniform"`
#'   preserves map-space weighting; `"occupancy"` weights occupied nodes by
#'   their training participant count.
#' @param screen_resamples Number of 90% subsample refits per map design.
#' @param seed_som Fixed SOM initialization seed used for every refit of a
#'   design.
#' @param stability_seed Seed for participant subsampling.
#' @param stability_cores Optional parallel workers for refits.
#' @param ... Additional arguments passed to [CreateClusterModel_SOM_MClust()].
#'
#' @return An object of class `SOMGridStabilityScreen` containing `decision_table`,
#'   `settings`, the source data, and fitting arguments for
#'   `ConfirmSOMGridStability()`.
#' @noRd
.ScreenSOMGridJointStability <- function(data,
    variables = NULL,
    k_range = 2:10,
    models = c(1, 2, 3),
    grid_candidates = NULL,
    topologies = c("hexagonal", "rectangular"),
    node_weighting = c("uniform", "occupancy"),
    screen_resamples = 10L,
    seed_som = 934521L,
    stability_seed = 934522L,
    stability_cores = NULL,
    min_participants_per_node = 10,
    ...) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is.numeric(screen_resamples) || length(screen_resamples) != 1L ||
      screen_resamples < 1L) {
    stop("screen_resamples must be a positive integer.")
  }
  if (!is.character(topologies) || !length(topologies) ||
      any(!topologies %in% c("hexagonal", "rectangular"))) {
    stop("topologies must contain one or both of 'hexagonal' and 'rectangular'.")
  }
  if (!is.character(node_weighting) || !length(node_weighting) ||
      any(!node_weighting %in% c("uniform", "occupancy"))) {
    stop("node_weighting must contain one or both of 'uniform' and 'occupancy'.")
  }

  df_scidr <- .AddClusterRowID(data)
  vars_complete <- if (is.null(variables)) {
    names(df_scidr)[vapply(df_scidr, is.numeric, logical(1))]
  } else {
    variables
  }
  if (!is.character(vars_complete) || !length(vars_complete) ||
      length(setdiff(vars_complete, names(df_scidr)))) {
    stop("variables must identify one or more columns in data.")
  }
  n_complete <- sum(stats::complete.cases(df_scidr[, vars_complete, drop = FALSE]))
  if (n_complete < 3L) stop("At least three complete participants are required.")
  grids <- .SOMGridCandidates(n_complete, grid_candidates,
    k_range = k_range, policy = "k_informed", models = models,
    p = length(vars_complete),
    min_participants_per_node = min_participants_per_node)
  requested_parameter_count <- .SOMGridMaxParameters(models,
    max(as.integer(k_range)), length(vars_complete))

  dots <- list(...)
  prohibited <- intersect(names(dots), c(
    "data", "variables", "method", "k_range", "models", "som_xdim",
    "som_ydim", "som_topo", "seed_som", "stability_resamples",
    "stability_candidates", "stability_seed", "stability_cores", "node_weighting"))
  if (length(prohibited)) {
    stop("Do not supply these arguments through ...: ",
      paste(prohibited, collapse = ", "), ".")
  }

  rows <- list()
  row_index <- 0L
  for (grid_index in seq_len(nrow(grids))) {
    for (topology in unique(topologies)) {
      for (weighting in unique(node_weighting)) {
      fit_args <- c(list(
        data = data, variables = variables, method = "exploratory",
        k_range = k_range, models = models,
        som_xdim = grids$som_xdim[[grid_index]],
        som_ydim = grids$som_ydim[[grid_index]], som_topo = topology,
        node_weighting = weighting,
        seed_som = seed_som, Relabel = FALSE,
        stability_resamples = as.integer(screen_resamples),
        stability_candidates = "best_per_k",
        stability_seed = stability_seed, stability_cores = stability_cores,
        min_participants_per_node = min_participants_per_node), dots)
      fitted <- tryCatch(do.call(CreateClusterModel_SOM_MClust, fit_args),
        error = function(e) e)
      if (inherits(fitted, "error")) {
        row_index <- row_index + 1L
        rows[[row_index]] <- dplyr::tibble(
          .ScreenRow = row_index, Topology = topology,
          NodeWeighting = weighting,
          SOM_XDim = grids$som_xdim[[grid_index]],
          SOM_YDim = grids$som_ydim[[grid_index]],
          NodeCount = grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]],
          ParticipantsPerNode = n_complete / (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]),
          OccupancyCompliant = n_complete / (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) >= min_participants_per_node,
          MaxMixtureParameters = requested_parameter_count,
          MixtureDimension = length(vars_complete),
          NodesPerParameter_MaxRequestedModel =
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) /
              requested_parameter_count,
          Overparameterized_MaxRequestedModel =
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) /
              requested_parameter_count < 1,
          SupportStatus_MaxRequestedModel = dplyr::case_when(
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) / requested_parameter_count < 1 ~ "over-parameterized",
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) / requested_parameter_count < 3 ~ "under-supported",
            TRUE ~ "supported"),
          Model = NA_integer_, Classes = NA_integer_, Eligible = FALSE,
          StabilityARI_Mean = NA_real_, StabilityARI_P05 = NA_real_,
          StabilityJaccard_Mean = NA_real_, StabilitySuccessRate = NA_real_,
          FitError = conditionMessage(fitted))
        next
      }
      fit_table <- fitted$ModelInfo_MClust$fit_table
      occupancy <- fitted$ModelInfo_SOM$SOMFit$occupancy_summary
      lpa_preprocess <- fitted$ModelInfo_MClust$diagnostics$lpa_preprocess
      distance_summary <- fitted$ModelInfo_SOM$SOMFit$distance_summary
      selected <- fit_table %>%
        dplyr::filter(!is.na(.data$StabilityARI_Mean)) %>%
        dplyr::mutate(
          Topology = topology,
          NodeWeighting = weighting,
          SOM_XDim = grids$som_xdim[[grid_index]],
          SOM_YDim = grids$som_ydim[[grid_index]],
          NodeCount = grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]],
          ParticipantsPerNode = n_complete / .data$NodeCount,
          OccupancyCompliant = .data$ParticipantsPerNode >= min_participants_per_node,
          MaxMixtureParameters = .SOMGridMaxParameters(models,
            max(as.integer(k_range)), lpa_preprocess$n_lpa_variables_used),
          MixtureDimension = lpa_preprocess$n_lpa_variables_used,
          NodesPerParameter_MaxRequestedModel = .data$NodeCount /
            .data$MaxMixtureParameters,
          Overparameterized_MaxRequestedModel =
            .data$NodesPerParameter_MaxRequestedModel < 1,
          SupportStatus_MaxRequestedModel = dplyr::case_when(
            .data$Overparameterized_MaxRequestedModel ~ "over-parameterized",
            .data$NodesPerParameter_MaxRequestedModel < 3 ~ "under-supported",
            TRUE ~ "supported"),
          FitSuccessRate = .data$StabilitySuccessRate,
          EmptyNodes = occupancy$n_empty_nodes,
          SingletonNodes = occupancy$n_singleton_nodes,
          NodeUtilization = (occupancy$n_nodes - occupancy$n_empty_nodes) /
            occupancy$n_nodes,
          TotalNodeWeight = lpa_preprocess$total_node_weight,
          EffectiveNodeWeight = lpa_preprocess$effective_node_weight,
          OccupancyMedian = occupancy$median_occupancy,
          OccupancyMaximum = occupancy$max_occupancy,
          SOMDistance_Mean = fitted$ModelInfo_SOM$SOMFit$overall_mean,
          SOMDistance_Median = distance_summary[["median"]],
          SOMDistance_P95 = distance_summary[["p95"]],
          TrainQuantizationError = fitted$ModelInfo_SOM$SOMFit$train_quant_error,
          FitError = NA_character_)
      if (!nrow(selected)) next
      selected$.ScreenRow <- seq.int(row_index + 1L, length.out = nrow(selected))
      row_index <- row_index + nrow(selected)
      rows <- c(rows, split(selected, seq_len(nrow(selected))))
      }
    }
  }
  decision_table <- dplyr::bind_rows(rows)
  if (nrow(decision_table) && "StabilityARI_Mean" %in% names(decision_table)) {
    decision_table <- decision_table %>%
      dplyr::mutate(
        StabilityRank = dplyr::dense_rank(dplyr::desc(.data$StabilityARI_Mean))) %>%
      dplyr::arrange(.data$StabilityRank, dplyr::desc(.data$StabilityARI_P05),
        .data$Topology, .data$NodeWeighting, .data$NodeCount, .data$Classes) %>%
      dplyr::group_by(.data$Topology, .data$NodeWeighting, .data$Classes) %>%
      dplyr::mutate(
        StabilityPlateau = .data$StabilityARI_Mean >= max(.data$StabilityARI_Mean, na.rm = TRUE) - 0.02 &
          .data$StabilityARI_P05 >= max(.data$StabilityARI_P05, na.rm = TRUE) - 0.05) %>%
      dplyr::ungroup()
  }
  out <- list(
    decision_table = decision_table,
    settings = list(
      scoring_cohort = "all_reference_participants",
      perturbation = "90% subsample without replacement",
      refits_requested = as.integer(screen_resamples),
      fixed_grid_policy = "candidate dimensions and topology fixed within refits",
      fixed_initialization_seed_policy = "fixed_to_reference_fit",
      projection_comparison = "refit projected to full reference cohort",
      ranking = "StabilityARI_Mean, then StabilityARI_P05",
      screen_scope = "joint",
      min_participants_per_node = min_participants_per_node,
      stability_plateau = "within 0.02 mean ARI and 0.05 ARI P05 within topology, weighting, and k",
      automatic_selection = FALSE),
    data = data,
    variables = variables,
    fit_arguments = dots,
    k_range = k_range,
    models = models,
    node_weighting = node_weighting,
    seed_som = seed_som,
    stability_seed = stability_seed,
    stability_cores = stability_cores
  )
  class(out) <- "SOMGridStabilityScreen"
  out
}

.SOMMapOnlyNodeClustering <- function(som_codes) {
  list(
    node_cluster = rep(1L, nrow(som_codes)),
    fit_table = data.frame(), ahp_best_row = NULL,
    recommendation = "Map-only SOM screen; Mclust was not fitted.",
    fit_plot = NULL
  )
}

.SOMMapProjection <- function(model, data) {
  df_scidr <- .AddClusterRowID(data)
  z_projection <- if (model$ZScoreType == "PreZScored") {
    df_scidr[, model$ZScoreVars, drop = FALSE]
  } else {
    SciDataReportR::ProjectZScore(
      data = df_scidr, variables = model$vars_used,
      parameters = model$ZScoreObject, ParameterInputType = "ZScoreObj",
      names_prefix = "Z_")$ZScores[, model$ZScoreVars, drop = FALSE]
  }
  complete <- stats::complete.cases(z_projection)
  node <- rep(NA_integer_, nrow(df_scidr))
  distance <- rep(NA_real_, nrow(df_scidr))
  if (any(complete)) {
    mapped <- kohonen::map(model$ModelInfo_SOM$som_model,
      newdata = as.matrix(z_projection[complete, , drop = FALSE]))
    node[complete] <- mapped$unit.classif
    distance[complete] <- sqrt(mapped$distances)
  }
  coordinates <- matrix(NA_real_, nrow(df_scidr), 2L)
  coordinates[complete, ] <- model$ModelInfo_SOM$som_grid$pts[node[complete], , drop = FALSE]
  dplyr::tibble(.row_id = df_scidr$.row_id, SOM_Node = node,
    SOM_Distance = distance, Grid_X = coordinates[, 1], Grid_Y = coordinates[, 2])
}

.SOMMapPairIndices <- function(n, pair_sample, seed) {
  set.seed(seed)
  pair_sample <- min(as.integer(pair_sample), n * (n - 1L) / 2L)
  first <- sample.int(n, pair_sample, replace = TRUE)
  second <- sample.int(n, pair_sample, replace = TRUE)
  same <- first == second
  while (any(same)) {
    second[same] <- sample.int(n, sum(same), replace = TRUE)
    same <- first == second
  }
  cbind(first, second)
}

.SOMMapComparison <- function(reference, refit, pairs, neighborhood_rows,
    neighborhood_k) {
  valid <- !is.na(reference$Grid_X) & !is.na(refit$Grid_X)
  ref_distance <- sqrt((reference$Grid_X[pairs[, 1]] - reference$Grid_X[pairs[, 2]])^2 +
    (reference$Grid_Y[pairs[, 1]] - reference$Grid_Y[pairs[, 2]])^2)
  refit_distance <- sqrt((refit$Grid_X[pairs[, 1]] - refit$Grid_X[pairs[, 2]])^2 +
    (refit$Grid_Y[pairs[, 1]] - refit$Grid_Y[pairs[, 2]])^2)
  pair_valid <- is.finite(ref_distance) & is.finite(refit_distance)
  topology <- if (sum(pair_valid) > 2L) {
    suppressWarnings(stats::cor(ref_distance[pair_valid], refit_distance[pair_valid],
      method = "spearman"))
  } else NA_real_
  sampled <- neighborhood_rows[valid[neighborhood_rows]]
  if (length(sampled) < 2L) return(list(topology = topology, neighborhood = NA_real_))
  k <- min(as.integer(neighborhood_k), length(sampled) - 1L)
  neighbor_overlap <- vapply(sampled, function(row) {
    reference_distance <- (reference$Grid_X[sampled] - reference$Grid_X[[row]])^2 +
      (reference$Grid_Y[sampled] - reference$Grid_Y[[row]])^2
    refit_distance <- (refit$Grid_X[sampled] - refit$Grid_X[[row]])^2 +
      (refit$Grid_Y[sampled] - refit$Grid_Y[[row]])^2
    reference_neighbors <- sampled[order(reference_distance, sampled)][seq_len(k + 1L)]
    refit_neighbors <- sampled[order(refit_distance, sampled)][seq_len(k + 1L)]
    reference_neighbors <- setdiff(reference_neighbors, row)
    refit_neighbors <- setdiff(refit_neighbors, row)
    length(intersect(reference_neighbors, refit_neighbors)) /
      length(union(reference_neighbors, refit_neighbors))
  }, numeric(1))
  list(topology = topology, neighborhood = mean(neighbor_overlap, na.rm = TRUE))
}

.SOMMapStabilityPlateau <- function(decision_table) {
  required <- c("TopologyPreservation_Mean", "ProjectedDistance_Mean",
    "SingletonNodeProportion")
  if (!all(required %in% names(decision_table))) return(decision_table)
  valid <- is.finite(decision_table$TopologyPreservation_Mean) &
    is.finite(decision_table$ProjectedDistance_Mean) &
    is.finite(decision_table$SingletonNodeProportion)
  if (!any(valid)) return(decision_table)
  best_topology <- max(decision_table$TopologyPreservation_Mean[valid])
  best_distance <- min(decision_table$ProjectedDistance_Mean[valid])
  decision_table %>%
    dplyr::mutate(
      MapStabilityPlateau = .data$TopologyPreservation_Mean >= best_topology - .02 &
        .data$ProjectedDistance_Mean <= best_distance * 1.05 &
        .data$SingletonNodeProportion <= .05,
      MapStabilityRank = dplyr::dense_rank(dplyr::desc(.data$TopologyPreservation_Mean))) %>%
    dplyr::arrange(.data$MapStabilityRank,
      dplyr::desc(.data$NeighborhoodOverlap_Mean))
}

.ScreenSOMMapStability <- function(data, variables = NULL, k_range = 2:10,
    models = c(1, 2, 3),
    grid_candidates = NULL, topologies = c("hexagonal", "rectangular"),
    screen_resamples = 3L, seed_som = 934521L, stability_seed = 934522L,
    pair_sample = 10000L, neighborhood_sample = 500L, neighborhood_k = 10L,
    min_participants_per_node = 10,
    ...) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is.numeric(screen_resamples) || length(screen_resamples) != 1L ||
      screen_resamples < 1L) stop("screen_resamples must be a positive integer.")
  df_scidr <- .AddClusterRowID(data)
  vars_complete <- if (is.null(variables)) names(df_scidr)[vapply(df_scidr, is.numeric, logical(1))] else variables
  if (!is.character(vars_complete) || !length(vars_complete) ||
      length(setdiff(vars_complete, names(df_scidr)))) stop("variables must identify columns in data.")
  reference_rows <- which(stats::complete.cases(df_scidr[, vars_complete, drop = FALSE]))
  if (length(reference_rows) < 3L) stop("At least three complete participants are required.")
  grids <- .SOMGridCandidates(length(reference_rows), grid_candidates,
    k_range = k_range, policy = "k_informed", models = models,
    p = length(vars_complete),
    min_participants_per_node = min_participants_per_node)
  dots <- list(...)
  prohibited <- intersect(names(dots), c("data", "variables", "method", "k_range", "models",
    "som_xdim", "som_ydim", "som_topo", "seed_som", "stability_resamples",
    "stability_candidates", "stability_seed", "stability_cores", "node_weighting", ".NodeClusterFn"))
  if (length(prohibited)) stop("Do not supply these arguments through ...: ", paste(prohibited, collapse = ", "), ".")
  pairs <- .SOMMapPairIndices(length(reference_rows), pair_sample, stability_seed)
  set.seed(stability_seed + 1L)
  neighborhood_rows <- sort(sample(reference_rows,
    min(as.integer(neighborhood_sample), length(reference_rows)), replace = FALSE))
  rows <- list()
  replicate_rows <- list()
  row_index <- 0L
  for (grid_index in seq_len(nrow(grids))) {
    for (topology in unique(topologies)) {
      fit_args <- c(list(data = data, variables = variables, method = "exploratory",
        k_range = 2, models = 1, som_xdim = grids$som_xdim[[grid_index]],
        som_ydim = grids$som_ydim[[grid_index]], som_topo = topology,
        seed_som = seed_som, Relabel = FALSE, stability_resamples = 0L,
        min_participants_per_node = min_participants_per_node,
        .NodeClusterFn = .SOMMapOnlyNodeClustering), dots)
      reference_model <- tryCatch(do.call(CreateClusterModel_SOM_MClust, fit_args), error = function(e) e)
      row_index <- row_index + 1L
      if (inherits(reference_model, "error")) {
        rows[[row_index]] <- dplyr::tibble(.ScreenRow = row_index, Topology = topology,
          SOM_XDim = grids$som_xdim[[grid_index]], SOM_YDim = grids$som_ydim[[grid_index]],
          NodeCount = grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]],
          ParticipantsPerNode = length(reference_rows) /
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]),
          OccupancyCompliant = length(reference_rows) /
            (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) >= min_participants_per_node,
          FitError = conditionMessage(reference_model))
        next
      }
      reference_projection <- .SOMMapProjection(reference_model, data)[reference_rows, , drop = FALSE]
      refits <- lapply(seq_len(screen_resamples), function(replicate) {
        set.seed(stability_seed + 100000L + replicate)
        sampled_rows <- sort(sample(reference_rows, max(2L, floor(length(reference_rows) * .90)), replace = FALSE))
        refit_args <- fit_args
        refit_args$data <- df_scidr[sampled_rows, , drop = FALSE]
        refit_model <- tryCatch(do.call(CreateClusterModel_SOM_MClust, refit_args), error = function(e) e)
        if (inherits(refit_model, "error")) return(list(status = "failed", error = conditionMessage(refit_model)))
        projection <- tryCatch(.SOMMapProjection(refit_model, data)[reference_rows, , drop = FALSE], error = function(e) e)
        if (inherits(projection, "error")) return(list(status = "failed", error = conditionMessage(projection)))
        comparison <- .SOMMapComparison(reference_projection, projection, pairs,
          match(neighborhood_rows, reference_rows), neighborhood_k)
        list(status = "success", topology = comparison$topology,
          neighborhood = comparison$neighborhood,
          distance_mean = mean(projection$SOM_Distance, na.rm = TRUE),
          distance_p95 = stats::quantile(projection$SOM_Distance, .95, na.rm = TRUE),
          participants_compared = nrow(projection),
          sampled_participants = length(sampled_rows))
      })
      replicate_rows[[row_index]] <- dplyr::bind_rows(lapply(seq_along(refits), function(replicate) {
        refit <- refits[[replicate]]
        dplyr::tibble(
          .ScreenRow = row_index, Topology = topology,
          SOM_XDim = grids$som_xdim[[grid_index]],
          SOM_YDim = grids$som_ydim[[grid_index]], Replicate = replicate,
          Status = refit$status,
          ParticipantsCompared = if (identical(refit$status, "success")) refit$participants_compared else NA_integer_,
          SampledParticipants = if (identical(refit$status, "success")) refit$sampled_participants else NA_integer_,
          TopologyPreservation = if (identical(refit$status, "success")) refit$topology else NA_real_,
          NeighborhoodOverlap = if (identical(refit$status, "success")) refit$neighborhood else NA_real_,
          ProjectedDistanceMean = if (identical(refit$status, "success")) refit$distance_mean else NA_real_,
          ProjectedDistanceP95 = if (identical(refit$status, "success")) refit$distance_p95 else NA_real_,
          Error = if (identical(refit$status, "success")) NA_character_ else refit$error)
      }))
      successes <- Filter(function(x) identical(x$status, "success"), refits)
      occupancy <- reference_model$ModelInfo_SOM$SOMFit$occupancy_summary
      distances <- reference_model$ModelInfo_SOM$SOMFit$distance_summary
      mixture_dimension <- reference_model$ModelInfo_MClust$diagnostics$lpa_preprocess$n_lpa_variables_used
      max_mixture_parameters <- .SOMGridMaxParameters(models, max(as.integer(k_range)),
        mixture_dimension)
      nodes_per_parameter <- (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) /
        max_mixture_parameters
      rows[[row_index]] <- dplyr::tibble(
        .ScreenRow = row_index, Topology = topology,
        SOM_XDim = grids$som_xdim[[grid_index]], SOM_YDim = grids$som_ydim[[grid_index]],
        NodeCount = grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]],
        ParticipantsPerNode = length(reference_rows) /
          (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]),
        OccupancyCompliant = length(reference_rows) /
          (grids$som_xdim[[grid_index]] * grids$som_ydim[[grid_index]]) >= min_participants_per_node,
        MaxMixtureParameters = max_mixture_parameters,
        MixtureDimension = mixture_dimension,
        NodesPerParameter_MaxRequestedModel = nodes_per_parameter,
        Overparameterized_MaxRequestedModel = nodes_per_parameter < 1,
        SupportStatus_MaxRequestedModel = dplyr::case_when(
          nodes_per_parameter < 1 ~ "over-parameterized",
          nodes_per_parameter < 3 ~ "under-supported",
          TRUE ~ "supported"),
        TopologyPreservation_Mean = if (length(successes)) mean(vapply(successes, `[[`, numeric(1), "topology"), na.rm = TRUE) else NA_real_,
        TopologyPreservation_P05 = if (length(successes)) stats::quantile(vapply(successes, `[[`, numeric(1), "topology"), .05, na.rm = TRUE) else NA_real_,
        NeighborhoodOverlap_Mean = if (length(successes)) mean(vapply(successes, `[[`, numeric(1), "neighborhood"), na.rm = TRUE) else NA_real_,
        ProjectedDistance_Mean = if (length(successes)) mean(vapply(successes, `[[`, numeric(1), "distance_mean"), na.rm = TRUE) else NA_real_,
        ProjectedDistance_P95 = if (length(successes)) mean(vapply(successes, `[[`, numeric(1), "distance_p95"), na.rm = TRUE) else NA_real_,
        EmptyNodes = occupancy$n_empty_nodes, SingletonNodes = occupancy$n_singleton_nodes,
        SingletonNodeProportion = occupancy$n_singleton_nodes / occupancy$n_nodes,
        NodeUtilization = (occupancy$n_nodes - occupancy$n_empty_nodes) / occupancy$n_nodes,
        SOMDistance_Mean = reference_model$ModelInfo_SOM$SOMFit$overall_mean,
        SOMDistance_P95 = distances[["p95"]], FitSuccessRate = length(successes) / screen_resamples,
        FitError = if (length(successes) == screen_resamples) NA_character_ else "One or more map refits failed.")
    }
  }
  decision_table <- dplyr::bind_rows(rows)
  decision_table <- .AddExpectedNodesPerProfile(decision_table, k_range)
  decision_table <- .SOMMapStabilityPlateau(decision_table)
  out <- list(decision_table = decision_table, settings = list(screen_scope = "map",
    resamples = as.integer(screen_resamples), scoring_cohort = "all_reference_participants",
    perturbation = "90% subsample without replacement", fixed_grid_policy = "candidate dimensions and topology fixed within refits",
    initialization_seed_policy = "fixed_to_reference_fit", pair_sample = pair_sample,
    neighborhood_sample = neighborhood_sample, neighborhood_k = neighborhood_k,
    mclust_fitted = FALSE,
    grid_design_policy = "k_informed_10_20_30_nodes_per_maximum_requested_profile_plus_parameter_floor",
    grid_design_k = max(as.integer(k_range)),
    min_nodes_per_parameter = 3,
    min_participants_per_node = min_participants_per_node,
    map_stability_plateau = "within 0.02 topology preservation, 5% projected distance, and 5% singleton nodes",
    automatic_selection = FALSE), map_replicates = dplyr::bind_rows(replicate_rows),
    data = data, variables = variables, fit_arguments = dots,
    seed_som = seed_som, stability_seed = stability_seed)
  class(out) <- "SOMGridStabilityScreen"
  out
}

.SOMGridCandidates <- function(n_complete, grid_candidates = NULL,
    k_range = 2:10, policy = "k_informed",
    models = c(1, 2, 3), p = NULL, min_nodes_per_parameter = 3,
    min_participants_per_node = 10) {
  policy <- match.arg(policy)
  if (is.null(grid_candidates)) {
    if (!is.numeric(k_range) || !length(k_range) || anyNA(k_range) ||
        any(k_range < 2L) || any(k_range != as.integer(k_range))) {
      stop("k_range must contain whole-number profile counts of at least 2.")
    }
    sides <- if (policy == "k_informed") {
      ceiling(sqrt(max(as.integer(k_range)) * c(10L, 20L, 30L)))
    } else {
      ceiling(sqrt(n_complete / c(100, 50, 25)))
    }
    parameter_floor_side <- integer()
    if (policy == "k_informed" && !is.null(p)) {
      max_parameters <- .SOMGridMaxParameters(models, max(as.integer(k_range)), p)
      parameter_floor_side <- ceiling(sqrt(min_nodes_per_parameter * max_parameters))
    }
    ## The fixed anchors are always shown. Generated ladder/floor candidates
    ## above the occupancy ceiling are excluded so an automatic screen cannot
    ## recommend a sparse map.
    anchors <- c(3L, 10L)
    occupancy_ceiling_side <- floor(sqrt(n_complete / min_participants_per_node))
    generated <- c(sides, parameter_floor_side)
    generated <- generated[generated <= occupancy_ceiling_side]
    sides <- sort(unique(as.integer(c(anchors, generated))))
    return(dplyr::tibble(som_xdim = sides, som_ydim = sides))
  }
  if (is.numeric(grid_candidates) && is.null(dim(grid_candidates))) {
    if (any(grid_candidates != as.integer(grid_candidates), na.rm = TRUE)) {
      stop("All candidate grid dimensions must be whole numbers.")
    }
    grids <- dplyr::tibble(som_xdim = as.integer(grid_candidates),
      som_ydim = as.integer(grid_candidates))
  } else if (is.data.frame(grid_candidates) &&
      all(c("som_xdim", "som_ydim") %in% names(grid_candidates))) {
    if (any(grid_candidates$som_xdim != as.integer(grid_candidates$som_xdim), na.rm = TRUE) ||
        any(grid_candidates$som_ydim != as.integer(grid_candidates$som_ydim), na.rm = TRUE)) {
      stop("All candidate grid dimensions must be whole numbers.")
    }
    grids <- dplyr::transmute(grid_candidates,
      som_xdim = as.integer(.data$som_xdim),
      som_ydim = as.integer(.data$som_ydim))
  } else {
    stop("grid_candidates must be numeric square side lengths or a data frame with som_xdim and som_ydim.")
  }
  grids <- dplyr::distinct(grids)
  if (!nrow(grids) || anyNA(grids) || any(grids$som_xdim < 2L) ||
      any(grids$som_ydim < 2L)) {
    stop("All candidate grid dimensions must be integers of at least 2.")
  }
  grids
}

.SOMGridMaxParameters <- function(models, k, p) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 1L) {
    stop("p must be a positive finite mixture dimension.")
  }
  max(vapply(as.integer(models), function(model) {
    .SOMMclustParameterCount(model, as.integer(p), as.integer(k))
  }, integer(1)))
}

.AddExpectedNodesPerProfile <- function(decision_table, k_range) {
  if (!nrow(decision_table) || !"NodeCount" %in% names(decision_table)) {
    return(decision_table)
  }
  k_values <- sort(unique(as.integer(k_range)))
  for (k in k_values) {
    decision_table[[paste0("ExpectedNodesPerProfile_K", k)]] <-
      decision_table$NodeCount / k
  }
  decision_table$GridDesignK <- max(k_values)
  decision_table$ExpectedNodesPerMaximumProfile <-
    decision_table$NodeCount / decision_table$GridDesignK
  decision_table
}

#' Confirm selected SOM grid stability candidates
#'
#' @description Refit exact user-shortlisted rows from
#'   `ScreenSOMGridStability(screen_scope = "joint")` with more 90% subsample
#'   refits. For the default map-only screen, choose and freeze a map first,
#'   then compare k with [CreateClusterModel_SOM_MClust()]. This function does
#'   not select candidates: every requested joint-screen row is returned as a
#'   separate, projectable SOM + Mclust model.
#'
#' @param screen A `SOMGridStabilityScreen` object.
#' @param candidates Values of `.ScreenRow`, or rows from the screen decision
#'   table containing that column.
#' @param confirm_resamples Number of stability refits, default 30.
#' @param distance_percentile Distance percentile threshold for participant
#'   influence flags.
#' @param inclusion_probability Inclusion-probability threshold for influence
#'   flags.
#'
#' @return A `SOMGridStabilityConfirmation` object with a decision table and
#'   named, projectable candidate models.
#' @export
ConfirmSOMGridStability <- function(screen,
    candidates,
    confirm_resamples = 30L,
    distance_percentile = 0.95,
    inclusion_probability = 0.60) {
  if (!inherits(screen, "SOMGridStabilityScreen")) {
    stop("screen must be returned by ScreenSOMGridStability().")
  }
  if (identical(screen$settings$screen_scope, "map")) {
    stop("ConfirmSOMGridStability() confirms screen_scope = 'joint' rows. Choose a map from the map screen, then compare fixed-grid k values with CreateClusterModel_SOM_MClust().")
  }
  if (!is.numeric(confirm_resamples) || length(confirm_resamples) != 1L ||
      confirm_resamples < 1L) stop("confirm_resamples must be a positive integer.")
  candidate_rows <- if (is.data.frame(candidates)) candidates else {
    dplyr::filter(screen$decision_table, .data$.ScreenRow %in% candidates)
  }
  if (!".ScreenRow" %in% names(candidate_rows)) {
    stop("candidates must be .ScreenRow values or decision-table rows containing .ScreenRow.")
  }
  candidate_rows <- dplyr::filter(candidate_rows,
    !is.na(.data$Model), !is.na(.data$Classes))
  if (!nrow(candidate_rows)) stop("No eligible candidate rows were selected.")

  models <- list()
  summaries <- list()
  for (i in seq_len(nrow(candidate_rows))) {
    candidate <- candidate_rows[i, , drop = FALSE]
    fit_args <- c(list(
      data = screen$data, variables = screen$variables, method = "finalize",
      final_k = candidate$Classes[[1]], final_model = candidate$Model[[1]],
      som_xdim = candidate$SOM_XDim[[1]], som_ydim = candidate$SOM_YDim[[1]],
      som_topo = candidate$Topology[[1]], seed_som = screen$seed_som,
      node_weighting = candidate$NodeWeighting[[1]],
      Relabel = FALSE, stability_resamples = as.integer(confirm_resamples),
      stability_candidates = "best_per_k", stability_seed = screen$stability_seed,
      stability_cores = screen$stability_cores), screen$fit_arguments)
    fitted <- do.call(CreateClusterModel_SOM_MClust, fit_args)
    fitted$Influence <- .SOMInfluenceReport(fitted, distance_percentile,
      inclusion_probability)
    name <- paste0("screen_", candidate$.ScreenRow[[1]])
    models[[name]] <- fitted
    confirmed <- fitted$ModelInfo_MClust$fit_table %>%
      dplyr::filter(.data$Model == candidate$Model[[1]],
        .data$Classes == candidate$Classes[[1]]) %>%
      dplyr::mutate(.ScreenRow = candidate$.ScreenRow[[1]],
        Topology = candidate$Topology[[1]], SOM_XDim = candidate$SOM_XDim[[1]],
        SOM_YDim = candidate$SOM_YDim[[1]], NodeCount = candidate$NodeCount[[1]],
        NodeWeighting = candidate$NodeWeighting[[1]],
        InfluenceCandidateCount = sum(fitted$Influence$Flag_InfluenceCandidate))
    summaries[[i]] <- confirmed
  }
  out <- list(
    decision_table = dplyr::bind_rows(summaries),
    models = models,
    settings = list(refits_requested = as.integer(confirm_resamples),
      scoring_cohort = "all_reference_participants",
      projection_comparison = "refit projected to full reference cohort",
      influence_distance_percentile = distance_percentile,
      influence_inclusion_probability = inclusion_probability,
      automatic_selection = FALSE),
    screen = screen
  )
  class(out) <- "SOMGridStabilityConfirmation"
  out
}

.SOMInfluenceReport <- function(model, distance_percentile = 0.95,
    inclusion_probability = 0.60) {
  participant <- model$Stability$participant_inclusion
  if (is.null(participant) || !nrow(participant)) {
    return(dplyr::tibble())
  }
  selected <- model$Specification$selected
  participant <- participant %>%
    dplyr::filter(.data$Model == selected$model, .data$Classes == selected$k) %>%
    dplyr::select(".row_id", "InclusionProbability", "SuccessfulRefits")
  individual <- model$ProbFit$individual %>%
    dplyr::select(dplyr::any_of(c(".row_id", "SOM_Distance",
      "SOMDist_percentile_train", "SOM_Node", "Cluster",
      "Top_Distance_Variable", "Top_Distance_Share")))
  dplyr::left_join(individual, participant, by = ".row_id") %>%
    dplyr::mutate(
      Flag_InfluenceCandidate = .data$SOMDist_percentile_train >= distance_percentile &
        .data$InclusionProbability <= inclusion_probability)
}

#' Refit an explicit participant-exclusion SOM sensitivity model
#'
#' @description Creates a separately labelled sensitivity model after a
#'   reviewer-specified exclusion.  It never changes the source model and
#'   records the excluded row IDs and stated reason in the returned object.
#'
#' @param object A confirmed or ordinary `Pipeline_SOM_MClust` object.
#' @param data The original analysis data.
#' @param exclude_row_ids Explicit `.row_id` values to exclude.
#' @param exclusion_reason Human-readable reason recorded in provenance.
#' @param stability_resamples Sensitivity-model refits, default 30.
#' @param ... Additional arguments passed to [CreateClusterModel_SOM_MClust()].
#'
#' @return A separate `Pipeline_SOM_MClust` sensitivity model.
#' @export
RefitSOMSensitivity <- function(object,
    data,
    exclude_row_ids,
    exclusion_reason,
    stability_resamples = 30L,
    ...) {
  if (!inherits(object, "Pipeline_SOM_MClust")) {
    stop("object must be a Pipeline_SOM_MClust model.")
  }
  if (missing(exclude_row_ids) || !length(exclude_row_ids)) {
    stop("exclude_row_ids must explicitly identify one or more participants.")
  }
  if (missing(exclusion_reason) || !is.character(exclusion_reason) ||
      length(exclusion_reason) != 1L || !nzchar(exclusion_reason)) {
    stop("exclusion_reason must be a non-empty character string.")
  }
  df_scidr <- .AddClusterRowID(data)
  unknown <- setdiff(exclude_row_ids, df_scidr$.row_id)
  if (length(unknown)) stop("Some exclude_row_ids are not present in data.")
  sensitivity_data <- dplyr::filter(df_scidr, !.data$.row_id %in% exclude_row_ids)
  if (!nrow(sensitivity_data)) stop("The requested exclusion removes every participant.")
  grid <- object$ModelInfo_SOM$som_grid
  selected <- object$Specification$selected
  fitted <- CreateClusterModel_SOM_MClust(
    data = sensitivity_data, variables = object$vars_used, method = "finalize",
    final_k = selected$k, final_model = selected$model,
    ClusterVariableName = object$ClusterVariableName,
    ZScoreType = object$ZScoreType, ZScoreObject = object$ZScoreObject,
    ZScoreVars = object$ZScoreVars,
    som_xdim = object$ModelInfo_SOM$som_grid_info$som_xdim_used,
    som_ydim = object$ModelInfo_SOM$som_grid_info$som_ydim_used,
    som_topo = as.character(grid$topo),
    som_neigh = as.character(grid$neighbourhood.fct),
    node_weighting = object$Specification$preprocessing$node_weighting,
    assignment = object$Specification$projection$assignment %||% "hard",
    soft_n_nodes = object$Specification$projection$soft_n_nodes %||% 3L,
    boundary_margin_threshold = object$Specification$projection$boundary_margin_threshold %||% 0.20,
    boundary_entropy_threshold = object$Specification$projection$boundary_entropy_threshold %||% 0.60,
    seed_som = object$Specification$seeds$seed_som,
    seed_lpa = object$Specification$seeds$seed_lpa,
    Relabel = FALSE, stability_resamples = stability_resamples, ...)
  fitted$Sensitivity <- list(
    source_model = "user-specified SOM sensitivity refit",
    excluded_row_ids = exclude_row_ids,
    exclusion_reason = exclusion_reason,
    n_excluded = length(exclude_row_ids),
    original_specification = object$Specification)
  fitted
}
