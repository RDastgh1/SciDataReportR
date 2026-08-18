#' Reusable numeric clustering models
#'
#' @description Fit projectable Mclust or K-means models. These functions keep
#' training scaling parameters and assign new observations without refitting.
#' @keywords internal
.PrepareClusterNumeric <- function(data, variables, Scaling) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (is.null(variables)) variables <- names(data)[vapply(data, is.numeric, logical(1))]
  if (!is.character(variables) || !length(variables)) stop("variables must be a non-empty character vector.")
  if (length(setdiff(variables, names(data)))) stop("data is missing required clustering variables.")
  if (!all(vapply(data[variables], is.numeric, logical(1)))) stop("All clustering variables must be numeric.")
  df_scidr <- data
  if (Scaling == "None" || Scaling == "PreZScored") {
    X <- as.matrix(df_scidr[variables])
    z_obj <- NULL
  } else {
    z_obj <- CreateZScoreObject(
      data = df_scidr, variables = variables, names_prefix = ".scidr_z_",
      RetainLabels = FALSE,
      center = Scaling %in% c("Center and Scale", "Center Only"),
      scale = Scaling %in% c("Center and Scale", "Scale Only")
    )
    X <- as.matrix(z_obj$ZScores)
    colnames(X) <- variables
  }
  list(data = df_scidr, X = X, complete_rows = stats::complete.cases(X),
       variables = variables, Scaling = Scaling, ZScoreObject = z_obj)
}

.ResolveNumericClusterScaling <- function(ZScoreType = NULL, Scaling = NULL) {
  choices <- c("Center and Scale", "Center Only", "Scale Only", "None", "PreZScored")
  if (is.null(ZScoreType) && is.null(Scaling)) return(choices[[1]])
  if (!is.null(ZScoreType) && (!is.character(ZScoreType) || length(ZScoreType) != 1L ||
      !ZScoreType %in% choices)) {
    stop("ZScoreType must be one of: ", paste(choices, collapse = ", "), ".")
  }
  if (!is.null(Scaling) && (!is.character(Scaling) || length(Scaling) != 1L ||
      !Scaling %in% choices)) {
    stop("Scaling must be one of: ", paste(choices, collapse = ", "), ".")
  }
  if (!is.null(ZScoreType) && !is.null(Scaling) && !identical(ZScoreType, Scaling)) {
    stop("ZScoreType and Scaling specify different preprocessing; supply only one or use matching values.")
  }
  if (is.null(ZScoreType)) Scaling else ZScoreType
}

.MclustModelNames <- c(`1` = "EEI", `2` = "VVI", `3` = "EEE")

.ResolveMclustModels <- function(models, argument = "models") {
  if (!is.numeric(models) || !length(models) || anyNA(models) ||
      any(models != as.integer(models)) || any(!as.integer(models) %in% 1:3)) {
    stop(argument, " must contain integer model IDs 1, 2, and/or 3.")
  }
  as.integer(models)
}

.ProjectClusterNumeric <- function(object, new_df) {
  vars <- object$vars_used
  if (!is.data.frame(new_df)) stop("new_df must be a data frame.")
  if (length(setdiff(vars, names(new_df)))) {
    stop("new_df is missing required clustering variable(s): ", paste(setdiff(vars, names(new_df)), collapse = ", "))
  }
  df_scidr <- new_df
  Scaling <- object$Preprocessing$Scaling
  if (Scaling == "None" || Scaling == "PreZScored") {
    X <- as.matrix(df_scidr[vars])
  } else {
    z <- ProjectZScore(data = df_scidr, variables = vars,
                       parameters = object$Preprocessing$ZScoreObject,
                       ParameterInputType = "ZScoreObj", names_prefix = ".scidr_z_")
    X <- as.matrix(z$ZScores)
    colnames(X) <- vars
  }
  list(data = df_scidr, X = X, complete_rows = stats::complete.cases(X),
       variables = vars)
}

.ClusterJaccard <- function(reference, resampled) {
  clusters <- sort(unique(reference[!is.na(reference)]))
  dplyr::bind_rows(lapply(clusters, function(this_cluster) {
    a <- reference == this_cluster
    candidates <- unique(resampled[!is.na(resampled)])
    j <- vapply(candidates, function(candidate) {
      b <- resampled == candidate
      sum(a & b, na.rm = TRUE) / sum(a | b, na.rm = TRUE)
    }, numeric(1))
    dplyr::tibble(Cluster = this_cluster, Jaccard = if (length(j)) max(j) else NA_real_)
  }))
}

.ClusterPartitionMetrics <- function(reference, resampled) {
  keep <- !is.na(reference) & !is.na(resampled)
  reference <- reference[keep]
  resampled <- resampled[keep]
  if (length(reference) < 2L) {
    return(c(VI = NA_real_, NMI = NA_real_, FowlkesMallows = NA_real_))
  }
  tab <- table(reference, resampled)
  n <- sum(tab)
  pij <- tab / n
  pi <- rowSums(pij)
  pj <- colSums(pij)
  nonzero <- pij > 0
  mutual_information <- sum(pij[nonzero] * log(
    pij[nonzero] / (pi[row(pij)[nonzero]] * pj[col(pij)[nonzero]])))
  entropy_reference <- -sum(pi[pi > 0] * log(pi[pi > 0]))
  entropy_resampled <- -sum(pj[pj > 0] * log(pj[pj > 0]))
  vi <- entropy_reference + entropy_resampled - 2 * mutual_information
  nmi <- if ((entropy_reference + entropy_resampled) > 0) {
    2 * mutual_information / (entropy_reference + entropy_resampled)
  } else 1
  pair_counts <- tab * (tab - 1) / 2
  true_pairs <- sum(rowSums(tab) * (rowSums(tab) - 1) / 2)
  predicted_pairs <- sum(colSums(tab) * (colSums(tab) - 1) / 2)
  shared_pairs <- sum(pair_counts)
  precision <- if (predicted_pairs > 0) shared_pairs / predicted_pairs else NA_real_
  recall <- if (true_pairs > 0) shared_pairs / true_pairs else NA_real_
  fm <- if (is.finite(precision) && is.finite(recall)) sqrt(precision * recall) else NA_real_
  c(VI = vi, NMI = nmi, FowlkesMallows = fm)
}

.ClusterMatchedAssignments <- function(reference, resampled, noise_label = NULL) {
  clusters <- sort(unique(reference[!is.na(reference)]))
  if (!is.null(noise_label)) clusters <- setdiff(clusters, noise_label)
  out <- rep(NA_integer_, length(reference))
  for (cluster in clusters) {
    candidates <- unique(resampled[!is.na(resampled)])
    if (!is.null(noise_label)) candidates <- setdiff(candidates, noise_label)
    if (!length(candidates)) next
    scores <- vapply(candidates, function(candidate) {
      a <- reference == cluster
      b <- resampled == candidate
      sum(a & b, na.rm = TRUE) / sum(a | b, na.rm = TRUE)
    }, numeric(1))
    out[reference == cluster] <- candidates[[which.max(scores)]]
  }
  out
}

.ClusterStabilityDiagnostics <- function(reference, assignments, valid_reference,
    noise_label = NULL, coassignment_limit = 2000L) {
  successful <- Filter(Negate(is.null), assignments)
  n_success <- length(successful)
  empty <- list(
    participant_inclusion = dplyr::tibble(), cluster_inclusion = dplyr::tibble(),
    coassignment = list(status = "not_available", reason = "No successful refits.",
      matrix = NULL, row_ids = integer()))
  if (!n_success || !length(valid_reference)) return(empty)
  reference_complete <- reference[valid_reference]
  clusters <- sort(unique(reference_complete[!is.na(reference_complete)]))
  if (!is.null(noise_label)) clusters <- setdiff(clusters, noise_label)
  inclusion_count <- integer(length(valid_reference))
  inclusion_denominator <- integer(length(valid_reference))
  for (result in successful) {
    result_complete <- result[valid_reference]
    matched <- .ClusterMatchedAssignments(reference_complete, result_complete, noise_label)
    eligible <- !is.na(matched) & !is.na(result_complete)
    inclusion_denominator[eligible] <- inclusion_denominator[eligible] + 1L
    inclusion_count[eligible] <- inclusion_count[eligible] +
      as.integer(result_complete[eligible] == matched[eligible])
  }
  participant_inclusion <- dplyr::tibble(
    RowIndex = valid_reference,
    Cluster = reference_complete,
    SuccessfulRefits = inclusion_denominator,
    InclusionProbability = ifelse(inclusion_denominator > 0,
      inclusion_count / inclusion_denominator, NA_real_))
  if (!is.null(noise_label)) participant_inclusion <- dplyr::filter(
    participant_inclusion, .data$Cluster != noise_label)
  cluster_inclusion <- if (nrow(participant_inclusion)) participant_inclusion %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      MeanInclusion = mean(.data$InclusionProbability, na.rm = TRUE),
      P05Inclusion = as.numeric(stats::quantile(.data$InclusionProbability, .05, na.rm = TRUE)),
      MinInclusion = min(.data$InclusionProbability, na.rm = TRUE),
      .groups = "drop") else dplyr::tibble()
  if (length(valid_reference) > coassignment_limit) {
    return(list(participant_inclusion = participant_inclusion,
      cluster_inclusion = cluster_inclusion,
      coassignment = list(status = "skipped", reason = paste0(
        "Complete-case n exceeds coassignment_limit (", coassignment_limit, ")."),
        matrix = NULL, row_ids = valid_reference)))
  }
  coassignment <- matrix(0, nrow = length(valid_reference), ncol = length(valid_reference))
  denominator <- matrix(0, nrow = length(valid_reference), ncol = length(valid_reference))
  for (result in successful) {
    labels <- result[valid_reference]
    eligible <- !is.na(labels)
    if (!is.null(noise_label)) eligible <- eligible & labels != noise_label
    if (!any(eligible)) next
    denominator[eligible, eligible] <- denominator[eligible, eligible] + 1
    coassignment[eligible, eligible] <- coassignment[eligible, eligible] +
      outer(labels[eligible], labels[eligible], FUN = "==")
  }
  coassignment <- ifelse(denominator > 0, coassignment / denominator, NA_real_)
  diag(coassignment) <- ifelse(diag(denominator) > 0, 1, NA_real_)
  list(participant_inclusion = participant_inclusion,
    cluster_inclusion = cluster_inclusion,
    coassignment = list(status = "available", reason = NA_character_,
      matrix = coassignment, row_ids = valid_reference))
}

.ClusterMethod <- function(method) {
  if (length(method) > 1L) method <- method[[1]]
  method <- match.arg(method, c("exploratory", "explore", "finalize"))
  if (identical(method, "explore")) "exploratory" else method
}

.ValidateClusterStability <- function(stability_resamples, stability_seed,
    stability_progress) {
  if (!is.numeric(stability_resamples) || length(stability_resamples) != 1L ||
      is.na(stability_resamples) || stability_resamples < 0 ||
      stability_resamples != as.integer(stability_resamples)) {
    stop("stability_resamples must be a single non-negative integer.")
  }
  if (!is.numeric(stability_seed) || length(stability_seed) != 1L ||
      is.na(stability_seed) || !is.finite(stability_seed)) {
    stop("stability_seed must be a single finite number.")
  }
  if (!is.logical(stability_progress) || length(stability_progress) != 1L ||
      is.na(stability_progress)) {
    stop("stability_progress must be TRUE or FALSE.")
  }
  as.integer(stability_resamples)
}

# Record enough context to reproduce a fitted phenotype model without making
# callers reverse-engineer it from nested engine objects.
.ClusterSpecification <- function(method, variables, seeds = list(),
    preprocessing = list(), candidate_grid = NULL, selected = list(),
    complete_rows = NULL, support = list()) {
  versions <- vapply(c("SciDataReportR", "mclust", "cluster", "dbscan", "poLCA",
                       "kohonen", "aweSOM"), function(package) {
    if (requireNamespace(package, quietly = TRUE))
      as.character(utils::packageVersion(package)) else NA_character_
  }, character(1))
  list(
    method = method, variables = variables, seeds = seeds,
    rng_kind = RNGkind(), preprocessing = preprocessing,
    candidate_grid = candidate_grid, selected = selected,
    fitted_at = Sys.time(), package_versions = versions,
    n_rows = if (is.null(complete_rows)) NA_integer_ else length(complete_rows),
    n_complete = if (is.null(complete_rows)) NA_integer_ else sum(complete_rows),
    support = support)
}

.ClusterProjectionContract <- function(individual, summary, by_cluster,
    out_of_support = dplyr::tibble(), policy = list(), plots = list()) {
  list(individual = individual, summary = summary, by_cluster = by_cluster,
       out_of_support = out_of_support, policy = policy, plots = plots,
       # Backward-compatible names retained until the next major release.
       fit_class = individual$Projection_Fit_Class)
}

.CombineClusterStabilities <- function(stabilities, resamples, seed) {
  if (!length(stabilities)) return(NULL)
  list(
    settings = list(resamples = as.integer(resamples), seed = seed,
      refit_scope = "full_pipeline"),
    replicates = dplyr::bind_rows(lapply(stabilities, `[[`, "replicates")),
    cluster_recovery = dplyr::bind_rows(lapply(stabilities, `[[`, "cluster_recovery")),
    summary = dplyr::bind_rows(lapply(stabilities, `[[`, "summary")),
    failures = dplyr::bind_rows(lapply(stabilities, `[[`, "failures")),
    participant_inclusion = dplyr::bind_rows(lapply(stabilities, `[[`, "participant_inclusion")),
    cluster_inclusion = dplyr::bind_rows(lapply(stabilities, `[[`, "cluster_inclusion")),
    coassignment = if (length(stabilities) == 1L) stabilities[[1]]$coassignment else
      lapply(stabilities, `[[`, "coassignment")
  )
}

.ClusterAHP <- function(fit_table, higher = character(), lower = character(), setting = "candidate") {
  scale_metric <- function(x, high) {
    valid <- is.finite(x)
    out <- rep(NA_real_, length(x))
    if (sum(valid) == 1L) out[valid] <- 1 else if (sum(valid) > 1L) {
      r <- range(x[valid]); out[valid] <- if (diff(r) == 0) 1 else (x[valid] - r[[1]]) / diff(r)
      if (!high) out[valid] <- 1 - out[valid]
    }
    out
  }
  metrics <- intersect(c(higher, lower), names(fit_table))
  for (metric in metrics) fit_table[[paste0(metric, "_scaled")]] <- scale_metric(fit_table[[metric]], metric %in% higher)
  scaled <- grep("_scaled$", names(fit_table), value = TRUE)
  fit_table$ahp_index <- if (length(scaled)) rowMeans(as.matrix(fit_table[scaled]), na.rm = TRUE) else NA_real_
  fit_table$ahp_index[is.nan(fit_table$ahp_index)] <- NA_real_
  fit_table$Recommended <- FALSE
  if (any(is.finite(fit_table$ahp_index))) fit_table$Recommended[which.max(fit_table$ahp_index)] <- TRUE
  recommendation <- if (any(fit_table$Recommended)) {
    row <- fit_table[which(fit_table$Recommended)[1], , drop = FALSE]
    identifiers <- intersect(c("Model", "Classes", "MinPts", "Epsilon"), names(row))
    description <- paste(vapply(identifiers, function(name) {
      paste0(name, " = ", format(row[[name]][[1]], trim = TRUE))
    }, character(1)), collapse = ", ")
    paste0("AHP-style review recommends ", setting,
      if (nzchar(description)) paste0(" (", description, ")") else "",
      ". Review this advisory choice alongside the candidate plots.")
  } else "No candidate had sufficient finite metrics for an AHP-style recommendation."
  list(fit_table = fit_table, AHP = list(ahp_index = fit_table$ahp_index,
    ahp_best_row = fit_table[fit_table$Recommended, , drop = FALSE], recommendation = recommendation))
}

# A plain bootstrap can drop a rare category entirely, and a model refit on
# that sample then has no coefficient for it, so projecting the original cohort
# fails and the replicate is wasted. Swap a carrying row back in for any level
# the draw missed; this keeps every replicate usable at a negligible cost to
# the resample.
.PreserveResampleLevels <- function(data, sampled, variables) {
  if (!length(variables)) return(sampled)
  for (variable in variables) {
    values <- data[[variable]]
    if (!is.factor(values) && !is.character(values)) next
    required <- if (is.factor(values)) levels(values) else
      unique(as.character(values[!is.na(values)]))
    present <- unique(as.character(values[sampled]))
    missing_levels <- setdiff(required, present)
    for (level in missing_levels) {
      donors <- which(as.character(values) == level)
      if (!length(donors)) next
      sampled[sample.int(length(sampled), 1L)] <- donors[
        sample.int(length(donors), 1L)]
    }
  }
  sampled
}

.ClusterBootstrapStability <- function(data, reference, fit_project,
    resamples = 0L, seed = 93422L, candidate = list(), noise_label = NULL,
    progress = FALSE, preserve_levels = character(), replacement = TRUE,
    resample_type = "bootstrap", coassignment_limit = 2000L,
    subsample_resamples = 0L, subsample_fraction = .80,
    prediction_resamples = 0L) {
  if (!is.numeric(resamples) || length(resamples) != 1L || is.na(resamples) ||
      resamples < 0 || resamples != as.integer(resamples)) {
    stop("stability_resamples must be a single non-negative integer.")
  }
  resamples <- as.integer(resamples)
  if (resamples < 1L) return(NULL)
  valid_reference <- which(!is.na(reference))
  replicate_rows <- vector("list", resamples)
  recovery_rows <- vector("list", resamples)
  assignments <- vector("list", resamples)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, resamples)
  for (replicate in seq_len(resamples)) {
    if (isTRUE(progress)) message("Stability bootstrap ", replicate, "/", resamples)
    set.seed(seeds[[replicate]])
    sample_n <- if (replacement) length(valid_reference) else max(2L,
      floor(length(valid_reference) * subsample_fraction))
    sampled <- sample(valid_reference, sample_n, replace = replacement)
    sampled <- .PreserveResampleLevels(data, sampled, preserve_levels)
    result <- tryCatch(fit_project(data[sampled, , drop = FALSE], data), error = function(e) e)
    if (inherits(result, "error")) {
      replicate_rows[[replicate]] <- dplyr::tibble(Replicate = replicate,
        ARI = NA_real_, NoiseSensitivity = NA_real_, NoiseSpecificity = NA_real_,
        Status = "failed", Error = conditionMessage(result))
      next
    }
    compare <- valid_reference[!is.na(result[valid_reference])]
    ari <- if (length(compare) > 1L && requireNamespace("mclust", quietly = TRUE))
      mclust::adjustedRandIndex(reference[compare], result[compare]) else NA_real_
    partition_metrics <- .ClusterPartitionMetrics(reference[compare], result[compare])
    noise_sensitivity <- noise_specificity <- NA_real_
    if (!is.null(noise_label) && length(compare)) {
      truth_noise <- reference[compare] == noise_label
      result_noise <- result[compare] == noise_label
      noise_sensitivity <- if (any(truth_noise)) mean(result_noise[truth_noise]) else NA_real_
      noise_specificity <- if (any(!truth_noise)) mean(!result_noise[!truth_noise]) else NA_real_
    }
    recovery <- .ClusterJaccard(reference[compare], result[compare])
    if (nrow(recovery)) recovery_rows[[replicate]] <- dplyr::mutate(recovery, Replicate = replicate)
    assignments[[replicate]] <- result
    replicate_rows[[replicate]] <- dplyr::tibble(Replicate = replicate, ARI = ari,
      VI = partition_metrics[["VI"]], NMI = partition_metrics[["NMI"]],
      FowlkesMallows = partition_metrics[["FowlkesMallows"]],
      NoiseSensitivity = noise_sensitivity, NoiseSpecificity = noise_specificity,
      Status = "success", Error = NA_character_)
  }
  replicates <- dplyr::bind_rows(replicate_rows)
  cluster_recovery <- dplyr::bind_rows(recovery_rows)
  successful <- dplyr::filter(replicates, .data$Status == "success")
  jaccard <- if (nrow(cluster_recovery)) cluster_recovery$Jaccard else NA_real_
  SafeMean <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  SafeQuantile <- function(x, probability) {
    if (any(is.finite(x))) as.numeric(stats::quantile(
      x, probability, na.rm = TRUE)) else NA_real_
  }
  summary <- dplyr::tibble(
    StabilitySuccessRate = mean(replicates$Status == "success"),
    StabilityARI_Mean = SafeMean(successful$ARI),
    StabilityARI_P05 = SafeQuantile(successful$ARI, .05),
    StabilityJaccard_Mean = SafeMean(jaccard),
    StabilityJaccard_Min = if (any(is.finite(jaccard))) min(jaccard, na.rm = TRUE) else NA_real_,
    NoiseSensitivity = SafeMean(successful$NoiseSensitivity),
    NoiseSpecificity = SafeMean(successful$NoiseSpecificity))
  summary$ReproducibilityScore <- rowMeans(summary[c("StabilityARI_Mean", "StabilityJaccard_Mean")], na.rm = TRUE)
  summary$ReproducibilityScore[is.nan(summary$ReproducibilityScore)] <- NA_real_
  diagnostics <- .ClusterStabilityDiagnostics(
    reference, assignments, valid_reference, noise_label, coassignment_limit)
  for (name in names(candidate)) {
    replicates[[name]] <- candidate[[name]]
    if (nrow(cluster_recovery)) cluster_recovery[[name]] <- candidate[[name]]
    summary[[name]] <- candidate[[name]]
  }
  prediction_strength <- NULL
  if (prediction_resamples > 0L && length(valid_reference) >= 4L) {
    prediction_rows <- vector("list", prediction_resamples)
    set.seed(seed + 1000003L)
    prediction_seeds <- sample.int(.Machine$integer.max, prediction_resamples)
    for (replicate in seq_len(prediction_resamples)) {
      set.seed(prediction_seeds[[replicate]])
      train <- sample(valid_reference, floor(length(valid_reference) / 2))
      test <- setdiff(valid_reference, train)
      projected <- tryCatch(fit_project(data[train, , drop = FALSE], data[test, , drop = FALSE]),
        error = function(e) e)
      independent <- tryCatch(fit_project(data[test, , drop = FALSE], data[test, , drop = FALSE]),
        error = function(e) e)
      if (inherits(projected, "error") || inherits(independent, "error")) {
        prediction_rows[[replicate]] <- dplyr::tibble(Replicate = replicate,
          Cluster = NA_integer_, PredictionStrength = NA_real_, EvaluablePairs = 0L,
          Status = "failed")
        next
      }
      predicted_labels <- projected[!is.na(projected)]
      independent_labels <- independent[!is.na(projected)]
      clusters <- unique(predicted_labels)
      if (!is.null(noise_label)) clusters <- setdiff(clusters, noise_label)
      prediction_rows[[replicate]] <- dplyr::bind_rows(lapply(clusters, function(cluster) {
        members <- which(predicted_labels == cluster)
        pairs <- if (length(members) > 1L) utils::combn(members, 2) else matrix(integer(), 2, 0)
        strength <- if (ncol(pairs)) mean(independent_labels[pairs[1, ]] == independent_labels[pairs[2, ]]) else NA_real_
        dplyr::tibble(Replicate = replicate, Cluster = cluster,
          PredictionStrength = strength, EvaluablePairs = ncol(pairs), Status = "success")
      }))
    }
    prediction_strength <- dplyr::bind_rows(prediction_rows)
  }
  subsample <- if (subsample_resamples > 0L) .ClusterBootstrapStability(
    data = data, reference = reference, fit_project = fit_project,
    resamples = subsample_resamples, seed = seed + 2000003L, candidate = candidate,
    noise_label = noise_label, progress = progress, preserve_levels = preserve_levels,
    replacement = FALSE, resample_type = "subsample",
    coassignment_limit = coassignment_limit, subsample_resamples = 0L,
    prediction_resamples = 0L, subsample_fraction = subsample_fraction) else NULL
  list(settings = list(resamples = resamples, seed = seed,
    refit_scope = "full_pipeline", resample_type = resample_type,
    coassignment_limit = coassignment_limit, noise_policy = if (is.null(noise_label))
      "all clusters included" else "noise excluded from inclusion and coassignment",
    subsample_resamples = subsample_resamples, subsample_fraction = subsample_fraction,
    prediction_resamples = prediction_resamples), replicates = replicates,
    cluster_recovery = cluster_recovery, summary = summary,
    failures = dplyr::filter(replicates, .data$Status != "success"),
    participant_inclusion = diagnostics$participant_inclusion,
    cluster_inclusion = diagnostics$cluster_inclusion,
    coassignment = diagnostics$coassignment, subsample = subsample,
    prediction_strength = prediction_strength)
}

.ClusterOutput <- function(prep, cluster, ClusterVariableName, individual, ModelInfo,
    Stability = NULL) {
  cluster_full <- rep(NA_integer_, nrow(prep$data))
  cluster_full[prep$complete_rows] <- as.integer(cluster)
  df_out <- prep$data
  if (ClusterVariableName %in% names(df_out)) message("Overwriting existing column '", ClusterVariableName, "'.")
  df_out[[ClusterVariableName]] <- cluster_full
  individual <- dplyr::mutate(individual, Cluster = cluster_full)
  list(DataWithClusters = df_out, ProbFit = list(individual = individual),
       ModelInfo = ModelInfo, Stability = Stability)
}

# Freeze the training distance reference so projections can be triaged against
# it later, and build the matching figures. This mirrors what the SOM pipeline
# stores in ModelInfo_SOM$SOMFit and $ProjectionReference.
.ClusterFitDiagnostics <- function(individual, distance_var, distance_label,
    high_quantile = 0.95, noise_label = 0L) {
  distances <- individual[[distance_var]]
  usable <- is.finite(distances)
  if (!any(usable)) {
    return(list(distance_reference = list(high_distance_cutoff = NA_real_),
                plots = list()))
  }
  high_cutoff <- as.numeric(
    stats::quantile(distances[usable], high_quantile, na.rm = TRUE))
  cluster_cutoffs <- individual %>%
    dplyr::filter(is.finite(.data[[distance_var]]), !is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data[[distance_var]]),
      median_distance = stats::median(.data[[distance_var]]),
      cluster_high_cutoff = as.numeric(
        stats::quantile(.data[[distance_var]], high_quantile)),
      prop_high_distance = mean(.data[[distance_var]] > high_cutoff),
      .groups = "drop")
  list(
    distance_reference = list(
      mean = mean(distances[usable]),
      median = stats::median(distances[usable]),
      sd = stats::sd(distances[usable]),
      p90 = as.numeric(stats::quantile(distances[usable], 0.90)),
      p95 = as.numeric(stats::quantile(distances[usable], 0.95)),
      p99 = as.numeric(stats::quantile(distances[usable], 0.99)),
      high_distance_quantile = high_quantile,
      high_distance_cutoff = high_cutoff),
    overall_fit_summary = dplyr::tibble(
      n = sum(usable),
      mean_distance = mean(distances[usable]),
      median_distance = stats::median(distances[usable]),
      sd_distance = stats::sd(distances[usable]),
      high_distance_cutoff = high_cutoff,
      prop_high_distance = mean(distances[usable] > high_cutoff)),
    cluster_fit_summary = cluster_cutoffs,
    plots = c(
      .ClusterDistancePlots(
        individual, distance_var, distance_label, cutoff = high_cutoff,
        prefix = "Training", noise_label = noise_label),
      list(cluster_fit_summary_plot = ggplot2::ggplot(
        dplyr::mutate(
          cluster_cutoffs,
          ClusterDisplay = .ClusterFactor(.data$Cluster, noise_label)),
        ggplot2::aes(x = .data$ClusterDisplay, y = .data$prop_high_distance)) +
        ggplot2::geom_col(fill = "#3B5BDB") +
        ggplot2::scale_y_continuous(labels = scales::label_percent()) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = "High-distance training cases by cluster",
          subtitle = paste0(
            "Cases beyond the overall ",
            sprintf("%.0f", 100 * high_quantile), "th percentile of ",
            tolower(distance_label)),
          x = "Cluster", y = "Proportion"))))
}

.ClusterProbFit <- function(individual, confidence_var, confidence_label,
    threshold = NULL, prefix = "Training", noise_label = 0L) {
  list(individual = individual,
       plots = .ClusterConfidencePlots(
         individual, confidence_var, confidence_label, threshold = threshold,
         prefix = prefix, noise_label = noise_label))
}

# Compare projected cases against the frozen training reference and triage them
# the same way ProjectCluster() does.
.ClusterProjectionFit <- function(individual, reference, distance_var,
    distance_label, confidence_var = NULL, confidence_threshold = NULL,
    training_cluster = NULL, noise_label = 0L) {
  high_cutoff <- if (is.null(reference$high_distance_cutoff)) NA_real_ else
    reference$high_distance_cutoff
  fit_class <- .ClusterProjectionFitClass(
    individual[[distance_var]],
    if (is.null(confidence_var)) NULL else individual[[confidence_var]],
    high_cutoff, confidence_threshold)
  usable <- is.finite(individual[[distance_var]])
  plots <- c(
    .ClusterDistancePlots(
      individual, distance_var, distance_label, cutoff = high_cutoff,
      prefix = "Projected", noise_label = noise_label),
    list(projection_fit_class_bar = .ClusterFitClassPlot(fit_class)))
  # Occupancy comparisons are intentionally not returned: differing cohort
  # composition is not, by itself, evidence of poor individual projection.
  overall_fit_summary <- dplyr::tibble(
      n = sum(usable),
      mean_distance = if (any(usable)) mean(individual[[distance_var]][usable]) else NA_real_,
      median_distance = if (any(usable)) stats::median(individual[[distance_var]][usable]) else NA_real_,
      high_distance_cutoff = high_cutoff,
      prop_high_distance = if (any(usable) && is.finite(high_cutoff)) {
        mean(individual[[distance_var]][usable] > high_cutoff)
      } else NA_real_,
      prop_good_fit = mean(fit_class == "Good fit", na.rm = TRUE),
      prop_uncertain = mean(fit_class == "Uncertain membership", na.rm = TRUE),
      prop_poor_fit = mean(fit_class == "Poor fit to training structure", na.rm = TRUE),
      prop_potential_novel = mean(fit_class == "Potential novel phenotype", na.rm = TRUE))
  cluster_fit_summary <- individual %>%
      dplyr::filter(!is.na(.data$Cluster), is.finite(.data[[distance_var]])) %>%
      dplyr::group_by(.data$Cluster) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_distance = mean(.data[[distance_var]]),
        prop_high_distance = if (is.finite(high_cutoff)) {
          mean(.data[[distance_var]] > high_cutoff)
        } else NA_real_,
        .groups = "drop")
  potential_novel_cases <- individual[
      !is.na(fit_class) & fit_class == "Potential novel phenotype", , drop = FALSE]
  individual$Projection_Fit_Class <- fit_class
  out <- .ClusterProjectionContract(
    individual = individual, summary = overall_fit_summary,
    by_cluster = cluster_fit_summary,
    out_of_support = dplyr::tibble(),
    policy = list(
      high_distance_quantile = reference$high_distance_quantile,
      high_distance_cutoff = high_cutoff,
      confidence_threshold = confidence_threshold),
    plots = plots)
  # Keep old names while downstream code migrates.
  out$distance_reference <- reference
  out$overall_fit_summary <- overall_fit_summary
  out$cluster_fit_summary <- cluster_fit_summary
  out$potential_novel_cases <- potential_novel_cases
  out
}

.MclustCandidate <- function(X, k, model) {
  # mclust::Mclust evaluates mclustBIC in its caller's environment. Keeping
  # this binding local avoids attaching an optional package at runtime.
  mclustBIC <- mclust::mclustBIC
  model <- .ResolveMclustModels(model, "model")
  model_name <- unname(.MclustModelNames[as.character(model)])
  fit <- tryCatch(mclust::Mclust(X, G = k, modelNames = model_name, verbose = FALSE), error = function(e) e)
  if (inherits(fit, "error")) return(list(fit = NULL, error = conditionMessage(fit)))
  if (is.null(fit) || is.null(fit$classification) || is.null(fit$G)) {
    return(list(fit = NULL, error = "Mclust returned no usable classification."))
  }
  entropy <- if (!is.null(fit$z)) {
    p <- pmax(fit$z, .Machine$double.eps)
    1 - (-sum(p * log(p)) / (nrow(p) * log(ncol(p))))
  } else NA_real_
  cluster_sizes <- tabulate(fit$classification, nbins = fit$G)
  list(fit = fit, error = NA_character_, row = dplyr::tibble(
    Model = model, ModelName = model_name, Classes = as.integer(fit$G), BIC = fit$bic,
    ICL = if (is.null(fit$icl)) NA_real_ else fit$icl,
    AIC = -2 * fit$loglik + 2 * fit$df, Entropy = entropy,
    MinClusterN = min(cluster_sizes),
    SizeBalance = min(cluster_sizes) / max(cluster_sizes),
    MaxUncertainty = max(fit$uncertainty, na.rm = TRUE)
  ))
}

#' Fit a projectable Gaussian-mixture clustering model
#'
#' @description Best for continuous measures when clinically meaningful groups
#' may differ in means, variance, or covariance and posterior uncertainty is useful.
#' @references Scrucca L, Fop M, Murphy TB, Raftery AE. mclust 5: Clustering,
#' classification and density estimation using Gaussian finite mixture models.
#' *J Stat Softw.* 2016;8(1):1-29.
#'
#' @param data Data frame containing numeric clustering variables.
#' @param variables Variables used for clustering.
#' @param method Either `"exploratory"` or `"finalize"`.
#' @param k_range Candidate cluster counts in exploratory mode.
#' @param models Numeric Mclust model IDs: `1` = EEI, `2` = VVI, `3` = EEE.
#' @param final_k,final_model Finalized cluster count and numeric model ID.
#' @param ZScoreType Frozen numeric preprocessing. `Scaling` is a compatibility alias.
#' @param Scaling Compatibility alias for `ZScoreType`.
#' @param ClusterVariableName Output cluster column name.
#' @param seed Random seed retained for reproducibility.
#' @param stability_resamples Number of bootstrap refits used to estimate
#'   candidate reproducibility. Use `0` to disable stability analysis.
#' @param stability_seed Seed controlling bootstrap resampling.
#' @param stability_progress Whether to print bootstrap progress messages.
#' @return A projectable mixture model. `ModelInfo$fit_table` contains BIC,
#'   ICL, entropy, uncertainty, and bootstrap stability metrics, and
#'   AIC/BIC/ICL are lower-is-better mixture information criteria; entropy is
#'   higher-is-better classification separation and maximum uncertainty is
#'   lower-is-better. Shared bootstrap fields are defined in the clustering
#'   reference vignette.
#'   `ModelInfo$AHP` the advisory recommendation. Figures sit beside what they
#'   describe: `fit_plot` reviews candidates; `ModelInfo$plots` holds `bic`,
#'   `icl`, `entropy`, `centre_heatmap`, `centre_profile`, `separation_map`,
#'   `classification_uncertainty`, and `profiles`;
#'   `ModelInfo$FitDiagnostics$plots` holds the Mahalanobis distance triad;
#'   `ProbFit$plots` holds posterior-confidence figures; and `Stability$plots`
#'   holds cluster-recovery and complementary stability diagnostics.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Numeric <- paste0("Var", 1:12)
#' review <- CreateClusterModel_MClust(
#'   df_Training, vars_Numeric, k_range = 2:4, models = 1,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$bic
#' model <- CreateClusterModel_MClust(
#'   df_Training, vars_Numeric, method = "finalize", final_k = 4,
#'   final_model = 1, stability_resamples = 2
#' )
#' model$ModelInfo$plots$centre_heatmap
#' model$ModelInfo$plots$separation_map
#' model$ProbFit$plots$confidence_density
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @export
CreateClusterModel_MClust <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10,
    models = c(1L, 2L, 3L), final_k = NULL, final_model = NULL,
    ZScoreType = NULL, Scaling = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, stability_resamples = 0L,
    stability_seed = seed + 1L, stability_progress = FALSE) {
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' is required.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  method <- .ClusterMethod(method)
  Scaling <- .ResolveNumericClusterScaling(ZScoreType, Scaling)
  models <- .ResolveMclustModels(models)
  prep <- .PrepareClusterNumeric(data, variables, Scaling)
  if (!any(prep$complete_rows)) stop("No complete rows available for clustering.")
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  if (method == "finalize") {
    if (is.null(final_k) || is.null(final_model)) stop("final_k and final_model are required for method = 'finalize'.")
    final_model <- .ResolveMclustModels(final_model, "final_model")
    k_range <- final_k; models <- final_model
  }
  set.seed(seed)
  candidates <- unlist(lapply(models, function(m) lapply(k_range, function(k) .MclustCandidate(X, k, m))), recursive = FALSE)
  good <- Filter(function(x) !is.null(x$fit), candidates)
  if (!length(good)) stop("No Mclust candidate could be estimated.")
  fit_table <- dplyr::bind_rows(lapply(good, `[[`, "row")) %>% dplyr::arrange(dplyr::desc(.data$BIC))
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(seq_along(good), function(i) {
      candidate <- good[[i]]; k <- candidate$row$Classes[[1]]; model <- candidate$row$Model[[1]]
      reference <- rep(NA_integer_, nrow(prep$data)); reference[prep$complete_rows] <- candidate$fit$classification
      .ClusterBootstrapStability(prep$data, reference, function(boot, original) {
        fitted <- CreateClusterModel_MClust(boot, prep$variables, method = "finalize",
          final_k = k, final_model = model, ZScoreType = Scaling, seed = seed,
          stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i, list(Model = model, Classes = k),
      progress = stability_progress)
    })
    stability_summary <- dplyr::bind_rows(lapply(stabilities, `[[`, "summary"))
    fit_table <- dplyr::left_join(fit_table, stability_summary, by = c("Model", "Classes"))
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, stability_seed)
  }
  review <- .ClusterAHP(fit_table, higher = c("BIC", "ICL", "Entropy", "MinClusterN", "SizeBalance", "ReproducibilityScore"), lower = "MaxUncertainty", setting = "Gaussian-mixture model")
  fit_table <- review$fit_table
  recommended <- review$AHP$ahp_best_row
  best_i <- if (nrow(recommended)) which(vapply(good, function(x)
    identical(as.character(x$row$Model[[1]]), as.character(recommended$Model[[1]])) &&
      x$row$Classes[[1]] == recommended$Classes[[1]], logical(1)))[1] else
    which.max(vapply(good, function(x) x$row$BIC, numeric(1)))
  best <- good[[best_i]]$fit
  probs <- best$z
  individual <- dplyr::tibble(
    PosteriorMax = rep(NA_real_, nrow(prep$data)),
    Uncertainty = rep(NA_real_, nrow(prep$data)),
    MahalanobisDistance = rep(NA_real_, nrow(prep$data))
  )
  individual$PosteriorMax[prep$complete_rows] <- apply(probs, 1, max)
  individual$Uncertainty[prep$complete_rows] <- best$uncertainty
  individual$MahalanobisDistance[prep$complete_rows] <-
    .MclustDistanceToComponent(X, best, best$classification)
  for (i in seq_len(ncol(probs))) {
    individual[[paste0("prob_", i)]] <- NA_real_
    individual[[paste0("prob_", i)]][prep$complete_rows] <- probs[, i]
  }
  review_space <- .ClusterReviewSpace(X)
  centers <- t(best$parameters$mean)
  colnames(centers) <- prep$variables
  variable_labels <- vapply(
    prep$variables, function(v) .ClusterVariableLabel(prep$data, v), character(1))
  ModelInfo <- list(
    mclust_model = best, fit_table = fit_table, AHP = review$AHP,
    ReviewSpace = review_space, centers = centers,
    final_k = best$G, final_model = as.integer(names(.MclustModelNames)[.MclustModelNames == best$modelName]),
    final_model_name = best$modelName)
  base <- .ClusterOutput(prep, best$classification, ClusterVariableName, individual,
    ModelInfo = ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual
  centre_scale <- .ClusterCentreScaleLabel(Scaling)
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, "MahalanobisDistance", "Distance to component centre")
  ModelInfo$plots <- .DropNullPlots(list(
    bic = .ClusterCandidateCurve(fit_table, "BIC", "Bayesian information criterion"),
    icl = .ClusterCandidateCurve(fit_table, "ICL", "Integrated complete-data likelihood"),
    entropy = .ClusterCandidateCurve(fit_table, "Entropy", "Classification entropy"),
    centre_heatmap = PlotClusterCentreHeatmap(
      centers, variable_labels, "Component mean profiles", centre_scale),
    centre_profile = PlotClusterCentreProfile(
      centers, variable_labels, "Component mean profiles", centre_scale),
    profiles = PlotClusterProfiles(
      base$DataWithClusters, ClusterVariableName, prep$variables)))
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(
    individual, "PosteriorMax", "Posterior probability for the assigned class")
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(method = method, vars_used = prep$variables, ClusterVariableName = ClusterVariableName,
                Preprocessing = list(ZScoreType = Scaling, Scaling = Scaling, ZScoreObject = prep$ZScoreObject),
                fit_plot = PlotClusterFitReview(fit_table)), base)
  out$Stability <- Stability
  out$ModelInfo_MClust <- out$ModelInfo
  out$Specification <- .ClusterSpecification(
    "Mclust", prep$variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, dplyr::select(fit_table, dplyr::any_of(c("Model", "Classes"))),
    list(k = best$G, model = ModelInfo$final_model, model_name = best$modelName), prep$complete_rows,
    list(distance_metric = "assigned-component Mahalanobis distance"))
  class(out) <- c("Pipeline_MClust", class(out)); out
}

.ClusterCentreScaleLabel <- function(Scaling) {
  if (identical(Scaling, "Center and Scale") || identical(Scaling, "PreZScored")) {
    "Standard deviations"
  } else if (identical(Scaling, "Center Only")) {
    "Deviation from mean"
  } else "Centre"
}

.DropNullPlots <- function(plots) plots[!vapply(plots, is.null, logical(1))]

# Distance from each case to the centre of its assigned mixture component, in
# that component's own covariance metric. This is the mixture-model analogue of
# the SOM distance-to-best-matching-unit and gives projections a frozen cutoff.
.MclustDistanceToComponent <- function(X, fit, classification) {
  means <- fit$parameters$mean
  variance <- fit$parameters$variance
  vapply(seq_len(nrow(X)), function(i) {
    component <- classification[[i]]
    centre <- means[, component]
    delta <- as.numeric(X[i, ]) - centre
    sigma <- tryCatch({
      if (!is.null(variance$sigma)) variance$sigma[, , component] else
        diag(variance$scale, nrow = length(delta))
    }, error = function(e) NULL)
    if (is.null(sigma)) return(sqrt(sum(delta^2)))
    solved <- tryCatch(solve(sigma, delta), error = function(e) NULL)
    if (is.null(solved)) sqrt(sum(delta^2)) else sqrt(max(0, sum(delta * solved)))
  }, numeric(1))
}

.MclustUncertaintyMap <- function(review_space, X, classification, uncertainty) {
  coordinates <- .ClusterReviewCoordinates(review_space, X)
  if (is.null(coordinates)) return(NULL)
  coordinates$Cluster <- .ClusterFactor(classification, NULL)
  coordinates$Uncertainty <- uncertainty
  axis_labels <- .ClusterReviewAxisLabels(review_space)
  ggplot2::ggplot(
    coordinates,
    ggplot2::aes(
      x = .data$ReviewDimension1, y = .data$ReviewDimension2,
      color = .data$Cluster, size = .data$Uncertainty)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(0.5, 4)) +
    .ClusterColourScale(coordinates$Cluster) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Classification uncertainty",
      subtitle = "Larger points sit closer to a component boundary",
      x = axis_labels[[1]], y = axis_labels[[2]],
      color = "Cluster", size = "Uncertainty")
}

#' Project data onto a Mclust cluster model
#' @param object A `CreateClusterModel_MClust()` result.
#' @param new_df New observations.
#' @param ClusterVariableName Optional output cluster column name.
#' @return A full-length projection object. `ProbFit$individual` holds posterior
#'   membership, uncertainty, and distance to the assigned component;
#'   `ProjectionFit` triages every projected case against the frozen training
#'   distance reference and holds the projected distance triad, the
#'   fit-class summaries and method-specific projection diagnostics;
#'   `ProbFit$plots` holds posterior-confidence figures.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Numeric <- paste0("Var", 1:12)
#' model <- CreateClusterModel_MClust(
#'   df_Training, vars_Numeric, method = "finalize", final_k = 4,
#'   final_model = 1, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$overall_fit_summary
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' projected$ProjectionFit$plots$training_vs_projected_cluster_occupancy
#' projected$ProjectionFit$plots$separation_map
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_MClust <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  prep <- .ProjectClusterNumeric(object, new_df)
  ind <- dplyr::tibble(
    PosteriorMax = rep(NA_real_, nrow(prep$data)),
    Uncertainty = rep(NA_real_, nrow(prep$data)),
    MahalanobisDistance = rep(NA_real_, nrow(prep$data)))
  cl <- rep(NA_integer_, nrow(prep$data))
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  if (any(prep$complete_rows)) {
    prediction <- stats::predict(object$ModelInfo$mclust_model, newdata = X)
    cl[prep$complete_rows] <- prediction$classification
    ind$PosteriorMax[prep$complete_rows] <- apply(prediction$z, 1, max)
    ind$Uncertainty[prep$complete_rows] <- if (is.null(prediction$uncertainty)) 1 - apply(prediction$z, 1, max) else prediction$uncertainty
    ind$MahalanobisDistance[prep$complete_rows] <- .MclustDistanceToComponent(
      X, object$ModelInfo$mclust_model, prediction$classification)
    for (i in seq_len(ncol(prediction$z))) { ind[[paste0("prob_", i)]] <- NA_real_; ind[[paste0("prob_", i)]][prep$complete_rows] <- prediction$z[, i] }
  }
  base <- .ClusterOutput(
    prep, cl[prep$complete_rows], ClusterVariableName, ind, object$ModelInfo)
  individual <- base$ProbFit$individual
  ProjectionFit <- .ClusterProjectionFit(
    individual, object$ModelInfo$FitDiagnostics$distance_reference,
    "MahalanobisDistance", "Distance to component centre",
    confidence_var = "PosteriorMax", confidence_threshold = 0.8,
    training_cluster = object$ProbFit$individual$Cluster)
  ProjectionFit$plots$separation_map <- .ClusterSeparationMap(
    object$ModelInfo$ReviewSpace, X, cl[prep$complete_rows],
    object$ModelInfo$centers, prefix = "Projected")
  out <- c(list(vars_used = object$vars_used, ClusterVariableName = ClusterVariableName,
                Preprocessing = object$Preprocessing), base)
  out$ProbFit <- .ClusterProbFit(
    individual, "PosteriorMax", "Posterior probability for the assigned class",
    threshold = 0.8, prefix = "Projected")
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- ProjectionFit$fit_class
  out$ProbFit$individual$Projection_Fit_Class <- ProjectionFit$fit_class
  class(out) <- c("Project_Mclust", class(out)); out
}

#' Fit a projectable K-means clustering model
#'
#' @description Best for approximately spherical, similarly sized groups in a
#' numeric feature space; use scaling unless all variables are commensurate.
#' @references MacQueen JB. Some methods for classification and analysis of
#' multivariate observations. In: *Proceedings of the Fifth Berkeley Symposium*; 1967.
#' @inheritParams CreateClusterModel_MClust
#' @param final_k Number of clusters for a finalized K-means solution.
#' @param nstart Number of random K-means starts.
#' @return A fitted model with `ModelInfo$fit_table` (WSS, between-cluster
#'   sum of squares, silhouette, Calinski-Harabasz index, and minimum cluster
#'   size) and bootstrap stability. WSS is lower-is-better; between-cluster
#'   sum of squares, silhouette, and Calinski-Harabasz are higher-is-better.
#'   Figures sit beside what they describe:
#'   `fit_plot` reviews candidates; `ModelInfo$plots` holds `elbow`,
#'   `silhouette` (the per-participant silhouette profile of the selected
#'   solution), `silhouette_by_k`, `calinski_harabasz`, `centre_heatmap`,
#'   `centre_profile`, `separation_map`, and `profiles`;
#'   `ModelInfo$FitDiagnostics$plots` holds the distance-to-centroid triad;
#'   `ProbFit$plots` holds assignment-margin figures; and `Stability$plots`
#'   holds cluster-recovery and complementary stability diagnostics.
#'   `ModelInfo$ReviewSpace` freezes the two-dimensional review space shared by
#'   the training and projection maps; it is diagnostic only and does not
#'   affect clustering.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Numeric <- paste0("Var", 1:12)
#' review <- CreateClusterModel_KMeans(
#'   df_Training, vars_Numeric, k_range = 2:5, nstart = 10,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$elbow
#' review$ModelInfo$plots$silhouette_by_k
#' model <- CreateClusterModel_KMeans(
#'   df_Training, vars_Numeric, method = "finalize", final_k = 4,
#'   nstart = 10, stability_resamples = 2
#' )
#' model$ModelInfo$plots$silhouette
#' model$ModelInfo$plots$centre_heatmap
#' model$ModelInfo$plots$separation_map
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @export
CreateClusterModel_KMeans <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10, final_k = NULL,
    ZScoreType = NULL, Scaling = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, nstart = 50L,
    stability_resamples = 0L, stability_seed = seed + 1L,
    stability_progress = FALSE) {
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  method <- .ClusterMethod(method)
  Scaling <- .ResolveNumericClusterScaling(ZScoreType, Scaling)
  prep <- .PrepareClusterNumeric(data, variables, Scaling)
  X <- prep$X[prep$complete_rows, , drop = FALSE]; if (!nrow(X)) stop("No complete rows available for clustering.")
  ks <- if (method == "finalize") { if (is.null(final_k)) stop("final_k is required for method = 'finalize'."); final_k } else k_range
  fits <- lapply(ks, function(k) { set.seed(seed + k); stats::kmeans(X, centers = k, nstart = nstart) })
  rows <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]; sil <- NA_real_
    if (requireNamespace("cluster", quietly = TRUE) && nrow(fit$centers) > 1 && nrow(X) > 2) sil <- mean(cluster::silhouette(fit$cluster, dist(X))[, "sil_width"])
    ch <- if (fit$tot.withinss > 0 && length(fit$size) > 1) {
      (fit$betweenss / (length(fit$size) - 1)) / (fit$tot.withinss / (nrow(X) - length(fit$size)))
    } else NA_real_
    dplyr::tibble(Classes = ks[[i]], WSS = fit$tot.withinss,
      BetweenSS = fit$betweenss, CalinskiHarabasz = ch,
      Silhouette = sil, MinClusterN = min(fit$size),
      SizeBalance = min(fit$size) / max(fit$size))
  }))
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(seq_along(ks), function(i) {
      k <- ks[[i]]; reference <- rep(NA_integer_, nrow(prep$data)); reference[prep$complete_rows] <- fits[[i]]$cluster
      .ClusterBootstrapStability(prep$data, reference, function(boot, original) {
        fitted <- CreateClusterModel_KMeans(boot, prep$variables, method = "finalize",
          final_k = k, ZScoreType = Scaling, seed = seed, nstart = nstart,
          stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + k, list(Classes = k), progress = stability_progress)
    })
    stability_summary <- dplyr::bind_rows(lapply(stabilities, `[[`, "summary"))
    rows <- dplyr::left_join(rows, stability_summary, by = "Classes")
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, stability_seed)
  }
  review <- .ClusterAHP(rows, higher = c("Silhouette", "CalinskiHarabasz", "MinClusterN", "SizeBalance", "ReproducibilityScore"), lower = "WSS", setting = "K-means k")
  rows <- review$fit_table
  best_i <- if (any(rows$Recommended)) which(rows$Recommended)[1] else if (all(is.na(rows$Silhouette))) which.min(rows$WSS) else which.max(rows$Silhouette)
  best <- fits[[best_i]]
  ind <- dplyr::tibble(DistanceToCentroid = rep(NA_real_, nrow(prep$data)), AssignmentMargin = rep(NA_real_, nrow(prep$data)))
  d <- as.matrix(stats::dist(rbind(best$centers, X)))[seq_len(nrow(best$centers)), nrow(best$centers) + seq_len(nrow(X)), drop = FALSE]
  dd <- t(d); ind$DistanceToCentroid[prep$complete_rows] <- apply(dd, 1, min)
  ind$AssignmentMargin[prep$complete_rows] <- apply(dd, 1, function(x) diff(sort(x))[1])
  silhouette_object <- NULL
  if (requireNamespace("cluster", quietly = TRUE) &&
      nrow(best$centers) > 1L && nrow(X) > 2L) {
    silhouette_object <- tryCatch(
      cluster::silhouette(best$cluster, stats::dist(X)), error = function(e) NULL)
  }
  review_space <- .ClusterReviewSpace(X)
  variable_labels <- vapply(
    prep$variables, function(v) .ClusterVariableLabel(prep$data, v), character(1))
  ModelInfo <- list(
    kmeans_model = best, fit_table = rows, AHP = review$AHP,
    ReviewSpace = review_space, centers = best$centers,
    silhouette = silhouette_object, final_k = nrow(best$centers))
  base <- .ClusterOutput(prep, best$cluster, ClusterVariableName, ind,
    ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual
  centre_scale <- .ClusterCentreScaleLabel(Scaling)
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, "DistanceToCentroid", "Distance to assigned centroid")
  ModelInfo$plots <- .DropNullPlots(list(
    elbow = .ClusterCandidateCurve(
      rows, "WSS", "Within-cluster sum of squares by k",
      "Within-cluster sum of squares"),
    silhouette = PlotClusterSilhouette(
      silhouette_object,
      paste0("Silhouette profile (k = ", nrow(best$centers), ")")),
    silhouette_by_k = .ClusterCandidateCurve(
      rows, "Silhouette", "Average silhouette width by k",
      "Average silhouette width"),
    calinski_harabasz = .ClusterCandidateCurve(
      rows, "CalinskiHarabasz", "Calinski-Harabasz index by k"),
    centre_heatmap = PlotClusterCentreHeatmap(
      best$centers, variable_labels, "Cluster centroid profiles", centre_scale),
    centre_profile = PlotClusterCentreProfile(
      best$centers, variable_labels, "Cluster centroid profiles", centre_scale),
    profiles = PlotClusterProfiles(
      base$DataWithClusters, ClusterVariableName, prep$variables)))
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(
    individual, "AssignmentMargin",
    "Distance margin to the second-nearest centroid")
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(method = method, vars_used = prep$variables, ClusterVariableName = ClusterVariableName,
                Preprocessing = list(ZScoreType = Scaling, Scaling = Scaling, ZScoreObject = prep$ZScoreObject),
                fit_plot = PlotClusterFitReview(rows)), base)
  out$Stability <- Stability
  out$ModelInfo_KMeans <- out$ModelInfo
  out$Specification <- .ClusterSpecification(
    "KMeans", prep$variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, dplyr::select(rows, dplyr::any_of("Classes")),
    list(k = nrow(best$centers)), prep$complete_rows,
    list(distance_metric = "Euclidean distance to frozen centroid"))
  class(out) <- c("Pipeline_KMeans", class(out)); out
}

#' Project data onto a K-means cluster model
#' @inheritParams ProjectCluster
#' @return A full-length projection object with nearest-centroid assignments,
#'   distances, and assignment margins. `ProjectionFit` triages every projected
#'   case against the frozen training distance reference and holds the projected
#'   distance triad, the frozen-centroid `separation_map`, the
#'   fit-class summaries and method-specific projection diagnostics.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize",
#'   final_k = 4, nstart = 10, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$overall_fit_summary
#' projected$ProjectionFit$plots$separation_map
#' projected$ProjectionFit$plots$training_vs_projected_cluster_occupancy
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_KMeans <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  prep <- .ProjectClusterNumeric(object, new_df); cl <- rep(NA_integer_, nrow(prep$data)); ind <- dplyr::tibble(DistanceToCentroid = rep(NA_real_, nrow(prep$data)), AssignmentMargin = rep(NA_real_, nrow(prep$data)))
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  centers <- object$ModelInfo$kmeans_model$centers
  if (any(prep$complete_rows)) {
    ds <- sapply(seq_len(nrow(centers)), function(i) rowSums((sweep(X, 2, centers[i, ], "-"))^2))
    if (is.null(dim(ds))) ds <- matrix(ds, ncol = 1)
    cl[prep$complete_rows] <- max.col(-ds)
    ind$DistanceToCentroid[prep$complete_rows] <- sqrt(apply(ds, 1, min))
    ind$AssignmentMargin[prep$complete_rows] <- apply(
      ds, 1, function(x) diff(sqrt(sort(x)))[1])
  }
  base <- .ClusterOutput(
    prep, cl[prep$complete_rows], ClusterVariableName, ind, object$ModelInfo)
  individual <- base$ProbFit$individual
  ProjectionFit <- .ClusterProjectionFit(
    individual, object$ModelInfo$FitDiagnostics$distance_reference,
    "DistanceToCentroid", "Distance to assigned centroid",
    training_cluster = object$ProbFit$individual$Cluster)
  ProjectionFit$plots$separation_map <- .ClusterSeparationMap(
    object$ModelInfo$ReviewSpace, X, cl[prep$complete_rows], centers,
    prefix = "Projected")
  out <- c(list(vars_used = object$vars_used, ClusterVariableName = ClusterVariableName, Preprocessing = object$Preprocessing), base)
  out$ProbFit <- .ClusterProbFit(
    individual, "AssignmentMargin",
    "Distance margin to the second-nearest centroid", prefix = "Projected")
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- ProjectionFit$fit_class
  out$ProbFit$individual$Projection_Fit_Class <- ProjectionFit$fit_class
  class(out) <- c("Project_KMeans", class(out)); out
}

.CreatePCACluster <- function(data, variables, algorithm, method, k_range,
    final_k, final_model, Scaling, ClusterVariableName, seed, pca_variance_threshold,
    nstart = NULL, models = c(1L, 2L, 3L),
    stability_resamples = 0L, stability_seed = seed + 1L,
    stability_progress = FALSE) {
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  prep <- .PrepareClusterNumeric(data, variables, Scaling)
  if (!any(prep$complete_rows)) stop("No complete rows available for PCA clustering.")
  df_PCA <- as.data.frame(prep$X[prep$complete_rows, , drop = FALSE])
  pca <- CreatePCAObject(
    data = df_PCA, VarsToReduce = colnames(df_PCA), backend = "prcomp",
    center = FALSE, scale = FALSE, minThresh = pca_variance_threshold,
    numComponents = NULL, MissingDataStrategy = "complete_cases", Relabel = FALSE)
  df_Scores <- as.data.frame(pca$Scores)
  score_vars <- names(df_Scores)

  FitScoreModel <- function(fit_method, candidate_k = final_k,
      candidate_model = final_model) {
    if (algorithm == "mclust") {
      CreateClusterModel_MClust(
        df_Scores, score_vars, method = fit_method, k_range = k_range,
        models = models, final_k = candidate_k,
        final_model = candidate_model, ZScoreType = "None",
        ClusterVariableName = ClusterVariableName, seed = seed, stability_resamples = 0L)
    } else {
      CreateClusterModel_KMeans(
        df_Scores, score_vars, method = fit_method, k_range = k_range,
        final_k = candidate_k, ZScoreType = "None", ClusterVariableName = ClusterVariableName,
        seed = seed, nstart = nstart, stability_resamples = 0L)
    }
  }

  review_model <- FitScoreModel(method)
  fit_table <- review_model$ModelInfo$fit_table
  Stability <- NULL
  if (stability_resamples > 0L) {
    candidate_rows <- if (algorithm == "mclust") {
      dplyr::distinct(fit_table, .data$Model, .data$Classes)
    } else dplyr::distinct(fit_table, .data$Classes)
    stabilities <- lapply(seq_len(nrow(candidate_rows)), function(i) {
      candidate_k <- candidate_rows$Classes[[i]]
      candidate_model <- if (algorithm == "mclust") candidate_rows$Model[[i]] else NULL
      reference_fit <- .CreatePCACluster(
        data, prep$variables, algorithm, "finalize", candidate_k,
        candidate_k, candidate_model, Scaling, ClusterVariableName, seed,
        pca_variance_threshold, nstart, models, 0L, stability_seed, FALSE)
      reference <- reference_fit$ProbFit$individual$Cluster
      .ClusterBootstrapStability(data, reference, function(boot, original) {
        fitted <- .CreatePCACluster(
          boot, prep$variables, algorithm, "finalize", candidate_k,
          candidate_k, candidate_model, Scaling, ClusterVariableName, seed,
          pca_variance_threshold, nstart, models, 0L,
          stability_seed, FALSE)
        .ProjectPCACluster(fitted, original, algorithm, ClusterVariableName)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i,
      candidate = if (algorithm == "mclust") {
        list(Model = candidate_model, Classes = candidate_k)
      } else list(Classes = candidate_k), progress = stability_progress)
    })
    Stability <- .CombineClusterStabilities(
      stabilities, stability_resamples, stability_seed)
    join_by <- if (algorithm == "mclust") c("Model", "Classes") else "Classes"
    fit_table <- dplyr::left_join(fit_table, Stability$summary, by = join_by)
  }

  review <- if (algorithm == "mclust") {
    .ClusterAHP(fit_table,
      higher = c("BIC", "ICL", "Entropy", "MinClusterN", "SizeBalance", "ReproducibilityScore"),
      lower = "MaxUncertainty", setting = "PCA plus Gaussian-mixture model")
  } else {
    .ClusterAHP(fit_table,
      higher = c("Silhouette", "CalinskiHarabasz", "MinClusterN", "SizeBalance", "ReproducibilityScore"),
      lower = "WSS", setting = "PCA plus K-means k")
  }
  recommended <- review$AHP$ahp_best_row
  selected_k <- recommended$Classes[[1]]
  selected_model <- if (algorithm == "mclust") recommended$Model[[1]] else NULL
  selected <- FitScoreModel("finalize", selected_k, selected_model)
  selected$method <- method
  selected$ModelInfo$fit_table <- review$fit_table
  selected$ModelInfo$AHP <- review$AHP
  selected$Stability <- Stability

  reduced_individual <- selected$ProbFit$individual
  indicator_names <- setdiff(names(reduced_individual), "Cluster")
  individual <- as.data.frame(lapply(reduced_individual[indicator_names], function(x) {
    if (is.logical(x)) rep(NA, nrow(prep$data)) else if (is.character(x)) {
      rep(NA_character_, nrow(prep$data))
    } else rep(NA_real_, nrow(prep$data))
  }))
  individual[prep$complete_rows, indicator_names] <- reduced_individual[indicator_names]
  base <- .ClusterOutput(
    prep, reduced_individual$Cluster, ClusterVariableName, individual,
    selected$ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual

  # The inner model already built its structure figures in PCA score space,
  # which is where this pipeline actually clusters. Add the reduction layer's
  # own figures and restate the profiles in the original measurement scale,
  # because that is what a reader has to interpret clinically.
  ModelInfo <- selected$ModelInfo
  distance_var <- if (algorithm == "mclust") "MahalanobisDistance" else
    "DistanceToCentroid"
  distance_label <- if (algorithm == "mclust") "Distance to component centre" else
    "Distance to assigned centroid"
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, distance_var, distance_label)
  score_plots <- selected$ModelInfo$plots[
    !names(selected$ModelInfo$plots) %in% "profiles"]
  ModelInfo$plots <- .DropNullPlots(c(
    list(scree = pca$p_scree, loadings = pca$Lollipop),
    score_plots,
    list(
      profiles = PlotClusterProfiles(
        base$DataWithClusters, ClusterVariableName, prep$variables))))
  ModelInfo$fit_table <- review$fit_table
  ModelInfo$AHP <- review$AHP
  base$ModelInfo <- ModelInfo
  base$ProbFit <- selected$ProbFit
  base$ProbFit$individual <- individual
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(
    method = method, vars_used = prep$variables, ClusterVariableName = ClusterVariableName,
    Preprocessing = list(
      ZScoreType = Scaling, Scaling = Scaling, ZScoreObject = prep$ZScoreObject, PCAObject = pca,
      PCVariables = score_vars), PCA = pca,
    fit_plot = PlotClusterFitReview(review$fit_table)), base)
  out$Stability <- Stability
  out[[if (algorithm == "mclust") "ModelInfo_MClust" else "ModelInfo_KMeans"]] <- out$ModelInfo
  out$ModelInfo_PCA <- list(PCAObject = pca, ZScoreType = Scaling, Scaling = Scaling,
    ZScoreObject = prep$ZScoreObject, PCVariables = score_vars)
  out$Specification <- .ClusterSpecification(
    if (algorithm == "mclust") "PCA + Mclust" else "PCA + KMeans",
    prep$variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, review$fit_table, list(), prep$complete_rows)
  class(out) <- c(
    if (algorithm == "mclust") "Pipeline_PCA_MClust" else "Pipeline_PCA_KMeans",
    class(selected))
  out
}

# Expand a reduced-space diagnostic table back to full length, preserving the
# type of each column so factors such as Projection_Fit_Class survive.
.PadIndividual <- function(reduced, indicator_names, n_rows, rows) {
  padded <- lapply(reduced[indicator_names], function(column) {
    filled <- if (is.factor(column)) {
      factor(rep(NA_character_, n_rows), levels = levels(column))
    } else if (is.logical(column)) {
      rep(NA, n_rows)
    } else if (is.character(column)) {
      rep(NA_character_, n_rows)
    } else rep(NA_real_, n_rows)
    filled[rows] <- column
    filled
  })
  dplyr::as_tibble(padded)
}

.ProjectPCACluster <- function(object, new_df, algorithm, ClusterVariableName) {
  prep <- .ProjectClusterNumeric(object, new_df)
  complete <- prep$complete_rows
  if (!any(complete)) {
    individual <- dplyr::tibble(Distance = rep(NA_real_, nrow(prep$data)))
    return(.ClusterOutput(
      prep, integer(), ClusterVariableName, individual, object$ModelInfo))
  }
  pca_projection <- ProjectPCA(
    data = as.data.frame(prep$X[complete, , drop = FALSE]),
    PCAInput = object$Preprocessing$PCAObject, InputType = "PCAObj",
    center = FALSE, scale = FALSE)
  base <- object
  base$vars_used <- object$Preprocessing$PCVariables
  base$Preprocessing <- list(ZScoreType = "None", Scaling = "None", ZScoreObject = NULL)
  class(base) <- if (algorithm == "mclust") "Pipeline_MClust" else "Pipeline_KMeans"
  projected_scores <- if (algorithm == "mclust") {
    ProjectCluster(base, pca_projection$Scores, ClusterVariableName)
  } else ProjectCluster(base, pca_projection$Scores, ClusterVariableName)
  reduced <- projected_scores$ProbFit$individual
  indicator_names <- setdiff(names(reduced), "Cluster")
  individual <- .PadIndividual(
    reduced, indicator_names, nrow(prep$data), complete)
  result <- .ClusterOutput(
    prep, reduced$Cluster, ClusterVariableName, individual, object$ModelInfo)
  individual <- result$ProbFit$individual
  out <- c(list(
    vars_used = object$vars_used, ClusterVariableName = ClusterVariableName,
    Preprocessing = object$Preprocessing), result)
  out$ProbFit <- projected_scores$ProbFit
  out$ProbFit$individual <- individual
  ProjectionFit <- projected_scores$ProjectionFit
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- individual$Projection_Fit_Class
  class(out) <- c(
    if (algorithm == "mclust") "Project_PCAMclust" else "Project_PCAKMeans",
    class(out))
  out
}

#' Fit PCA followed by Mclust
#'
#' @description Best for correlated continuous measures when clustering a
#' lower-dimensional, frozen PCA representation is preferable to raw features.
#' @references Jolliffe IT, Cadima J. Principal component analysis: a review
#' and recent developments. *Phil Trans R Soc A.* 2016;374:20150202.
#' @inheritParams CreateClusterModel_MClust
#' @param pca_variance_threshold Cumulative variance retained by the existing PCA workflow.
#' @return A frozen PCA and Gaussian-mixture pipeline with full-pipeline
#'   bootstrap stability. `ModelInfo$plots` leads with the reduction layer's
#'   `scree` and `loadings`, followed by the mixture-model structure figures in
#'   score space and `profiles` restated in the original measurement scale.
#'   Candidate metrics and the advisory recommendation are in
#'   `ModelInfo$fit_table` and `ModelInfo$AHP`; BIC/ICL/AIC, entropy, and
#'   uncertainty use the same interpretation as Mclust.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' vars_Numeric <- paste0("Var", 1:12)
#' review <- CreateClusterModel_PCA_MClust(
#'   df_Training, vars_Numeric, k_range = 2:4, models = 1,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$scree
#' review$ModelInfo$plots$separation_map
#' review$ModelInfo$plots$profiles
#' }
#' @export
CreateClusterModel_PCA_MClust <- function(data, variables = NULL, method = c("exploratory", "finalize"),
    k_range = 2:10, models = c(1L, 2L, 3L), final_k = NULL, final_model = NULL,
    ZScoreType = NULL, Scaling = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, pca_variance_threshold = 0.85,
    stability_resamples = 0L, stability_seed = seed + 1L, stability_progress = FALSE) {
  resolved_scaling <- .ResolveNumericClusterScaling(ZScoreType, Scaling)
  .CreatePCACluster(data, variables, "mclust", .ClusterMethod(method), k_range, final_k, final_model,
                    resolved_scaling, ClusterVariableName, seed, pca_variance_threshold, models = .ResolveMclustModels(models),
                    stability_resamples = stability_resamples, stability_seed = stability_seed,
                    stability_progress = stability_progress)
}

#' Project data onto a PCA + Mclust model
#' @inheritParams ProjectCluster
#' @return A full-length projection through the frozen scaling, PCA, and
#'   Mclust layers, with posterior membership in `ProbFit` and the projection
#'   triage in `ProjectionFit`.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_PCA_MClust(
#'   df_Training, paste0("Var", 1:12), method = "finalize",
#'   final_k = 4, final_model = 1, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$overall_fit_summary
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_PCA_MClust <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  .ProjectPCACluster(object, new_df, "mclust", ClusterVariableName)
}

#' Fit PCA followed by K-means
#'
#' @description Best for high-dimensional correlated continuous measures with
#' compact clusters in PCA score space.
#' @references Jolliffe IT, Cadima J. *Phil Trans R Soc A.* 2016;374:20150202.
#' @inheritParams CreateClusterModel_KMeans
#' @param final_k Number of clusters for a finalized PCA + K-means solution.
#' @param pca_variance_threshold Cumulative variance retained by the existing PCA workflow.
#' @return A frozen PCA and K-means pipeline with full-pipeline bootstrap
#'   stability. `ModelInfo$plots` leads with the reduction layer's `scree` and
#'   `loadings`, followed by the K-means structure figures in score space
#'   (including the per-participant `silhouette` profile) and `profiles`
#'   restated in the original measurement scale. WSS/BSS, silhouette, and
#'   Calinski-Harabasz use the same interpretation as K-means.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' vars_Numeric <- paste0("Var", 1:12)
#' review <- CreateClusterModel_PCA_KMeans(
#'   df_Training, vars_Numeric, k_range = 2:5, nstart = 10,
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#'
#' # The reduction layer comes first: how many components were retained, and
#' # what each one is actually made of.
#' review$ModelInfo$plots$scree
#' review$ModelInfo$plots$loadings
#'
#' # Then the clustering in that score space.
#' review$ModelInfo$plots$silhouette
#' }
#' @export
CreateClusterModel_PCA_KMeans <- function(data, variables = NULL, method = c("exploratory", "finalize"),
    k_range = 2:10, final_k = NULL,
    ZScoreType = NULL, Scaling = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, nstart = 50L, pca_variance_threshold = 0.85,
    stability_resamples = 0L, stability_seed = seed + 1L, stability_progress = FALSE) {
  resolved_scaling <- .ResolveNumericClusterScaling(ZScoreType, Scaling)
  .CreatePCACluster(data, variables, "kmeans", .ClusterMethod(method), k_range, final_k, NULL,
                    resolved_scaling, ClusterVariableName, seed, pca_variance_threshold, nstart,
                    stability_resamples = stability_resamples, stability_seed = stability_seed,
                    stability_progress = stability_progress)
}

#' Project data onto a PCA + K-means model
#' @inheritParams ProjectCluster
#' @return A full-length projection through frozen scaling, PCA, and centroid
#'   layers, with distances and margins in `ProbFit` and the projection triage
#'   in `ProjectionFit`.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_PCA_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize",
#'   final_k = 4, nstart = 10, stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_PCA_KMeans <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  .ProjectPCACluster(object, new_df, "kmeans", ClusterVariableName)
}

#' Fit a projectable HDBSCAN model
#'
#' @description Best for irregularly shaped numeric clusters, variable density,
#' and data where a meaningful noise/outlier group is expected.
#' @references McInnes L, Healy J, Astels S. hdbscan: Hierarchical density based
#' clustering. *J Open Source Softw.* 2017;2(11):205.
#' @inheritParams CreateClusterModel_MClust
#' @param minPts_range Candidate minimum-points settings in exploratory mode.
#' @param cluster_selection_epsilon_range Candidate extraction epsilon values.
#' @param final_minPts,final_cluster_selection_epsilon Finalized density settings.
#' @return A fitted HDBSCAN model with cluster/noise assignments, membership
#'   probabilities, outlier scores, frozen nearest-core support thresholds, and
#'   bootstrap ARI/Jaccard and noise-recovery metrics. Persistence is
#'   higher-is-better, noise proportion is lower-is-better, and the extracted
#'   class count is data-derived; membership probability and outlier score are
#'   assignment diagnostics rather than candidate fit metrics. Figures sit beside what
#'   they describe: `fit_plot` reviews the density grid;
#'   `ModelInfo$plots` holds `density_review`, `persistence`,
#'   `cluster_persistence`, `separation_map`, and `profiles`; `ModelInfo$FitDiagnostics$plots` holds the
#'   core-distance triad and `outlier_score_by_cluster`; `ProbFit$plots` holds
#'   membership-probability figures; and `Stability$plots` holds bootstrap
#'   agreement, per-cluster recovery, and noise recovery.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' review <- CreateClusterModel_HDBSCAN(
#'   df_Training, c("DensityX", "DensityY"), minPts_range = c(6, 10),
#'   cluster_selection_epsilon_range = c(0, 0.05),
#'   stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$density_review
#' model <- CreateClusterModel_HDBSCAN(
#'   df_Training, c("DensityX", "DensityY"), method = "finalize",
#'   final_minPts = 10, final_cluster_selection_epsilon = 0,
#'   stability_resamples = 2
#' )
#' model$ModelInfo$plots$separation_map
#' model$ModelInfo$plots$cluster_persistence
#' model$ModelInfo$FitDiagnostics$plots$outlier_score_by_cluster
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$nearest_core_support
#' }
#' @export
CreateClusterModel_HDBSCAN <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), minPts_range = 2:10,
    cluster_selection_epsilon_range = c(0, 0.05, 0.10), final_minPts = NULL,
    final_cluster_selection_epsilon = NULL,
    ZScoreType = NULL, Scaling = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, stability_resamples = 0L,
    stability_seed = seed + 1L, stability_progress = FALSE) {
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  method <- .ClusterMethod(method)
  Scaling <- .ResolveNumericClusterScaling(ZScoreType, Scaling)
  prep <- .PrepareClusterNumeric(data, variables, Scaling)
  X <- prep$X[prep$complete_rows, , drop = FALSE]; if (!nrow(X)) stop("No complete rows available for clustering.")
  if (method == "finalize" && (is.null(final_minPts) || is.null(final_cluster_selection_epsilon))) stop("final_minPts and final_cluster_selection_epsilon are required for method = 'finalize'.")
  if (Scaling == "None" && method == "exploratory" && missing(cluster_selection_epsilon_range)) stop("Provide a scale-appropriate cluster_selection_epsilon_range when Scaling = 'None'.")
  grid <- if (method == "finalize") expand.grid(MinPts = final_minPts, Epsilon = final_cluster_selection_epsilon) else expand.grid(MinPts = minPts_range, Epsilon = cluster_selection_epsilon_range)
  fits <- lapply(seq_len(nrow(grid)), function(i) dbscan::hdbscan(X, minPts = grid$MinPts[[i]], cluster_selection_epsilon = grid$Epsilon[[i]]))
  rows <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    sizes <- if (any(f$cluster > 0)) table(f$cluster[f$cluster > 0]) else integer()
    dplyr::tibble(
      Classes = length(sizes), MinPts = grid$MinPts[[i]],
      Epsilon = grid$Epsilon[[i]],
      Persistence = if (length(f$cluster_scores)) mean(f$cluster_scores, na.rm = TRUE) else NA_real_,
      NoiseProportion = mean(f$cluster == 0),
      MeanMembershipProbability = mean(f$membership_prob, na.rm = TRUE),
      MinClusterN = if (length(sizes)) min(sizes) else 0L,
      SizeBalance = if (length(sizes)) min(sizes) / max(sizes) else 0)
  }))
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(seq_len(nrow(grid)), function(i) {
      reference <- rep(NA_integer_, nrow(prep$data)); reference[prep$complete_rows] <- fits[[i]]$cluster
      .ClusterBootstrapStability(prep$data, reference, function(boot, original) {
        fitted <- CreateClusterModel_HDBSCAN(boot, prep$variables, method = "finalize",
          final_minPts = grid$MinPts[[i]], final_cluster_selection_epsilon = grid$Epsilon[[i]],
          ZScoreType = Scaling, seed = seed, stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i,
      list(MinPts = grid$MinPts[[i]], Epsilon = grid$Epsilon[[i]]), noise_label = 0L,
      progress = stability_progress)
    })
    stability_summary <- dplyr::bind_rows(lapply(stabilities, `[[`, "summary"))
    rows <- dplyr::left_join(rows, stability_summary, by = c("MinPts", "Epsilon"))
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, stability_seed)
  }
  review <- .ClusterAHP(rows, higher = c("Persistence", "MeanMembershipProbability", "MinClusterN", "SizeBalance", "ReproducibilityScore", "NoiseSensitivity", "NoiseSpecificity"), lower = "NoiseProportion", setting = "HDBSCAN density setting")
  rows <- review$fit_table; best_i <- if (any(rows$Recommended)) which(rows$Recommended)[1] else 1L
  fit <- fits[[best_i]]; minPts <- grid$MinPts[[best_i]]
  thresholds <- tapply(fit$coredist[fit$cluster > 0], fit$cluster[fit$cluster > 0], stats::quantile, probs = .95, na.rm = TRUE)
  ind <- dplyr::tibble(
    MembershipProbability = rep(NA_real_, nrow(prep$data)),
    OutlierScore = rep(NA_real_, nrow(prep$data)),
    CoreDistance = rep(NA_real_, nrow(prep$data)))
  ind$MembershipProbability[prep$complete_rows] <- fit$membership_prob
  ind$OutlierScore[prep$complete_rows] <- fit$outlier_scores
  ind$CoreDistance[prep$complete_rows] <- fit$coredist
  hdb_review <- tidyr::pivot_longer(
    rows, dplyr::all_of(intersect(
      c("Persistence", "MeanMembershipProbability", "NoiseProportion", "MinClusterN"),
      names(rows))), names_to = "Metric", values_to = "Value")
  p_density_review <- ggplot2::ggplot(
    hdb_review,
    ggplot2::aes(x = .data$MinPts, y = .data$Value,
      color = factor(.data$Epsilon), group = .data$Epsilon)) +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::facet_wrap(
      ~Metric, scales = "free_y",
      labeller = ggplot2::labeller(Metric = .ClusterMetricLabel)) +
    .SciDataColourScale() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Density-setting review",
      subtitle = "Extraction behaviour across the minPts and epsilon grid",
      x = "Minimum points", color = "Epsilon", y = "Fit metric")
  p_cluster_persistence <- if (length(fit$cluster_scores)) {
    ggplot2::ggplot(
      dplyr::tibble(
        Cluster = .ClusterFactor(seq_along(fit$cluster_scores), NULL),
        Persistence = as.numeric(fit$cluster_scores)),
      ggplot2::aes(x = .data$Cluster, y = .data$Persistence)) +
      ggplot2::geom_col(fill = "#3B5BDB") + ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Cluster persistence",
        subtitle = "Stability of each extracted cluster in the condensed tree",
        x = "Cluster", y = "Persistence")
  } else NULL
  review_space <- .ClusterReviewSpace(X)
  ModelInfo <- list(
    hdbscan_model = fit, core_points = X, core_clusters = fit$cluster,
    support_thresholds = thresholds, ReviewSpace = review_space,
    fit_table = rows, AHP = review$AHP,
    final_minPts = minPts, final_cluster_selection_epsilon = grid$Epsilon[[best_i]])
  base <- .ClusterOutput(prep, fit$cluster, ClusterVariableName, ind,
    ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, "CoreDistance", "Core distance")
  ModelInfo$FitDiagnostics$plots$outlier_score_by_cluster <-
    PlotClusterDiagnostic(
      individual, "OutlierScore", "Outlier score by cluster")
  ModelInfo$plots <- .DropNullPlots(list(
    density_review = p_density_review,
    persistence = .ClusterCandidateCurve(
      rows, "Persistence", "Mean cluster persistence by minPts",
      x = "MinPts"),
    cluster_persistence = p_cluster_persistence,
    profiles = PlotClusterProfiles(
      dplyr::filter(base$DataWithClusters, .data[[ClusterVariableName]] > 0),
      ClusterVariableName, prep$variables)))
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(
    individual, "MembershipProbability", "Membership probability")
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(method = method, vars_used = prep$variables, ClusterVariableName = ClusterVariableName,
                Preprocessing = list(ZScoreType = Scaling, Scaling = Scaling, ZScoreObject = prep$ZScoreObject),
                fit_plot = PlotClusterFitReview(rows, x = "MinPts")), base)
  out$Stability <- Stability
  out$ModelInfo_HDBSCAN <- out$ModelInfo
  out$Specification <- .ClusterSpecification(
    "HDBSCAN", prep$variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, grid, list(minPts = minPts,
      epsilon = grid$Epsilon[[best_i]]), prep$complete_rows,
    list(projection = "conservative nearest-training-support assignment"))
  class(out) <- c("Pipeline_HDBSCAN", class(out)); out
}

#' Project data onto an HDBSCAN model using conservative nearest-core support
#' @inheritParams ProjectCluster
#' @return A full-length conservative projection containing cluster/noise
#'   assignments, nearest-core distance, and support flags. `ProjectionFit`
#'   holds the projected distance triad, the frozen-space `separation_map`, the
#'   fit-class summaries and the
#'   method-specific `nearest_core_support` plot.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' model <- CreateClusterModel_HDBSCAN(
#'   df_Training, c("DensityX", "DensityY"), method = "finalize",
#'   final_minPts = 10, final_cluster_selection_epsilon = 0,
#'   stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$plots$nearest_core_support
#' projected$ProjectionFit$plots$training_vs_projected_cluster_occupancy
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_HDBSCAN <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  prep <- .ProjectClusterNumeric(object, new_df); cl <- rep(NA_integer_, nrow(prep$data)); ind <- dplyr::tibble(NearestCoreDistance = rep(NA_real_, nrow(prep$data)), InTrainingSupport = rep(NA, nrow(prep$data)))
  if (any(prep$complete_rows)) {
    core_rows <- which(object$ModelInfo$core_clusters > 0)
    if (!length(core_rows)) {
      cl[prep$complete_rows] <- 0L
      ind$InTrainingSupport[prep$complete_rows] <- FALSE
    } else for (row_i in which(prep$complete_rows)) {
      ds <- sqrt(rowSums((sweep(object$ModelInfo$core_points[core_rows, , drop = FALSE], 2, prep$X[row_i, ], "-"))^2))
      nearest <- core_rows[which.min(ds)]; candidate <- object$ModelInfo$core_clusters[nearest]
      threshold <- object$ModelInfo$support_thresholds[[as.character(candidate)]]
      ind$NearestCoreDistance[row_i] <- min(ds)
      ind$InTrainingSupport[row_i] <- is.finite(threshold) && min(ds) <= threshold
      cl[row_i] <- if (isTRUE(ind$InTrainingSupport[row_i])) candidate else 0L
    }
  }
  support_data <- dplyr::filter(
    dplyr::mutate(ind, Cluster = cl), !is.na(.data$InTrainingSupport))
  p_support <- ggplot2::ggplot(
    support_data,
    ggplot2::aes(x = .data$InTrainingSupport, fill = .data$InTrainingSupport)) +
    ggplot2::geom_bar() + .SciDataFillScale() + ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(labels = c(`FALSE` = "Out of support", `TRUE` = "In support")) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = "Nearest-core support against the frozen training clusters",
      subtitle = paste0(
        "Cases beyond a cluster's frozen core-distance threshold are assigned ",
        "to noise rather than forced into a phenotype"),
      x = "Nearest-core support", y = "Participants")
  base <- .ClusterOutput(
    prep, cl[prep$complete_rows], ClusterVariableName, ind, object$ModelInfo)
  individual <- base$ProbFit$individual
  ProjectionFit <- .ClusterProjectionFit(
    individual, object$ModelInfo$FitDiagnostics$distance_reference,
    "NearestCoreDistance", "Distance to the nearest training core point",
    training_cluster = object$ProbFit$individual$Cluster)
  ProjectionFit$plots$nearest_core_support <- p_support
  ProjectionFit$plots$separation_map <- .ClusterSeparationMap(
    object$ModelInfo$ReviewSpace, prep$X[prep$complete_rows, , drop = FALSE],
    cl[prep$complete_rows], prefix = "Projected")
  out <- c(list(vars_used = object$vars_used, ClusterVariableName = ClusterVariableName, Preprocessing = object$Preprocessing), base)
  out$ProbFit <- .ClusterProbFit(
    individual, "NearestCoreDistance",
    "Distance to the nearest training core point", prefix = "Projected")
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- ProjectionFit$fit_class
  out$ProbFit$individual$Projection_Fit_Class <- ProjectionFit$fit_class
  class(out) <- c("Project_HDBSCAN", class(out)); out
}

.GowerToMedoids <- function(data, medoids, variables, ranges, levels) {
  out <- matrix(NA_real_, nrow(data), nrow(medoids))
  for (i in seq_len(nrow(data))) for (j in seq_len(nrow(medoids))) {
    parts <- vapply(variables, function(v) {
      x <- data[[v]][i]; y <- medoids[[v]][j]
      if (is.na(x) || is.na(y)) return(NA_real_)
      if (is.numeric(data[[v]])) { r <- ranges[[v]]; if (is.na(r) || r == 0) 0 else abs(x - y) / r } else as.numeric(as.character(x) != as.character(y))
    }, numeric(1))
    out[i, j] <- mean(parts, na.rm = TRUE)
  }
  out
}

#' Fit a projectable Gower-distance PAM model for mixed clinical data
#'
#' @description Best for mixed continuous, binary, ordinal, and nominal clinical
#' measures where medoid exemplars are more interpretable than centroids.
#' @references Gower JC. A general coefficient of similarity and some of its
#' properties. *Biometrics.* 1971;27:857-871. Kaufman L, Rousseeuw PJ.
#' *Finding Groups in Data*. Wiley; 1990.
#' @param data Data frame containing numeric, logical, factor, or ordered variables.
#' @param variables Variables used for clustering.
#' @param method Either `"exploratory"`/`"explore"` or `"finalize"`.
#' @param k_range Candidate medoid counts in exploratory mode.
#' @param final_k Final medoid count in finalized mode.
#' @param seed Random seed.
#' @param stability_resamples Number of bootstrap stability refits.
#' @param stability_seed Seed controlling bootstrap resampling.
#' @param stability_progress Whether to print bootstrap progress messages.
#' @param ClusterVariableName Output cluster column name.
#' @return A fitted PAM model with frozen medoids, numeric ranges, categorical
#'   levels, silhouette and medoid-distance metrics in `ModelInfo$fit_table`,
#'   and bootstrap stability. Mean silhouette is higher-is-better and mean
#'   assigned-medoid Gower distance is lower-is-better (a dissimilarity, not a
#'   probability). Figures sit beside what they describe:
#'   `fit_plot` reviews candidates; `ModelInfo$plots` holds `silhouette` (the
#'   per-participant profile from the selected PAM solution),
#'   `silhouette_by_k`, `gower_map`, `medoid_table`,
#'   `categorical_composition`, and `profiles`;
#'   `ModelInfo$FitDiagnostics$plots` holds the medoid-distance triad;
#'   `ProbFit$plots` holds assignment-margin figures; and `Stability$plots`
#'   holds cluster-recovery and complementary stability diagnostics.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Mixed <- c("Var1", "Var2", "Var3", "CatVar1", "CatVar2", "CatVar3")
#' review <- CreateClusterModel_Gower_PAM(
#'   df_Training, vars_Mixed, k_range = 2:5, stability_resamples = 2
#' )
#' review$ModelInfo$fit_table
#' review$ModelInfo$AHP$recommendation
#' review$fit_plot
#' review$ModelInfo$plots$silhouette_by_k
#' model <- CreateClusterModel_Gower_PAM(
#'   df_Training, vars_Mixed, method = "finalize", final_k = 4,
#'   stability_resamples = 2
#' )
#' model$ModelInfo$plots$silhouette
#' model$ModelInfo$plots$gower_map
#' model$ModelInfo$plots$categorical_composition
#' model$ModelInfo$medoids
#' projected <- ProjectCluster(model, df_Projection)
#' projected$ProjectionFit$plots$projection_fit_class_bar
#' }
#' @export
CreateClusterModel_Gower_PAM <- function(data, variables = NULL,
    method = c("exploratory", "finalize"), k_range = 2:10, final_k = NULL,
    ClusterVariableName = "Cluster", seed = 93421L, stability_resamples = 0L,
    stability_seed = seed + 1L, stability_progress = FALSE) {
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Package 'cluster' is required.")
  stability_resamples <- .ValidateClusterStability(
    stability_resamples, stability_seed, stability_progress)
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (is.null(variables)) variables <- names(data)
  if (length(setdiff(variables, names(data)))) stop("data is missing required clustering variables.")
  df_scidr <- data
  method <- .ClusterMethod(method)
  complete <- stats::complete.cases(df_scidr[variables]); if (sum(complete) < 2) stop("Fewer than two complete rows are available.")
  training <- df_scidr[complete, variables, drop = FALSE]
  ks <- if (method == "finalize") { if (is.null(final_k)) stop("final_k is required for method = 'finalize'."); final_k } else k_range
  ks <- sort(unique(as.integer(ks[ks >= 2 & ks < nrow(training)])))
  if (!length(ks)) stop("k_range must contain a value from 2 to n_complete_rows - 1.")
  fits <- lapply(ks, function(k) cluster::pam(cluster::daisy(training, metric = "gower"), k = k, diss = TRUE))
  rows <- dplyr::bind_rows(lapply(seq_along(fits), function(i) dplyr::tibble(
    Classes = ks[[i]], Silhouette = fits[[i]]$silinfo$avg.width,
    MinClusterN = min(fits[[i]]$clusinfo[, "size"]),
    SizeBalance = min(fits[[i]]$clusinfo[, "size"]) /
      max(fits[[i]]$clusinfo[, "size"]),
    MeanMedoidDistance = fits[[i]]$objective[[2]])))
  Stability <- NULL
  if (stability_resamples > 0L) {
    stabilities <- lapply(seq_along(ks), function(i) {
      reference <- rep(NA_integer_, nrow(df_scidr)); reference[complete] <- fits[[i]]$clustering
      .ClusterBootstrapStability(df_scidr, reference, function(boot, original) {
        fitted <- CreateClusterModel_Gower_PAM(boot, variables, method = "finalize",
          final_k = ks[[i]], seed = seed, stability_resamples = 0L)
        ProjectCluster(fitted, original)$ProbFit$individual$Cluster
      }, stability_resamples, stability_seed + i, list(Classes = ks[[i]]),
      progress = stability_progress, preserve_levels = variables)
    })
    stability_summary <- dplyr::bind_rows(lapply(stabilities, `[[`, "summary"))
    rows <- dplyr::left_join(rows, stability_summary, by = "Classes")
    Stability <- .CombineClusterStabilities(stabilities, stability_resamples, stability_seed)
  }
  review <- .ClusterAHP(rows, higher = c("Silhouette", "MinClusterN", "SizeBalance", "ReproducibilityScore"), lower = "MeanMedoidDistance", setting = "Gower/PAM k")
  rows <- review$fit_table; best_i <- if (any(rows$Recommended)) which(rows$Recommended)[1] else 1L
  fit <- fits[[best_i]]; k <- ks[[best_i]]
  ranges <- lapply(training, function(x) if (is.numeric(x)) diff(range(x, na.rm = TRUE)) else NA_real_)
  medoids <- training[fit$id.med, , drop = FALSE]
  prep <- list(data = df_scidr, complete_rows = complete, variables = variables)
  d <- .GowerToMedoids(training, medoids, variables, ranges, lapply(training, levels))
  ind <- dplyr::tibble(DistanceToMedoid = rep(NA_real_, nrow(prep$data)), AssignmentMargin = rep(NA_real_, nrow(prep$data)))
  ind$DistanceToMedoid[complete] <- apply(d, 1, min); ind$AssignmentMargin[complete] <- apply(d, 1, function(x) diff(sort(x))[1])
  categorical <- variables[!vapply(training, is.numeric, logical(1))]
  numeric_variables <- variables[vapply(training, is.numeric, logical(1))]
  p_categorical <- PlotClusterComposition(training, categorical, fit$clustering)
  p_categorical_by_cluster <- PlotClusterComposition(
    training, categorical, fit$clustering, facet_by = "cluster")
  p_categorical_enrichment <- PlotClusterComposition(
    training, categorical, fit$clustering, style = "enrichment")
  # The Gower configuration is frozen so the training map stays the reference
  # picture; projected cases are reviewed against the medoids instead, because
  # classical scaling has no exact out-of-sample extension.
  gower_coordinates <- stats::cmdscale(
    cluster::daisy(training, metric = "gower"), k = 2L)
  gower_assignment <- as.data.frame(gower_coordinates)
  names(gower_assignment) <- c("GowerDimension1", "GowerDimension2")
  gower_assignment[[ClusterVariableName]] <- fit$clustering
  p_gower_map <- PlotClusterMap(
    gower_assignment, "GowerDimension1", "GowerDimension2", ClusterVariableName,
    title = "Training cluster review map",
    subtitle = "Classical scaling of the training Gower dissimilarities",
    xlab = "Gower dimension 1", ylab = "Gower dimension 2")
  ModelInfo <- list(
    pam_model = fit, medoids = medoids, ranges = ranges,
    levels = lapply(training, levels), GowerCoordinates = gower_assignment,
    silhouette = fit$silinfo$widths, fit_table = rows, AHP = review$AHP,
    final_k = k)
  base <- .ClusterOutput(prep, fit$clustering, ClusterVariableName, ind,
    ModelInfo, Stability = Stability)
  individual <- base$ProbFit$individual
  ModelInfo$FitDiagnostics <- .ClusterFitDiagnostics(
    individual, "DistanceToMedoid", "Gower distance to assigned medoid")
  ModelInfo$plots <- .DropNullPlots(list(
    silhouette = PlotClusterSilhouette(
      fit$silinfo$widths, paste0("Silhouette profile (k = ", k, ")")),
    silhouette_by_k = .ClusterCandidateCurve(
      rows, "Silhouette", "Average silhouette width by k",
      "Average silhouette width"),
    medoid_distance_by_k = .ClusterCandidateCurve(
      rows, "MeanMedoidDistance", "Mean distance to medoid by k"),
    categorical_composition = p_categorical,
    categorical_composition_by_cluster = p_categorical_by_cluster,
    categorical_enrichment = p_categorical_enrichment,
    profiles = if (length(numeric_variables)) {
      PlotClusterProfiles(base$DataWithClusters, ClusterVariableName, numeric_variables)
    } else p_categorical))
  base$ModelInfo <- ModelInfo
  base$ProbFit <- .ClusterProbFit(
    individual, "AssignmentMargin",
    "Gower distance margin to the second-nearest medoid")
  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)
  out <- c(list(method = method, vars_used = variables, ClusterVariableName = ClusterVariableName, Preprocessing = list(Scaling = "Gower", ZScoreObject = NULL), fit_plot = PlotClusterFitReview(rows)), base)
  out$Stability <- Stability
  out$ModelInfo_GowerPAM <- out$ModelInfo
  out$Specification <- .ClusterSpecification(
    "GowerPAM", variables, list(seed = seed, stability_seed = stability_seed),
    out$Preprocessing, dplyr::select(rows, dplyr::any_of("Classes")),
    list(k = k), complete,
    list(distance_metric = "Gower distance to frozen medoid", levels = ModelInfo$levels))
  class(out) <- c("Pipeline_Gower_PAM", class(out)); out
}

#' Project data onto a Gower + PAM model
#' @inheritParams ProjectCluster
#' @return A full-length projection using frozen Gower ranges, levels, and
#'   medoids. `ProbFit` holds the distance to the assigned medoid and the
#'   margin to the second-nearest medoid; `ProjectionFit` triages every case
#'   against the frozen training distance reference and holds the projected
#'   distance triad, a `medoid_map` placing cases by their two nearest medoid
#'   distances, fit-class summaries, and the
#'   fit-class bar.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' df_Projection <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Projection")
#' vars_Mixed <- c("Var1", "Var2", "Var3", "CatVar1", "CatVar2", "CatVar3")
#' model <- CreateClusterModel_Gower_PAM(
#'   df_Training, vars_Mixed, method = "finalize", final_k = 4,
#'   stability_resamples = 2
#' )
#' projected <- ProjectCluster(model, df_Projection)
#' head(projected$ProbFit$individual)
#' projected$ProjectionFit$plots$medoid_map
#' projected$ProjectionFit$plots$training_vs_projected_cluster_occupancy
#' }
#' @noRd
#' @export
ProjectCluster.Pipeline_Gower_PAM <- function(object, new_df, ClusterVariableName = object$ClusterVariableName, ...) {
  if (length(setdiff(object$vars_used, names(new_df)))) stop("new_df is missing required clustering variables.")
  out_of_support <- dplyr::bind_rows(lapply(object$vars_used, function(variable) {
    known <- object$ModelInfo$levels[[variable]]
    if (is.null(known) || is.numeric(new_df[[variable]])) return(dplyr::tibble())
    observed <- as.character(new_df[[variable]])
    bad <- !is.na(observed) & !observed %in% known
    if (!any(bad)) return(dplyr::tibble())
    dplyr::tibble(RowID = which(bad), Variable = variable,
      Value = observed[bad], Issue = "category absent from training")
  }))
  new_df_supported <- new_df
  if (nrow(out_of_support)) {
    for (variable in unique(out_of_support$Variable)) {
      rows <- out_of_support$RowID[out_of_support$Variable == variable]
      new_df_supported[[variable]][rows] <- NA
    }
    warning("Projected Gower/PAM data contain categories absent from training; affected rows are out of support.", call. = FALSE)
  }
  prep <- list(data = new_df_supported, complete_rows = stats::complete.cases(new_df_supported[object$vars_used]),
               variables = object$vars_used)
  d <- .GowerToMedoids(prep$data[prep$complete_rows, object$vars_used, drop = FALSE], object$ModelInfo$medoids, object$vars_used, object$ModelInfo$ranges, object$ModelInfo$levels)
  cl <- rep(NA_integer_, nrow(prep$data)); cl[prep$complete_rows] <- max.col(-d)
  ind <- dplyr::tibble(DistanceToMedoid = rep(NA_real_, nrow(prep$data)), AssignmentMargin = rep(NA_real_, nrow(prep$data))); ind$DistanceToMedoid[prep$complete_rows] <- apply(d, 1, min); ind$AssignmentMargin[prep$complete_rows] <- apply(d, 1, function(x) diff(sort(x))[1])
  # Classical scaling cannot place new cases in the frozen training
  # configuration, so projected cases are reviewed in medoid-distance space:
  # how close each case sits to its own medoid versus its closest rival.
  sorted_distances <- t(apply(d, 1, sort))
  distance_plot_data <- data.frame(
    NearestMedoidDistance = sorted_distances[, 1],
    SecondNearestMedoidDistance = sorted_distances[, min(2L, ncol(sorted_distances))],
    Cluster = cl[prep$complete_rows])
  p_medoid_map <- PlotClusterMap(
    distance_plot_data, "NearestMedoidDistance", "SecondNearestMedoidDistance",
    "Cluster", title = "Projected cluster review map",
    subtitle = paste0(
      "Frozen medoid-distance space; points near the diagonal sit between ",
      "two phenotypes"),
    xlab = "Gower distance to assigned medoid",
    ylab = "Gower distance to second-nearest medoid") +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", color = "grey40")
  base <- .ClusterOutput(
    prep, cl[prep$complete_rows], ClusterVariableName, ind, object$ModelInfo)
  individual <- base$ProbFit$individual
  ProjectionFit <- .ClusterProjectionFit(
    individual, object$ModelInfo$FitDiagnostics$distance_reference,
    "DistanceToMedoid", "Gower distance to assigned medoid",
    training_cluster = object$ProbFit$individual$Cluster)
  ProjectionFit$out_of_support <- out_of_support
  ProjectionFit$plots$medoid_map <- p_medoid_map
  out <- c(list(vars_used = object$vars_used, ClusterVariableName = ClusterVariableName, Preprocessing = object$Preprocessing), base)
  out$ProbFit <- .ClusterProbFit(
    individual, "AssignmentMargin",
    "Gower distance margin to the second-nearest medoid", prefix = "Projected")
  out$ProjectionFit <- ProjectionFit
  out$DataWithClusters$Projection_Fit_Class <- ProjectionFit$fit_class
  out$ProbFit$individual$Projection_Fit_Class <- ProjectionFit$fit_class
  class(out) <- c("Project_GowerPAM", class(out)); out
}
