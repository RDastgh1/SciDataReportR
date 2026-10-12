#' Compare projectable clustering methods with common participant subsamples
#'
#' @description Fits each requested candidate on the same full reference cohort,
#'   then on identical 90% participant subsamples. Every refit is projected to
#'   the full complete reference cohort before participant-level agreement is
#'   calculated. This makes ARI comparable across clustering methods rather than
#'   comparing only participants retained by a particular refit.
#'
#' @details Candidate rows use `Method` values `"som_mclust"`, `"mclust"`,
#'   `"kmeans"`, `"hdbscan"`, or `"ward_rf"`. K-based candidates require
#'   `Classes`; mixture candidates also require `Model`. SOM candidates use
#'   `NodeWeighting`, `SOM_XDim`, `SOM_YDim`, and `SOM_Topology`; HDBSCAN uses
#'   `MinPts` and `Epsilon`; and Ward-plus-RF uses `WardN` and `RF_Trees`.
#'   Full benchmarks near 20,000 participants can be slow and memory-intensive.
#'   HDBSCAN and Ward + RF allocate pairwise structures and are automatically
#'   serial under the memory policy; use a small, reviewed candidate table for
#'   large cohorts.
#'
#' @param data A data frame containing the clustering variables.
#' @param variables Numeric clustering variables.
#' @param candidates A candidate specification data frame. `NULL` constructs the
#'   planned SOM/Mclust, Mclust, K-means, HDBSCAN, and Ward-plus-random-forest
#'   candidate set. See `DefaultClusterBenchmarkCandidates()`.
#' @param stability_resamples Number of common 90% subsample refits. Defaults to
#'   10 for screening.
#' @param stability_seed Seed used to generate the shared participant draws.
#' @param stability_cores Optional number of parallel workers.
#' @param subsample_fraction Proportion of complete reference participants
#'   retained in each refit.
#' @param keep_models Logical; retain full-data reference models in the screen
#'   result. Defaults to `FALSE` to avoid holding one large model per candidate.
#'   Confirmation always retains its selected projectable models.
#'
#' @return A `ClusterMethodStabilityBenchmark` object with `decision_table`,
#'   `replicates`, shared `subsamples`, `participant_inclusion`, and full-data
#'   `models`. With `keep_models = FALSE`, `models` is an empty list.
#' @export
CompareClusterMethodStability <- function(data, variables,
    candidates = NULL, stability_resamples = 10L, stability_seed = 934522L,
    stability_cores = NULL, subsample_fraction = .90, keep_models = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is.character(variables) || !length(variables) ||
      length(setdiff(variables, names(data)))) {
    stop("variables must identify one or more columns in data.")
  }
  if (!is.numeric(stability_resamples) || length(stability_resamples) != 1L ||
      stability_resamples < 1L || stability_resamples != as.integer(stability_resamples)) {
    stop("stability_resamples must be a positive integer.")
  }
  if (!is.numeric(subsample_fraction) || length(subsample_fraction) != 1L ||
      !is.finite(subsample_fraction) || subsample_fraction <= 0 ||
      subsample_fraction >= 1) {
    stop("subsample_fraction must be a single value strictly between zero and one.")
  }
  if (!is.logical(keep_models) || length(keep_models) != 1L || is.na(keep_models)) {
    stop("keep_models must be TRUE or FALSE.")
  }
  df_scidr <- .AddClusterRowID(data)
  complete_rows <- which(stats::complete.cases(df_scidr[, variables, drop = FALSE]))
  if (length(complete_rows) < 3L) stop("At least three complete participants are required.")
  candidates <- .ValidateClusterBenchmarkCandidates(candidates)

  set.seed(stability_seed)
  subsamples <- lapply(seq_len(stability_resamples), function(replicate) {
    sampled <- sort(sample(complete_rows,
      max(2L, floor(length(complete_rows) * subsample_fraction)), replace = FALSE))
    dplyr::tibble(Replicate = replicate, .row_id = df_scidr$.row_id[sampled])
  })
  df_subsamples <- dplyr::bind_rows(subsamples)

  models <- list()
  reference <- list()
  rows <- list()
  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, , drop = FALSE]
    fitted <- tryCatch(.FitBenchmarkCandidate(df_scidr, variables, candidate,
      seed = stability_seed + i), error = function(e) e)
    candidate_id <- candidate$CandidateID[[1]]
    if (inherits(fitted, "error")) {
      rows[[i]] <- .BenchmarkFailureRow(candidate, conditionMessage(fitted))
      next
    }
    projected <- .ProjectBenchmarkCandidate(fitted, df_scidr, candidate)
    labels <- projected$Cluster[match(df_scidr$.row_id, projected$.row_id)]
    non_noise <- sort(unique(labels[!is.na(labels) & labels != 0L]))
    applicable <- candidate$Method[[1]] != "hdbscan" || length(non_noise) >= 2L
    if (isTRUE(keep_models)) models[[candidate_id]] <- fitted
    reference[[candidate_id]] <- list(labels = labels, applicable = applicable,
      native = .BenchmarkNativeMetrics(fitted, candidate),
      non_noise_clusters = length(non_noise))
    rm(projected)
    if (!isTRUE(keep_models)) rm(fitted)
  }

  run_candidate <- function(candidate_index) {
    candidate <- candidates[candidate_index, , drop = FALSE]
    candidate_id <- candidate$CandidateID[[1]]
    ref <- reference[[candidate_id]]
    if (is.null(ref)) return(NULL)
    memory_method <- .BenchmarkMemoryMethod(candidate$Method[[1]])
    som_nodes <- if (identical(candidate$Method[[1]], "som_mclust")) {
      as.integer(candidate$SOM_XDim[[1]] * candidate$SOM_YDim[[1]])
    } else NULL
    worker_settings <- .ResolveClusterStabilityCores(
      stability_cores, stability_resamples, n = length(complete_rows),
      p = length(variables), method = memory_method, som_nodes = som_nodes,
      ward_n = if (identical(candidate$Method[[1]], "ward_rf")) candidate$WardN[[1]] else NULL,
      subsample_fraction = subsample_fraction)
    if (!isTRUE(ref$applicable)) {
      return(list(replicates = dplyr::tibble(
        CandidateID = candidate_id, Replicate = seq_len(stability_resamples),
        Status = "not_applicable", ARI = NA_real_, ParticipantsCompared = length(complete_rows),
        Error = "Fewer than two non-noise reference clusters.", JaccardMean = NA_real_),
        assignments = list(), worker_settings = worker_settings))
    }
    results <- .ClusterParallelMap(seq_len(stability_resamples), function(replicate) {
      sampled_ids <- subsamples[[replicate]]$.row_id
      sampled_rows <- match(sampled_ids, df_scidr$.row_id)
      fitted <- tryCatch(.FitBenchmarkCandidate(df_scidr[sampled_rows, , drop = FALSE],
        variables, candidate, seed = stability_seed + candidate_index * 10000L + replicate),
        error = function(e) e)
      if (inherits(fitted, "error")) {
        return(list(replicate = dplyr::tibble(CandidateID = candidate_id,
          Replicate = replicate, Status = "failed", ARI = NA_real_,
          ParticipantsCompared = NA_integer_, Error = conditionMessage(fitted)), assignment = NULL))
      }
      projected <- tryCatch(.ProjectBenchmarkCandidate(fitted, df_scidr, candidate),
        error = function(e) e)
      if (inherits(projected, "error")) {
        rm(projected, fitted)
        gc(FALSE)
        return(list(replicate = dplyr::tibble(CandidateID = candidate_id,
          Replicate = replicate, Status = "failed", ARI = NA_real_,
          ParticipantsCompared = NA_integer_, Error = conditionMessage(projected)), assignment = NULL))
      }
      refit <- projected$Cluster[match(df_scidr$.row_id, projected$.row_id)]
      valid <- !is.na(ref$labels) & !is.na(refit)
      ari <- if (sum(valid) > 1L) mclust::adjustedRandIndex(ref$labels[valid], refit[valid]) else NA_real_
      jaccard <- .ClusterJaccard(ref$labels[valid], refit[valid],
        noise_label = if (candidate$Method[[1]] == "hdbscan") 0L else NULL)
      rm(projected, fitted)
      gc(FALSE)
      list(replicate = dplyr::tibble(CandidateID = candidate_id,
        Replicate = replicate, Status = "success", ARI = ari,
        ParticipantsCompared = sum(valid), Error = NA_character_,
        JaccardMean = if (nrow(jaccard)) mean(jaccard$Jaccard, na.rm = TRUE) else NA_real_),
      assignment = refit)
    }, cores = worker_settings$resolved)
    list(replicates = dplyr::bind_rows(lapply(results, `[[`, "replicate")),
      assignments = lapply(results, `[[`, "assignment"),
      worker_settings = worker_settings)
  }
  benchmark_results <- lapply(seq_len(nrow(candidates)), run_candidate)

  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, , drop = FALSE]
    candidate_id <- candidate$CandidateID[[1]]
    ref <- reference[[candidate_id]]
    if (is.null(ref)) next
    result <- benchmark_results[[i]]
    replicates <- result$replicates
    successful <- dplyr::filter(replicates, .data$Status == "success")
    assignments <- Filter(Negate(is.null), result$assignments)
    diagnostics <- .ClusterStabilityDiagnostics(ref$labels, assignments,
      valid_reference = complete_rows, row_id = df_scidr$.row_id,
      noise_label = if (candidate$Method[[1]] == "hdbscan") 0L else NULL)
    inclusion <- diagnostics$participant_inclusion %>%
      dplyr::mutate(CandidateID = candidate_id)
    reference_sizes <- table(ref$labels[complete_rows])
    rows[[i]] <- dplyr::bind_cols(candidate,
      dplyr::tibble(Eligibility = if (isTRUE(ref$applicable)) "eligible" else "not_applicable",
        NonNoiseClusters = ref$non_noise_clusters,
        SmallestClusterProportion = if (length(reference_sizes)) min(reference_sizes) / sum(reference_sizes) else NA_real_,
        StabilityARI_Mean = .BenchmarkSafeMean(successful$ARI),
        StabilityARI_P05 = .BenchmarkSafeQuantile(successful$ARI, .05),
        StabilityJaccard_Mean = .BenchmarkSafeMean(successful$JaccardMean),
        StabilitySuccessRate = mean(replicates$Status == "success"),
        ResolvedWorkers = result$worker_settings$resolved,
        WorkerPolicy = result$worker_settings$worker_policy,
        EstimatedBytesPerWorker = result$worker_settings$bytes_per_worker,
        EstimatedTotalBytes = result$worker_settings$estimated_total_bytes,
        SafeMemoryBudgetBytes = result$worker_settings$safe_memory_budget_bytes,
        Persistence = .BenchmarkMetricValue(ref$native, "Persistence"),
        NoiseProportion = .BenchmarkMetricValue(ref$native, "NoiseProportion"),
        MeanMembershipProbability = .BenchmarkMetricValue(ref$native, "MeanMembershipProbability"),
        NativeMetrics = list(ref$native),
        ParticipantInclusion = list(inclusion)))
  }
  decision_table <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(StabilityRank = dplyr::if_else(.data$Eligibility == "eligible",
      dplyr::dense_rank(dplyr::desc(.data$StabilityARI_Mean)), NA_integer_)) %>%
    dplyr::arrange(.data$StabilityRank, dplyr::desc(.data$StabilityARI_P05),
      dplyr::desc(.data$StabilityJaccard_Mean))
  participant_inclusion <- dplyr::bind_rows(Filter(is.data.frame,
    lapply(decision_table$ParticipantInclusion, identity)))
  decision_table <- dplyr::select(decision_table, -dplyr::all_of("ParticipantInclusion"))
  replicate_tables <- lapply(Filter(Negate(is.null), benchmark_results), `[[`, "replicates")
  out <- list(decision_table = decision_table,
    replicates = dplyr::bind_rows(replicate_tables),
    subsamples = df_subsamples, participant_inclusion = participant_inclusion,
    models = models, candidates = candidates, data = data, variables = variables,
    settings = list(scoring_cohort = "all_reference_participants",
      perturbation = "90% subsample without replacement",
      refits_requested = as.integer(stability_resamples),
      subsample_fraction = subsample_fraction, stability_seed = stability_seed,
      stability_cores = if (is.null(stability_cores)) "automatic" else as.integer(stability_cores),
      stability_cores_resolved = max(vapply(benchmark_results, function(result) {
        if (is.null(result)) 1L else result$worker_settings$resolved
      }, integer(1))),
      worker_estimates = stats::setNames(lapply(benchmark_results, function(result) {
        if (is.null(result)) NULL else result$worker_settings
      }), candidates$CandidateID),
      keep_models = keep_models,
      ranking = "StabilityARI_Mean, then StabilityARI_P05 and StabilityJaccard_Mean",
      automatic_selection = FALSE))
  class(out) <- "ClusterMethodStabilityBenchmark"
  out
}

#' Confirm method-specific stability leaders
#'
#' @param benchmark A `ClusterMethodStabilityBenchmark` object.
#' @param candidates Optional candidate IDs or decision-table rows. If omitted,
#'   confirms the ARI-leading eligible candidate within each method variant.
#' @param confirm_resamples Number of common refits, default 30.
#'
#' @return A `ClusterMethodStabilityConfirmation` object.
#' @export
ConfirmClusterMethodStability <- function(benchmark, candidates = NULL,
    confirm_resamples = 30L) {
  if (!inherits(benchmark, "ClusterMethodStabilityBenchmark")) {
    stop("benchmark must be returned by CompareClusterMethodStability().")
  }
  if (!is.numeric(confirm_resamples) || length(confirm_resamples) != 1L ||
      confirm_resamples < 1L || confirm_resamples != as.integer(confirm_resamples)) {
    stop("confirm_resamples must be a positive integer.")
  }
  table <- dplyr::filter(benchmark$decision_table, .data$Eligibility == "eligible")
  if (is.null(candidates)) {
    table <- table %>% dplyr::group_by(.data$MethodVariant) %>%
      dplyr::slice_min(.data$StabilityRank, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else if (is.data.frame(candidates)) {
    table <- dplyr::semi_join(table, candidates, by = "CandidateID")
  } else {
    table <- dplyr::filter(table, .data$CandidateID %in% candidates)
  }
  if (!nrow(table)) stop("No eligible candidates were selected for confirmation.")
  confirmed <- CompareClusterMethodStability(benchmark$data, benchmark$variables,
    candidates = dplyr::select(table, dplyr::all_of(names(benchmark$candidates))),
    stability_resamples = confirm_resamples,
    stability_seed = benchmark$settings$stability_seed,
    stability_cores = if (identical(benchmark$settings$stability_cores, "automatic")) NULL else benchmark$settings$stability_cores,
    subsample_fraction = benchmark$settings$subsample_fraction,
    keep_models = TRUE)
  out <- list(decision_table = confirmed$decision_table, replicates = confirmed$replicates,
    subsamples = confirmed$subsamples, participant_inclusion = confirmed$participant_inclusion,
    models = confirmed$models, settings = confirmed$settings, screen = benchmark)
  class(out) <- "ClusterMethodStabilityConfirmation"
  out
}

#' Construct default common-method benchmark candidates
#'
#' @return A candidate specification tibble for
#'   `CompareClusterMethodStability()`.
#' @export
DefaultClusterBenchmarkCandidates <- function() {
  ks <- 3:5
  som <- tidyr::expand_grid(Method = "som_mclust", NodeWeighting = c("uniform", "occupancy"),
    Classes = ks) %>% dplyr::mutate(Model = 2L, SOM_XDim = 10L, SOM_YDim = 10L,
      SOM_Topology = "hexagonal")
  mclust <- tidyr::expand_grid(Method = "mclust", Model = 1:3, Classes = ks)
  kmeans <- dplyr::tibble(Method = "kmeans", Classes = ks)
  hdbscan <- dplyr::tibble(Method = "hdbscan", MinPts = c(20L, 50L, 100L),
    Epsilon = c(0, .05, .10))
  ward <- dplyr::tibble(Method = "ward_rf", Classes = ks, WardN = 6000L, RF_Trees = 500L)
  dplyr::bind_rows(som, mclust, kmeans, hdbscan, ward) %>%
    dplyr::mutate(CandidateID = paste0(.data$Method, "_", dplyr::row_number())) %>%
    .ValidateClusterBenchmarkCandidates()
}

.BenchmarkCandidateColumns <- c("CandidateID", "Method", "MethodVariant", "Classes", "Model",
  "NodeWeighting", "SOM_XDim", "SOM_YDim", "SOM_Topology", "Assignment",
  "SoftNodes", "MinPts", "Epsilon", "WardN", "RF_Trees")

.BenchmarkMemoryMethod <- function(method) {
  switch(method,
    som_mclust = "som_mclust",
    mclust = "mclust",
    kmeans = "kmeans",
    hdbscan = "hdbscan",
    ward_rf = "ward_rf",
    "generic")
}

.ValidateClusterBenchmarkCandidates <- function(candidates) {
  if (is.null(candidates)) candidates <- DefaultClusterBenchmarkCandidates()
  if (!is.data.frame(candidates) || !all(c("CandidateID", "Method") %in% names(candidates))) {
    stop("candidates must be a data frame with CandidateID and Method columns.")
  }
  candidates <- dplyr::as_tibble(candidates)
  if (anyDuplicated(candidates$CandidateID) || anyNA(candidates$CandidateID)) {
    stop("CandidateID values must be present and unique.")
  }
  supported <- c("som_mclust", "mclust", "kmeans", "hdbscan", "ward_rf")
  if (any(!candidates$Method %in% supported)) {
    stop("Method must be one of: ", paste(supported, collapse = ", "), ".")
  }
  for (column in setdiff(.BenchmarkCandidateColumns, names(candidates))) candidates[[column]] <- NA
  candidates <- candidates[, .BenchmarkCandidateColumns]
  candidates$MethodVariant <- dplyr::case_when(
    candidates$Method == "som_mclust" ~ paste0("SOM + Mclust (", candidates$NodeWeighting, ")"),
    candidates$Method == "ward_rf" ~ "Ward + RF surrogate",
    candidates$Method == "mclust" ~ "Mclust",
    candidates$Method == "kmeans" ~ "K-means",
    TRUE ~ "HDBSCAN")
  required_k <- candidates$Method %in% c("som_mclust", "mclust", "kmeans", "ward_rf")
  if (any(is.na(candidates$Classes[required_k]))) stop("Classes is required for k-based candidates.")
  if (any(candidates$Method %in% c("som_mclust", "mclust") & is.na(candidates$Model))) {
    stop("Model is required for SOM + Mclust and Mclust candidates.")
  }
  if (any(candidates$Method == "hdbscan" & is.na(candidates$MinPts))) stop("MinPts is required for HDBSCAN candidates.")
  candidates
}

.FitBenchmarkCandidate <- function(data, variables, candidate, seed) {
  method <- candidate$Method[[1]]
  if (method == "som_mclust") {
    assignment <- if (is.na(candidate$Assignment[[1]])) "soft" else candidate$Assignment[[1]]
    soft_n_nodes <- if (is.na(candidate$SoftNodes[[1]])) 3L else as.integer(candidate$SoftNodes[[1]])
    return(CreateClusterModel_SOM_MClust(data, variables,
    method = "finalize", final_k = candidate$Classes[[1]], final_model = candidate$Model[[1]],
    som_xdim = candidate$SOM_XDim[[1]], som_ydim = candidate$SOM_YDim[[1]],
    som_topo = candidate$SOM_Topology[[1]], node_weighting = candidate$NodeWeighting[[1]],
    assignment = assignment, soft_n_nodes = soft_n_nodes, lpa_min_node_n = 1L,
    seed_som = seed, Relabel = FALSE, stability_resamples = 0L))
  }
  if (method == "mclust") return(CreateClusterModel_MClust(data, variables,
    method = "finalize", final_k = candidate$Classes[[1]], final_model = candidate$Model[[1]],
    seed = seed, stability_resamples = 0L))
  if (method == "kmeans") return(.FitBenchmarkKMeans(data, variables, candidate, seed))
  if (method == "hdbscan") return(CreateClusterModel_HDBSCAN(data, variables,
    method = "finalize", final_minPts = candidate$MinPts[[1]],
    final_cluster_selection_epsilon = candidate$Epsilon[[1]], seed = seed,
    stability_resamples = 0L))
  .FitWardRFSurrogate(data, variables, candidate, seed)
}

.FitBenchmarkKMeans <- function(data, variables, candidate, seed) {
  # The ordinary constructor computes an O(n^2) full silhouette diagnostic.
  # Benchmark refits need projected assignments, centroids, and WSS, not that
  # diagnostic, so this deliberately uses the same scaled-space fit directly.
  prep <- .PrepareClusterNumeric(data, variables, "Center and Scale")
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  set.seed(seed)
  fit <- stats::kmeans(X, centers = candidate$Classes[[1]], nstart = 50L)
  structure(list(vars_used = prep$variables,
    Preprocessing = list(Scaling = prep$Scaling, ZScoreObject = prep$ZScoreObject),
    centers = fit$centers, fit = fit,
    native_metrics = dplyr::tibble(WSS = fit$tot.withinss,
      CalinskiHarabasz = (fit$betweenss / (length(fit$size) - 1L)) /
        (fit$tot.withinss / (nrow(X) - length(fit$size))),
      MinClusterN = min(fit$size))), class = "BenchmarkKMeans")
}

.ProjectBenchmarkKMeans <- function(model, data) {
  prep <- .ProjectClusterNumeric(model, data)
  cluster <- rep(NA_integer_, nrow(prep$data))
  X <- prep$X[prep$complete_rows, , drop = FALSE]
  if (nrow(X)) {
    distances <- sapply(seq_len(nrow(model$centers)), function(cluster_index) {
      rowSums((sweep(X, 2L, model$centers[cluster_index, ], "-"))^2)
    })
    if (is.null(dim(distances))) distances <- matrix(distances, ncol = 1L)
    cluster[prep$complete_rows] <- max.col(-distances)
  }
  dplyr::tibble(.row_id = prep$data$.row_id, Cluster = cluster)
}

.FitWardRFSurrogate <- function(data, variables, candidate, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' is required for ward_rf.")
  df_scidr <- .AddClusterRowID(data)
  complete_rows <- which(stats::complete.cases(df_scidr[, variables, drop = FALSE]))
  X <- scale(df_scidr[complete_rows, variables, drop = FALSE])
  set.seed(seed)
  ward_rows <- sort(sample(seq_len(nrow(X)), min(candidate$WardN[[1]], nrow(X)), replace = FALSE))
  hclust_fn <- if (requireNamespace("fastcluster", quietly = TRUE)) fastcluster::hclust else stats::hclust
  labels <- stats::cutree(hclust_fn(stats::dist(X[ward_rows, , drop = FALSE]), method = "ward.D2"),
    k = candidate$Classes[[1]])
  forest <- ranger::ranger(x = as.data.frame(X[ward_rows, , drop = FALSE]), y = factor(labels),
    num.trees = candidate$RF_Trees[[1]], seed = seed)
  structure(list(center = attr(X, "scaled:center"), scale_sd = attr(X, "scaled:scale"),
    forest = forest, variables = variables, ward_n = length(ward_rows), ward_seed = seed,
    rf_seed = seed, DataWithClusters = .ProjectWardRFSurrogate(list(
      center = attr(X, "scaled:center"), scale_sd = attr(X, "scaled:scale"), forest = forest,
      variables = variables), df_scidr)), class = "WardRFSurrogate")
}

.ProjectWardRFSurrogate <- function(model, data) {
  df_scidr <- .AddClusterRowID(data)
  X <- sweep(sweep(as.matrix(df_scidr[, model$variables, drop = FALSE]), 2L, model$center, "-"),
    2L, model$scale_sd, "/")
  cluster <- rep(NA_integer_, nrow(df_scidr))
  complete <- stats::complete.cases(X)
  if (any(complete)) cluster[complete] <- as.integer(as.character(
    predict(model$forest, as.data.frame(X[complete, , drop = FALSE]))$predictions))
  dplyr::tibble(.row_id = df_scidr$.row_id, Cluster = cluster)
}

.ProjectBenchmarkCandidate <- function(model, data, candidate) {
  if (candidate$Method[[1]] == "ward_rf") return(.ProjectWardRFSurrogate(model, data))
  if (candidate$Method[[1]] == "kmeans" && inherits(model, "BenchmarkKMeans")) {
    return(.ProjectBenchmarkKMeans(model, data))
  }
  projected <- ProjectCluster(model, data)$DataWithClusters
  df_scidr <- .AddClusterRowID(projected)
  dplyr::select(df_scidr, .data$.row_id, .data$Cluster)
}

.BenchmarkNativeMetrics <- function(model, candidate) {
  if (candidate$Method[[1]] == "ward_rf") return(dplyr::tibble(
    WardN = model$ward_n, WardSeed = model$ward_seed, RFSeed = model$rf_seed,
    ProjectionPolicy = "Ward.D2 labels projected by random-forest surrogate"))
  if (candidate$Method[[1]] == "kmeans" && inherits(model, "BenchmarkKMeans")) {
    return(model$native_metrics)
  }
  table <- if (candidate$Method[[1]] == "som_mclust") model$ModelInfo_MClust$fit_table else model$ModelInfo$fit_table
  if (is.null(table) || !nrow(table)) return(dplyr::tibble())
  dplyr::select(table[1, , drop = FALSE], dplyr::any_of(c("AIC", "BIC", "Entropy",
    "WSS", "Silhouette", "CalinskiHarabasz", "Persistence", "NoiseProportion",
    "MeanMembershipProbability", "MinClusterN")))
}

.BenchmarkFailureRow <- function(candidate, error) {
  dplyr::bind_cols(candidate, dplyr::tibble(Eligibility = "failed",
    NonNoiseClusters = NA_integer_, SmallestClusterProportion = NA_real_,
    StabilityARI_Mean = NA_real_, StabilityARI_P05 = NA_real_,
    StabilityJaccard_Mean = NA_real_, StabilitySuccessRate = 0,
    NativeMetrics = list(dplyr::tibble()), ParticipantInclusion = list(dplyr::tibble()),
    FitError = error))
}

.BenchmarkSafeMean <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
.BenchmarkSafeQuantile <- function(x, probability) if (any(is.finite(x))) as.numeric(stats::quantile(x, probability, na.rm = TRUE)) else NA_real_
.BenchmarkMetricValue <- function(metrics, name) {
  if (!is.data.frame(metrics) || !name %in% names(metrics) || !nrow(metrics)) return(NA_real_)
  as.numeric(metrics[[name]][[1]])
}
