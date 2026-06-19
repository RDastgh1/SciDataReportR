
#' Project new data onto an existing SOM clinical phenotype space
#'
#' @description
#' Train once, project many: a reusable unsupervised clinical phenotyping
#' framework that learns phenotype structure in a training cohort and projects
#' new participants into the fixed phenotype space while quantifying membership
#' uncertainty and projection fit. Given a fitted \code{CreateSOMClusterModel()}
#' object, project a new data frame onto:
#' \itemize{
#'   \item The same Z-score scaling (via SciDataReportR::ProjectZScore()) when
#'         the training object used computed or projected Z-scores.
#'   \item The same SOM (kohonen), using \code{kohonen::map()}.
#'   \item The same node-level latent profiles and posterior probabilities.
#' }
#'
#' If the training object used \code{ZScoreType = "PreZScored"}, this function
#' expects the new data already contains the same Z-score columns (by name)
#' stored in \code{object$ZScoreVars}. No re-zscoring is performed.
#'
#' For each projected case, the function:
#' \itemize{
#'   \item Assigns a SOM node and distance to BMU.
#'   \item Maps node-level cluster/posterior probabilities to the individual.
#'   \item Computes an overall SOM-distance z-score and flags:
#'     \itemize{
#'       \item \code{Flag_SOMDist_overallHigh} - distance > overall training p95.
#'       \item \code{Flag_SOMDist_clusterHigh} - distance > cluster-specific
#'             training p95 for that cluster.
#'     }
#' }
#'
#' Stable row id:
#' - \code{.scidr_rowid} is added to \code{new_df} (if not already present) and
#'   carried into \code{df_with_clusters} and \code{ProbFit$individual}.
#' - \code{ProbFit$individual$RowID} is set equal to \code{.scidr_rowid}.
#'
#' Missing data:
#' \itemize{
#'   \item Only rows with complete Z-scores are mapped to the SOM.
#'   \item Rows with missing Z-scores receive NA for SOM_Node, SOM_Distance,
#'         and the cluster label.
#' }
#'
#' @param object A SOM cluster model object from [CreateSOMClusterModel()].
#' @param new_df Data frame of new cases to project.
#' @param ClusterName Optional name for the cluster column; defaults to
#'   \code{object$ClusterName}. If that column already exists in \code{new_df}
#'   it is overwritten (with a message).
#' @param high_dist_quantile Numeric value between 0 and 1 used to define
#'   high SOM-distance flags from the training distance distribution. Default
#'   is \code{0.95}.
#' @param low_prob_threshold Numeric posterior probability threshold used to
#'   flag uncertain phenotype membership. Default is \code{0.70}.
#'
#' @return A list of class \code{"Project_SOMClust"} with components:
#'   \itemize{
#'     \item \code{vars_used}, \code{ClusterName}, \code{complete_rows}
#'     \item \code{df_with_clusters}: \code{new_df} with \code{.scidr_rowid} and
#'           only the cluster column appended.
#'     \item \code{SOMProj}: list with training and projected distance
#'           summaries, cluster-level flag summaries, comparison, and plots.
#'     \item \code{ProbFit}: list with \code{node} (training node-level
#'           posterior info), \code{individual} (projection-level info
#'           including distance flags and z-scores), and probability plots.
#'     \item \code{ModelInfo_SOM}, \code{ModelInfo_MClust}: references to the
#'           original model objects for convenience.
#'   }
#'
#' @export
ProjectSOMCluster <- function(
    object,
    new_df,
    ClusterName = NULL,
    high_dist_quantile = 0.95,
    low_prob_threshold = 0.70
) {

  if (!inherits(object, "Pipeline_SOMClust")) {
    stop("object must be a Pipeline_SOMClust result.")
  }

  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr' and 'ggplot2' are required.")
  }


  if (!is.numeric(high_dist_quantile) || length(high_dist_quantile) != 1 ||
      is.na(high_dist_quantile) || high_dist_quantile <= 0 || high_dist_quantile >= 1) {
    stop("high_dist_quantile must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(low_prob_threshold) || length(low_prob_threshold) != 1 ||
      is.na(low_prob_threshold) || low_prob_threshold < 0 || low_prob_threshold > 1) {
    stop("low_prob_threshold must be a single numeric value between 0 and 1.")
  }

  vars_used <- object$vars_used
  if (is.null(ClusterName)) {
    ClusterName <- object$ClusterName
  }


  missing_vars <- setdiff(vars_used, names(new_df))
  if (length(missing_vars) > 0) {
    stop(
      "new_df is missing required original training variable(s): ",
      paste(missing_vars, collapse = ", "),
      ". These variables are required so projected cohorts can be scaled and checked against the training range."
    )
  }

  training_variable_summary <- object$ModelInfo_SOM$PhenotypeReference$training_variable_summary
  if (is.null(training_variable_summary)) {
    training_variable_summary <- object$ModelInfo_SOM$training_variable_summary
  }
  if (!is.null(training_variable_summary)) {
    out_of_range_summary <- training_variable_summary %>%
      dplyr::filter(.data$Variable %in% names(new_df)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n_below_training_min = if (!is.numeric(new_df[[Variable]]) || is.na(min)) NA_integer_ else sum(new_df[[Variable]] < min, na.rm = TRUE),
        n_above_training_max = if (!is.numeric(new_df[[Variable]]) || is.na(max)) NA_integer_ else sum(new_df[[Variable]] > max, na.rm = TRUE),
        n_out_of_range = if (is.na(n_below_training_min) || is.na(n_above_training_max)) NA_integer_ else n_below_training_min + n_above_training_max,
        prop_out_of_range = if (is.na(n_out_of_range)) NA_real_ else n_out_of_range / sum(!is.na(new_df[[Variable]])),
        projected_min = if (!is.numeric(new_df[[Variable]]) || all(is.na(new_df[[Variable]]))) NA_real_ else min(new_df[[Variable]], na.rm = TRUE),
        projected_max = if (!is.numeric(new_df[[Variable]]) || all(is.na(new_df[[Variable]]))) NA_real_ else max(new_df[[Variable]], na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  } else {
    out_of_range_summary <- dplyr::tibble(
      Variable = vars_used,
      n_below_training_min = NA_integer_,
      n_above_training_max = NA_integer_,
      n_out_of_range = NA_integer_,
      prop_out_of_range = NA_real_,
      projected_min = NA_real_,
      projected_max = NA_real_
    )
  }

  out_of_range_vars <- out_of_range_summary %>%
    dplyr::filter(!is.na(.data$n_out_of_range), .data$n_out_of_range > 0) %>%
    dplyr::pull(.data$Variable)

  if (length(out_of_range_vars) > 0) {
    warning(
      "Projected data contain values outside the training range for: ",
      paste(out_of_range_vars, collapse = ", "),
      ". Review SOMProj$out_of_range_summary."
    )
  }

  # Stable row id ----------------------------------------------------------

  if (!".scidr_rowid" %in% names(new_df)) {
    new_df_scidr <- new_df %>%
      dplyr::mutate(.scidr_rowid = dplyr::row_number())
  } else {
    new_df_scidr <- new_df
    if (any(is.na(new_df_scidr$.scidr_rowid))) {
      stop("new_df has .scidr_rowid but it contains missing values.")
    }
  }

  # Determine Z-score columns used in training -----------------------------

  ZScoreType_train <- object$ZScoreType
  ZScoreVars_used  <- object$ZScoreVars

  if (is.null(ZScoreVars_used)) {
    ZScoreVars_used <- paste0("Z_", vars_used)
  }

  # Z-score projection onto new_df -----------------------------------------

  if (!is.null(ZScoreType_train) && ZScoreType_train == "PreZScored") {

    missing_z <- setdiff(ZScoreVars_used, names(new_df_scidr))
    if (length(missing_z) > 0) {
      stop("new_df is missing required pre-zscored columns: ",
           paste(missing_z, collapse = ", "))
    }

    z_df_new <- new_df_scidr[, ZScoreVars_used, drop = FALSE]

    is_num_z <- vapply(z_df_new, is.numeric, logical(1))
    if (!all(is_num_z)) {
      stop("All PreZScored columns must be numeric. Non-numeric: ",
           paste(names(z_df_new)[!is_num_z], collapse = ", "))
    }

    complete_rows <- stats::complete.cases(z_df_new)
    if (!any(complete_rows)) {
      stop("No complete rows in new_df pre-zscored columns; cannot map to SOM.")
    }

    if (!identical(names(z_df_new), ZScoreVars_used)) {
      stop("Projected Z-score variable ordering does not match the training object.")
    }

    if (!identical(names(z_df_new), ZScoreVars_used)) {
      stop("Projected Z-score variable ordering does not match the training object.")
    }

    zmat_new <- as.matrix(z_df_new[complete_rows, , drop = FALSE])

  } else {

    if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
      stop("SciDataReportR must be installed.")
    }

    Z_obj <- object$ZScoreObject
    if (is.null(Z_obj) || !"ZScoreObj" %in% class(Z_obj)) {
      stop("object$ZScoreObject must be a valid ZScoreObj from CreateSOMClusterModel.")
    }

    z_proj <- SciDataReportR::ProjectZScore(
      df                 = new_df_scidr,
      variables          = vars_used,
      parameters         = Z_obj,
      ParameterInputType = "ZScoreObj",
      names_prefix       = "Z_"
    )

    z_df_new <- z_proj$ZScores

    missing_z <- setdiff(ZScoreVars_used, names(z_df_new))
    if (length(missing_z) > 0) {
      stop("Projected Z-scores are missing required columns: ",
           paste(missing_z, collapse = ", "))
    }

    z_df_new <- z_df_new[, ZScoreVars_used, drop = FALSE]

    complete_rows <- stats::complete.cases(z_df_new)
    if (!any(complete_rows)) {
      stop("No complete rows in new_df after Z-score projection; cannot map to SOM.")
    }

    zmat_new <- as.matrix(z_df_new[complete_rows, , drop = FALSE])
  }

  # Retrieve SOM / distance baselines from training ------------------------

  som_model <- object$ModelInfo_SOM$som_model

  SOMFit_train        <- object$ModelInfo_SOM$SOMFit
  dist_summary_train  <- SOMFit_train$distance_summary
  overall_mean        <- SOMFit_train$overall_mean
  overall_sd          <- SOMFit_train$overall_sd
  dist_by_cluster_tr  <- SOMFit_train$dist_by_cluster
  flag_by_cluster_tr  <- SOMFit_train$flag_by_cluster

  training_individual <- object$ProbFit$individual
  train_dist <- training_individual$SOM_Distance[!is.na(training_individual$SOM_Distance)]
  overall_high_cutoff <- as.numeric(stats::quantile(train_dist, high_dist_quantile, na.rm = TRUE))

  dist_by_cluster_cutoffs <- training_individual %>%
    dplyr::filter(!is.na(.data$SOM_Distance), !is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      train_cluster_high_cutoff = as.numeric(stats::quantile(.data$SOM_Distance, high_dist_quantile, na.rm = TRUE)),
      .groups = "drop"
    )

  # Map new data to SOM ----------------------------------------------------

  mapping <- kohonen::map(som_model, newdata = zmat_new)

  SOM_Node_new <- rep(NA_integer_, nrow(new_df_scidr))
  SOM_Node_new[complete_rows] <- mapping$unit.classif

  SOM_Dist_new <- rep(NA_real_, nrow(new_df_scidr))
  SOM_Dist_new[complete_rows] <- mapping$distances

  # Node-level info from training -----------------------------------------

  node_df <- object$ProbFit$node
  if (is.null(node_df)) {
    stop("object$ProbFit$node is missing; run CreateSOMClusterModel first.")
  }

  # Individual table with distance flags -----------------------------------

  individual_tbl <- dplyr::tibble(
    .scidr_rowid = new_df_scidr$.scidr_rowid,
    RowID        = new_df_scidr$.scidr_rowid,
    SOM_Node     = SOM_Node_new,
    SOM_Distance = SOM_Dist_new
  ) %>%
    dplyr::left_join(node_df, by = c("SOM_Node" = "NodeID")) %>%
    dplyr::mutate(
      SOMDist_z_overall = if (is.finite(overall_sd) && overall_sd > 0)
        (SOM_Distance - overall_mean) / overall_sd else NA_real_,
      SOMDist_percentile_train = vapply(SOM_Distance, function(x) {
        if (is.na(x) || length(train_dist) == 0) NA_real_ else mean(train_dist <= x, na.rm = TRUE)
      }, numeric(1)),
      Flag_SOMDist_overallHigh =
        !is.na(SOM_Distance) & !is.na(overall_high_cutoff) & SOM_Distance > overall_high_cutoff
    ) %>%
    dplyr::left_join(
      dist_by_cluster_tr %>%
        dplyr::select(
          Cluster,
          train_mean_dist = mean,
          train_sd_dist   = sd,
          train_p90       = p90,
          train_p95       = p95,
          train_p99       = p99
        ),
      by = "Cluster"
    ) %>%
    dplyr::left_join(dist_by_cluster_cutoffs, by = "Cluster") %>%
    dplyr::mutate(
      Flag_SOMDist_clusterHigh =
        !is.na(SOM_Distance) & !is.na(train_cluster_high_cutoff) & SOM_Distance > train_cluster_high_cutoff,
      Projection_Fit_Class = dplyr::case_when(
        is.na(SOM_Distance) ~ NA_character_,
        (Flag_SOMDist_overallHigh | Flag_SOMDist_clusterHigh) &
          (is.na(max_prob) | max_prob < low_prob_threshold) ~ "Potential novel phenotype",
        Flag_SOMDist_overallHigh | Flag_SOMDist_clusterHigh ~ "Poor SOM fit",
        !is.na(max_prob) & max_prob < low_prob_threshold ~ "Uncertain membership",
        TRUE ~ "Good fit"
      )
    )

  # df_with_clusters for projected data ------------------------------------

  df_with_clusters <- new_df_scidr

  if (ClusterName %in% names(df_with_clusters)) {
    message("Column '", ClusterName, "' already exists in new_df and will be overwritten.")
  }

  df_with_clusters[[ClusterName]] <- individual_tbl$Cluster

  # SOMProj / distance diagnostics for projected cases ---------------------

  somproj_tbl    <- individual_tbl
  somproj_non_na <- somproj_tbl[!is.na(somproj_tbl$SOM_Distance), ]

  dist_proj_summary <- stats::quantile(
    somproj_non_na$SOM_Distance,
    probs = c(0, 0.25, 0.5, 0.75, 0.95, 1),
    na.rm = TRUE
  )
  names(dist_proj_summary) <- c("min", "q25", "median", "q75", "p95", "max")

  flag_by_cluster_proj <- somproj_non_na %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      n_proj            = dplyr::n(),
      prop_overall_high = mean(Flag_SOMDist_overallHigh),
      prop_cluster_high = mean(Flag_SOMDist_clusterHigh),
      .groups = "drop"
    )

  distance_flag_comparison <- flag_by_cluster_proj %>%
    dplyr::left_join(
      flag_by_cluster_tr %>%
        dplyr::rename(
          n_train                 = n,
          train_prop_overall_high = prop_overall_high
        ),
      by = "Cluster"
    ) %>%
    dplyr::mutate(
      diff_prop_overall_high =
        prop_overall_high - train_prop_overall_high
    )

  max_train_prop <- max(flag_by_cluster_tr$prop_overall_high, na.rm = TRUE)
  max_proj_prop  <- max(flag_by_cluster_proj$prop_overall_high, na.rm = TRUE)

  if (is.finite(max_proj_prop) && max_proj_prop > max_train_prop) {
    warning(
      "At least one projected cluster has a higher proportion of high-distance ",
      "cases than any cluster in training (proj max = ",
      sprintf("%.1f%%", 100 * max_proj_prop),
      ", train max = ",
      sprintf("%.1f%%", 100 * max_train_prop),
      "). Review flags Flag_SOMDist_overallHigh / Flag_SOMDist_clusterHigh ",
      "in ProbFit$individual and SOMProj$distance_flag_comparison."
    )
  }


  train_node_occupancy <- SOMFit_train$node_occupancy %>%
    dplyr::mutate(cohort = "Training")

  proj_node_occupancy <- dplyr::tibble(
    NodeID = seq_len(nrow(object$ModelInfo_SOM$som_codes))
  ) %>%
    dplyr::left_join(
      as.data.frame(table(somproj_non_na$SOM_Node)) %>%
        dplyr::transmute(
          NodeID = as.integer(as.character(Var1)),
          n = as.integer(Freq)
        ),
      by = "NodeID"
    ) %>%
    dplyr::mutate(
      n = dplyr::if_else(is.na(.data$n), 0L, .data$n),
      cohort = "Projected"
    )

  train_cluster_occupancy <- training_individual %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::count(.data$Cluster, name = "n") %>%
    dplyr::mutate(cohort = "Training")

  proj_cluster_occupancy <- somproj_non_na %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::count(.data$Cluster, name = "n") %>%
    dplyr::mutate(cohort = "Projected")

  all_nodes <- sort(unique(c(train_node_occupancy$NodeID, proj_node_occupancy$NodeID)))
  train_node_p <- train_node_occupancy %>%
    dplyr::right_join(dplyr::tibble(NodeID = all_nodes), by = "NodeID") %>%
    dplyr::mutate(n = dplyr::if_else(is.na(.data$n), 0L, .data$n)) %>%
    dplyr::arrange(.data$NodeID) %>%
    dplyr::pull(.data$n)
  proj_node_p <- proj_node_occupancy %>%
    dplyr::right_join(dplyr::tibble(NodeID = all_nodes), by = "NodeID") %>%
    dplyr::mutate(n = dplyr::if_else(is.na(.data$n), 0L, .data$n)) %>%
    dplyr::arrange(.data$NodeID) %>%
    dplyr::pull(.data$n)
  train_node_p <- train_node_p / sum(train_node_p)
  proj_node_p <- proj_node_p / sum(proj_node_p)
  node_m <- 0.5 * (train_node_p + proj_node_p)
  node_js <- 0.5 * sum(ifelse(train_node_p > 0, train_node_p * log2(train_node_p / node_m), 0), na.rm = TRUE) +
    0.5 * sum(ifelse(proj_node_p > 0, proj_node_p * log2(proj_node_p / node_m), 0), na.rm = TRUE)

  all_clusters <- sort(unique(c(train_cluster_occupancy$Cluster, proj_cluster_occupancy$Cluster)))
  train_cluster_p <- train_cluster_occupancy %>%
    dplyr::right_join(dplyr::tibble(Cluster = all_clusters), by = "Cluster") %>%
    dplyr::mutate(n = dplyr::if_else(is.na(.data$n), 0L, .data$n)) %>%
    dplyr::arrange(.data$Cluster) %>%
    dplyr::pull(.data$n)
  proj_cluster_p <- proj_cluster_occupancy %>%
    dplyr::right_join(dplyr::tibble(Cluster = all_clusters), by = "Cluster") %>%
    dplyr::mutate(n = dplyr::if_else(is.na(.data$n), 0L, .data$n)) %>%
    dplyr::arrange(.data$Cluster) %>%
    dplyr::pull(.data$n)
  train_cluster_p <- train_cluster_p / sum(train_cluster_p)
  proj_cluster_p <- proj_cluster_p / sum(proj_cluster_p)
  cluster_m <- 0.5 * (train_cluster_p + proj_cluster_p)
  cluster_js <- 0.5 * sum(ifelse(train_cluster_p > 0, train_cluster_p * log2(train_cluster_p / cluster_m), 0), na.rm = TRUE) +
    0.5 * sum(ifelse(proj_cluster_p > 0, proj_cluster_p * log2(proj_cluster_p / cluster_m), 0), na.rm = TRUE)

  training_mean_distance <- mean(train_dist, na.rm = TRUE)
  projected_mean_distance <- mean(somproj_non_na$SOM_Distance, na.rm = TRUE)
  training_median_distance <- stats::median(train_dist, na.rm = TRUE)
  projected_median_distance <- stats::median(somproj_non_na$SOM_Distance, na.rm = TRUE)

  distance_ratio <- projected_mean_distance / training_mean_distance
  phenotype_drift_index <- abs(log(distance_ratio))

  ProjectionDiagnostics <- dplyr::tibble(
    metric = c(
      "n_total",
      "n_complete",
      "n_missing",
      "training_mean_som_distance",
      "projected_mean_som_distance",
      "training_median_som_distance",
      "projected_median_som_distance",
      "mean_distance_ratio",
      "high_distance_burden",
      "cluster_high_distance_burden",
      "Phenotype Drift Index"
    ),
    value = c(
      nrow(new_df_scidr),
      sum(complete_rows),
      sum(!complete_rows),
      training_mean_distance,
      projected_mean_distance,
      training_median_distance,
      projected_median_distance,
      distance_ratio,
      mean(somproj_non_na$Flag_SOMDist_overallHigh, na.rm = TRUE),
      mean(somproj_non_na$Flag_SOMDist_clusterHigh, na.rm = TRUE),
      phenotype_drift_index
    )
  )

  TransportabilityDiagnostics <- list(
    summary = dplyr::tibble(
      metric = c(
        "js_divergence_node_occupancy",
        "js_divergence_cluster_occupancy"
      ),
      value = c(
        node_js,
        cluster_js
      )
    ),
    node_occupancy_train = train_node_occupancy,
    node_occupancy_proj = proj_node_occupancy,
    cluster_occupancy_train = train_cluster_occupancy,
    cluster_occupancy_proj = proj_cluster_occupancy
  )

  distance_compare_tbl <- dplyr::bind_rows(
    training_individual %>%
      dplyr::filter(!is.na(.data$SOM_Distance)) %>%
      dplyr::transmute(cohort = "Training", SOM_Distance = .data$SOM_Distance, Cluster = .data$Cluster),
    somproj_non_na %>%
      dplyr::transmute(cohort = "Projected", SOM_Distance = .data$SOM_Distance, Cluster = .data$Cluster)
  )

  cluster_occupancy_compare <- dplyr::bind_rows(train_cluster_occupancy, proj_cluster_occupancy) %>%
    dplyr::group_by(.data$cohort) %>%
    dplyr::mutate(prop = .data$n / sum(.data$n)) %>%
    dplyr::ungroup()

  p_dist_hist <- ggplot2::ggplot(
    somproj_non_na,
    ggplot2::aes(x = SOM_Distance)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "SOM distances (projected cases)",
      x     = "Distance to BMU",
      y     = "Count"
    )

  p_dist_box <- ggplot2::ggplot(
    somproj_non_na,
    ggplot2::aes(x = factor(Cluster), y = SOM_Distance)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "SOM distances by cluster (projected cases)",
      x     = "Cluster",
      y     = "Distance to BMU"
    )


  p_train_proj_dist_density <- ggplot2::ggplot(
    distance_compare_tbl,
    ggplot2::aes(x = SOM_Distance, fill = cohort)
  ) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Training vs projected SOM distance",
      x = "Distance to BMU",
      y = "Density",
      fill = "Cohort"
    )

  p_train_proj_dist_ecdf <- ggplot2::ggplot(
    distance_compare_tbl,
    ggplot2::aes(x = .data$SOM_Distance, colour = .data$cohort)
  ) +
    ggplot2::stat_ecdf() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Training vs projected SOM distance ECDF",
      x = "Distance to BMU",
      y = "Empirical cumulative probability",
      colour = "Cohort"
    )

  qq_probs <- seq(0.01, 0.99, by = 0.01)
  qq_df <- dplyr::tibble(
    train_quantile = as.numeric(stats::quantile(train_dist, probs = qq_probs, na.rm = TRUE)),
    projected_quantile = as.numeric(stats::quantile(somproj_non_na$SOM_Distance, probs = qq_probs, na.rm = TRUE))
  )

  p_train_proj_dist_qq <- ggplot2::ggplot(
    qq_df,
    ggplot2::aes(x = .data$train_quantile, y = .data$projected_quantile)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Training vs projected SOM distance QQ plot",
      x = "Training SOM distance quantile",
      y = "Projected SOM distance quantile"
    )

  p_cluster_occupancy <- ggplot2::ggplot(
    cluster_occupancy_compare,
    ggplot2::aes(x = factor(Cluster), y = prop, fill = cohort)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Training vs projected cluster occupancy",
      x = "Cluster",
      y = "Proportion",
      fill = "Cohort"
    )

  somproj_plot_non_na <- somproj_non_na %>%
    dplyr::filter(!is.na(.data$Cluster), !is.na(.data$max_prob))

  p_projected_posterior_by_cluster <- ggplot2::ggplot(
    somproj_plot_non_na,
    ggplot2::aes(x = .data$max_prob)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Projected posterior probability by assigned cluster",
      x = "Max posterior probability",
      y = "Count"
    )

  cluster_fit_summary <- somproj_non_na %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data$SOM_Distance, na.rm = TRUE),
      median_distance = stats::median(.data$SOM_Distance, na.rm = TRUE),
      mean_probability = mean(.data$prob_assigned, na.rm = TRUE),
      median_probability = stats::median(.data$prob_assigned, na.rm = TRUE),
      mean_distance_percentile = mean(.data$SOMDist_percentile_train, na.rm = TRUE),
      median_distance_percentile = stats::median(.data$SOMDist_percentile_train, na.rm = TRUE),
      prop_high_distance = mean(.data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh, na.rm = TRUE),
      prop_low_probability = mean(.data$prob_assigned < low_prob_threshold, na.rm = TRUE),
      prop_poor_fit = mean(.data$Projection_Fit_Class == "Poor SOM fit", na.rm = TRUE),
      prop_potential_novel = mean(.data$Projection_Fit_Class == "Potential novel phenotype", na.rm = TRUE),
      .groups = "drop"
    )

  poor_fit_by_cluster <- cluster_fit_summary %>%
    dplyr::transmute(
      Cluster = .data$Cluster,
      n = .data$n,
      prop_poor_or_high_dist = .data$prop_high_distance,
      prop_potential_novel = .data$prop_potential_novel
    )

  p_projection_fit_class <- ggplot2::ggplot(
    somproj_non_na %>% dplyr::filter(!is.na(.data$Projection_Fit_Class)),
    ggplot2::aes(x = .data$Projection_Fit_Class)
  ) +
    ggplot2::geom_bar() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Projected cases by projection fit class",
      x = "Projection fit class",
      y = "Count"
    )

  p_poor_fit_by_cluster <- ggplot2::ggplot(
    poor_fit_by_cluster,
    ggplot2::aes(x = factor(Cluster), y = prop_poor_or_high_dist)
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Poor-fit or high-distance projected cases by cluster",
      x = "Cluster",
      y = "Proportion"
    )

  SOMProj <- list(
    distance_summary_train   = dist_summary_train,
    distance_summary_proj    = dist_proj_summary,
    flag_by_cluster_train    = flag_by_cluster_tr,
    flag_by_cluster_proj     = flag_by_cluster_proj,
    distance_flag_comparison = distance_flag_comparison,
    ProjectionDiagnostics    = ProjectionDiagnostics,
    TransportabilityDiagnostics = TransportabilityDiagnostics,
    out_of_range_summary     = out_of_range_summary,
    node_occupancy_train     = train_node_occupancy,
    node_occupancy_proj      = proj_node_occupancy,
    cluster_occupancy_train  = train_cluster_occupancy,
    cluster_occupancy_proj   = proj_cluster_occupancy,
    cluster_occupancy_compare = cluster_occupancy_compare,
    poor_fit_by_cluster      = poor_fit_by_cluster,
    cluster_fit_summary      = cluster_fit_summary,
    table                    = somproj_tbl,
    plots                    = list(
      distance_hist           = p_dist_hist,
      distance_by_cluster_box = p_dist_box,
      training_vs_projected_distance_density = p_train_proj_dist_density,
      training_vs_projected_distance_ecdf = p_train_proj_dist_ecdf,
      training_vs_projected_distance_qq = p_train_proj_dist_qq,
      training_vs_projected_cluster_occupancy = p_cluster_occupancy,
      projected_posterior_by_cluster = p_projected_posterior_by_cluster,
      projection_fit_class_bar = p_projection_fit_class,
      poor_fit_by_cluster = p_poor_fit_by_cluster
    )
  )

  # ProbFit for projected cases --------------------------------------------

  indiv_non_na <- individual_tbl[!is.na(individual_tbl$max_prob), ]
  indiv_non_na$Cluster <- factor(indiv_non_na$Cluster)

  individual_MaxProbBoxplot <- ggplot2::ggplot(
    indiv_non_na,
    ggplot2::aes(x = Cluster, y = max_prob)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Individual max posterior probability (projected cases)",
      x     = "Cluster",
      y     = "Max posterior probability"
    )

  individual_ProbAssignedDensity <- ggplot2::ggplot(
    indiv_non_na,
    ggplot2::aes(x = prob_assigned)
  ) +
    ggplot2::geom_density(fill = "darkblue", alpha = 0.7) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Density of prob_assigned by class (projected cases)",
      x     = "Posterior probability for assigned class",
      y     = "Density"
    )

  ProbFit <- list(
    node       = node_df,
    individual = individual_tbl,
    plots      = list(
      individual_MaxProbBoxplot      = individual_MaxProbBoxplot,
      individual_ProbAssignedDensity = individual_ProbAssignedDensity
    )
  )

  out <- list(
    vars_used        = vars_used,
    ClusterName      = ClusterName,
    complete_rows    = complete_rows,
    df_with_clusters = df_with_clusters,
    SOMProj          = SOMProj,
    out_of_range_summary = out_of_range_summary,
    ProbFit          = ProbFit,
    ModelInfo_SOM    = object$ModelInfo_SOM,
    ModelInfo_MClust = object$ModelInfo_MClust
  )

  class(out) <- c("Project_SOMClust", class(out))
  out
}

#' Project new data onto an existing SOM clinical phenotype space
#'
#' Compatibility alias for [ProjectSOMCluster()]. Prefer `ProjectSOMCluster()`
#' in new code.
#'
#' @param ... Arguments passed to [ProjectSOMCluster()].
#' @return The same projection object returned by [ProjectSOMCluster()].
#' @export
Project_SOMClust <- function(...) {
  ProjectSOMCluster(...)
}
