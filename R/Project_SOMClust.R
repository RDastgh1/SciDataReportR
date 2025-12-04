#' Project new data onto an existing SOM + cluster solution
#'
#' @description
#' Given a fitted \code{Pipeline_SOMClust} object, project a new data frame onto:
#' \itemize{
#'   \item The same Z-score scaling (via SciDataReportR::Project_ZScore()).
#'   \item The same SOM (kohonen), using \code{kohonen::map()}.
#'   \item The same node-level latent profiles and posterior probabilities.
#' }
#'
#' For each projected case, the function:
#' \itemize{
#'   \item Assigns a SOM node and distance to BMU.
#'   \item Maps node-level cluster/posterior probabilities to the individual.
#'   \item Computes an overall SOM-distance z-score and flags:
#'     \itemize{
#'       \item \code{Flag_SOMDist_overallHigh} – distance > overall training p95.
#'       \item \code{Flag_SOMDist_clusterHigh} – distance > cluster-specific
#'             training p95 for that cluster.
#'     }
#' }
#'
#' Cluster-level summaries compare the fraction of high-distance projected
#' cases in each cluster to the fraction seen in training. A warning is issued
#' only if \emph{any} projected cluster has a higher proportion of
#' high-distance cases than any cluster in training.
#'
#' Missing data:
#' \itemize{
#'   \item Z-scores are computed for all rows, but only rows with complete
#'         Z-scores are mapped to the SOM.
#'   \item Rows with missing Z-scores receive NA for SOM_Node, SOM_Distance,
#'         and the cluster label.
#' }
#'
#' @param object A \code{Pipeline_SOMClust} object from a previous run.
#' @param new_df Data frame of new cases to project.
#' @param ClusterName Optional name for the cluster column; defaults to
#'   \code{object$ClusterName}. If that column already exists in \code{new_df}
#'   it is overwritten (with a message).
#'
#' @return A list of class \code{"Project_SOMClust"} with components:
#'   \itemize{
#'     \item \code{vars_used}, \code{ClusterName}, \code{complete_rows}
#'     \item \code{df_with_clusters}: \code{new_df} with only the cluster
#'           column appended.
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
Project_SOMClust <- function(
    object,
    new_df,
    ClusterName = NULL
) {

  if (!inherits(object, "Pipeline_SOMClust")) {
    stop("object must be a Pipeline_SOMClust result.")
  }

  if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
    stop("SciDataReportR must be installed.")
  }
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr' and 'ggplot2' are required.")
  }

  vars_used <- object$vars_used
  if (is.null(ClusterName)) {
    ClusterName <- object$ClusterName
  }

  # Z-score projection onto new_df -----------------------------------------

  Z_obj <- object$ZScoreObject
  if (is.null(Z_obj) || !"ZScoreObj" %in% class(Z_obj)) {
    stop("object$ZScoreObject must be a valid ZScoreObj from Pipeline_SOMClust.")
  }

  z_proj <- SciDataReportR::Project_ZScore(
    df                 = new_df,
    variables          = vars_used,
    parameters         = Z_obj,
    ParameterInputType = "ZScoreObj",
    names_prefix       = "Z_"
  )

  z_df_new <- z_proj$ZScores
  complete_rows <- stats::complete.cases(z_df_new)
  if (!any(complete_rows)) {
    stop("No complete rows in new_df after Z-score projection; cannot map to SOM.")
  }

  zmat_new <- as.matrix(z_df_new[complete_rows, , drop = FALSE])

  # Retrieve SOM / distance baselines from training ------------------------

  som_model <- object$ModelInfo_SOM$som_model

  SOMFit_train        <- object$ModelInfo_SOM$SOMFit
  dist_summary_train  <- SOMFit_train$distance_summary
  overall_p95         <- SOMFit_train$overall_p95
  overall_mean        <- SOMFit_train$overall_mean
  overall_sd          <- SOMFit_train$overall_sd
  dist_by_cluster_tr  <- SOMFit_train$dist_by_cluster
  flag_by_cluster_tr  <- SOMFit_train$flag_by_cluster

  # Map new data to SOM ----------------------------------------------------

  mapping <- kohonen::map(som_model, newdata = zmat_new)

  SOM_Node_new <- rep(NA_integer_, nrow(new_df))
  SOM_Node_new[complete_rows] <- mapping$unit.classif

  SOM_Dist_new <- rep(NA_real_, nrow(new_df))
  SOM_Dist_new[complete_rows] <- mapping$distances

  # Node-level info from training -----------------------------------------

  node_df <- object$ProbFit$node
  if (is.null(node_df)) {
    stop("object$ProbFit$node is missing; run Pipeline_SOMClust first.")
  }

  # Individual table with distance flags -----------------------------------

  individual_tbl <- dplyr::tibble(
    RowID        = seq_len(nrow(new_df)),
    SOM_Node     = SOM_Node_new,
    SOM_Distance = SOM_Dist_new
  ) %>%
    dplyr::left_join(node_df, by = c("SOM_Node" = "NodeID")) %>%
    dplyr::mutate(
      SOMDist_z_overall = if (is.finite(overall_sd) && overall_sd > 0)
        (SOM_Distance - overall_mean) / overall_sd else NA_real_,
      Flag_SOMDist_overallHigh =
        !is.na(SOM_Distance) & !is.na(overall_p95) & SOM_Distance > overall_p95
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
    dplyr::mutate(
      Flag_SOMDist_clusterHigh =
        !is.na(SOM_Distance) & !is.na(train_p95) & SOM_Distance > train_p95
    )

  # df_with_clusters for projected data ------------------------------------

  df_with_clusters <- new_df

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
          n_train                = n,
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

  SOMProj <- list(
    distance_summary_train   = dist_summary_train,
    distance_summary_proj    = dist_proj_summary,
    flag_by_cluster_train    = flag_by_cluster_tr,
    flag_by_cluster_proj     = flag_by_cluster_proj,
    distance_flag_comparison = distance_flag_comparison,
    table                    = somproj_tbl,
    plots                    = list(
      distance_hist           = p_dist_hist,
      distance_by_cluster_box = p_dist_box
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
    ProbFit          = ProbFit,
    ModelInfo_SOM    = object$ModelInfo_SOM,
    ModelInfo_MClust = object$ModelInfo_MClust
  )

  class(out) <- c("Project_SOMClust", class(out))
  out
}
