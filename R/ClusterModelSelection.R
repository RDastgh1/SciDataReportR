#' Model-selection table for clustering candidates
#'
#' @description One row per candidate model with the criteria used to choose
#' the number of clusters, formatted for reports and slides. It shows the same
#' metrics as [PlotClusterFitReview()], in the same order and labelled with the
#' same preferred directions, together with eligibility, the composite (AHP)
#' recommendation and, for SOM + Mclust, the participant cluster sizes each
#' candidate gives.
#' @param object A fitted model from a `CreateClusterModel_*()` function, or
#'   its candidate fit table.
#' @param metrics Metric columns to show, in column order. When `NULL`, the
#'   pipeline's review metrics are used, as in [PlotClusterFitReview()].
#' @details Numbers are rounded for display; use the fit table
#'   (`ModelInfo_MClust$fit_table` for SOM + Mclust, `ModelInfo$fit_table`
#'   for other pipelines) for exact values.
#'
#'   For SOM + Mclust, `Boundary %` is the share of participants flagged as
#'   boundary participants under soft assignment, and `Cluster sizes` lists the
#'   number of participants in each cluster, largest first, under the model's
#'   `assignment` rule. `Eligible` is false for candidates with fewer SOM nodes
#'   than mixture parameters, or whose smallest cluster is below
#'   `min_cluster_prop`; these cannot be recommended. The `AHP pick` column
#'   marks the candidate recommended by the composite index, which is computed
#'   from z-scores across the candidates in the table (see
#'   [CreateClusterModel_SOM_MClust()]).
#' @return A data frame with one row per candidate, sorted by number of
#'   clusters and then model.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), k_range = 2:4
#' )
#' ClusterModelSelectionTable(model)
#' }
#' @export
ClusterModelSelectionTable <- function(object, metrics = NULL) {
  fit_table <- as.data.frame(.ClusterFitTableFrom(object))
  if (!nrow(fit_table)) stop("The fit table has no candidates.")
  if (is.null(metrics)) metrics <- .ClusterReviewMetrics(fit_table)
  metrics <- metrics[metrics %in% names(fit_table)]
  metrics <- metrics[vapply(fit_table[metrics], function(values) {
    is.numeric(values) && any(is.finite(values))
  }, logical(1))]

  sort_keys <- intersect(c("Classes", "Model", "MinPts", "Epsilon"), names(fit_table))
  if (length(sort_keys)) {
    fit_table <- fit_table[do.call(order, unname(as.list(fit_table[sort_keys]))), , drop = FALSE]
  }
  out <- data.frame(row.names = seq_len(nrow(fit_table)))
  if ("ModelName" %in% names(fit_table)) {
    out$Model <- fit_table$ModelName
  } else if ("Model" %in% names(fit_table)) {
    out$Model <- fit_table$Model
  }
  for (identifier in c("Classes", "MinPts", "Epsilon")) {
    if (identifier %in% names(fit_table)) {
      out[[if (identifier == "Classes") "k" else identifier]] <- fit_table[[identifier]]
    }
  }
  headers <- .ClusterMetricDirectionLabel(metrics, fit_table)
  for (i in seq_along(metrics)) {
    out[[headers[[i]]]] <- .ClusterFormatMetric(fit_table[[metrics[[i]]]], metrics[[i]])
  }
  if ("ClusterSizes" %in% names(fit_table)) {
    out$`Cluster sizes` <- fit_table$ClusterSizes
  }
  if ("Eligible" %in% names(fit_table)) {
    out$Eligible <- ifelse(fit_table$Eligible %in% TRUE, "✓", "✗")
  }
  if ("Recommended" %in% names(fit_table) && nrow(fit_table) > 1L) {
    out$`AHP pick` <- ifelse(fit_table$Recommended %in% TRUE, "◀", "")
  }
  rownames(out) <- NULL
  out
}

# Display formatting for one metric column: shares as percentages, p-values
# with a floor, large information criteria to one decimal, and the rest to two.
.ClusterFormatMetric <- function(values, metric) {
  out <- rep("", length(values))
  finite <- is.finite(values)
  if (!any(finite)) return(out)
  if (identical(metric, "BoundaryPct")) {
    out[finite] <- sprintf("%.1f%%", 100 * values[finite])
  } else if (identical(metric, "BLRTPValue")) {
    out[finite] <- ifelse(values[finite] < 0.001, "<0.001", sprintf("%.3f", values[finite]))
  } else if (all(values[finite] == round(values[finite]))) {
    out[finite] <- format(round(values[finite]), trim = TRUE, big.mark = ",")
  } else if (max(abs(values[finite])) >= 100) {
    out[finite] <- sprintf("%.1f", values[finite])
  } else {
    out[finite] <- sprintf("%.2f", values[finite])
  }
  out
}
