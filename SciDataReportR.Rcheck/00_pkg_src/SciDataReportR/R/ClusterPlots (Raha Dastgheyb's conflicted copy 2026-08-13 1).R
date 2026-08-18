#' Shared clustering figures
#'
#' @description Figure builders shared by every clustering pipeline in this
#' package. Each `Create*ClusterModel()` and `Project*Cluster()` function stores
#' its figures beside the object they describe, following the layout established
#' by [CreateClusterModel_SOM_MClust()]:
#'
#' \itemize{
#'   \item `fit_plot` reviews candidate solutions across the search grid.
#'   \item `ModelInfo$plots` describes the structure of the selected solution.
#'   \item `ModelInfo$FitDiagnostics$plots` describes how well individual
#'     training cases sit inside that structure.
#'   \item `ProbFit$plots` describes membership confidence.
#'   \item `Stability$plots` describes bootstrap reproducibility.
#'   \item `ProjectionFit$plots` describes how projected cases compare with the
#'     frozen training reference.
#' }
#'
#' There is deliberately no flat `plots` list: a figure lives next to the table
#' it summarises so that the object stays self-describing.
#' @keywords internal
#' @name ClusterPlots
NULL

.ClusterVariableLabel <- function(data, variable, codebook = NULL) {
  if (!is.null(codebook) && all(c("Variable", "Label") %in% names(codebook))) {
    match_row <- codebook$Label[codebook$Variable == variable]
    if (length(match_row) && !is.na(match_row[[1]]) && nzchar(match_row[[1]])) {
      return(as.character(match_row[[1]]))
    }
  }
  if (!variable %in% names(data)) return(variable)
  variable_label <- attr(data[[variable]], "label", exact = TRUE)
  if (is.null(variable_label) || !length(variable_label) ||
      is.na(variable_label[[1]]) || !nzchar(as.character(variable_label[[1]]))) {
    variable
  } else {
    as.character(variable_label[[1]])
  }
}

.ClusterMetricLabel <- function(values) {
  labels <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", values)
  labels <- gsub("_", " ", labels)
  gsub("A[hH][pP]", "AHP", labels)
}

# Clusters are integers internally, with 0 reserved for the density-based noise
# group. Display always separates noise so a reader never mistakes it for a
# phenotype.
.ClusterFactor <- function(cluster, noise_label = 0L) {
  values <- as.character(cluster)
  levels_present <- sort(unique(values[!is.na(values)]))
  if (!is.null(noise_label) && as.character(noise_label) %in% levels_present) {
    levels_present <- c(
      setdiff(levels_present, as.character(noise_label)),
      as.character(noise_label))
  }
  labels <- levels_present
  if (!is.null(noise_label)) {
    labels[levels_present == as.character(noise_label)] <- "Noise"
  }
  factor(values, levels = levels_present, labels = labels)
}

.ClusterOccupancyTable <- function(cluster, noise_label = 0L) {
  df_Cluster <- dplyr::tibble(Cluster = .ClusterFactor(cluster, noise_label))
  df_Cluster %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::count(.data$Cluster, name = "n") %>%
    dplyr::mutate(
      Proportion = .data$n / sum(.data$n),
      Label = paste0(.data$n, " (", sprintf("%.0f%%", 100 * .data$Proportion), ")"))
}

# One colour assignment for every cluster figure, so a cluster keeps the same
# colour across the occupancy, silhouette, map, and composition views. A
# density-based noise group is held out in grey rather than taking a cluster
# colour of its own.
#
# Noise is matched by prefix, not equality: PlotClusterMap() rewrites its
# levels to "Noise (n = 12, 10%)" before the scale sees them, and an equality
# test would silently give noise a real cluster colour and shift every other
# cluster by one.
.ClusterPalette <- function(cluster) {
  cluster_levels <- levels(as.factor(cluster))
  is_noise <- startsWith(cluster_levels, "Noise")
  fixed <- stats::setNames(
    rep("grey60", sum(is_noise)), cluster_levels[is_noise])
  .SciDataNamedValues(cluster_levels, fixed = if (length(fixed)) fixed)
}

# Passing a column that does not exist yields NULL, and an empty palette then
# fails much later inside ggplot with "Insufficient values in manual scale",
# which points at the scale rather than at the caller. Fail here instead.
.ClusterScaleValues <- function(cluster, call_name) {
  values <- .ClusterPalette(cluster)
  if (!length(values)) {
    stop(
      call_name, "() received no cluster levels. This usually means the ",
      "cluster column named in the call does not exist on that data frame.",
      call. = FALSE
    )
  }
  values
}

.ClusterFillScale <- function(cluster, ...) {
  ggplot2::scale_fill_manual(
    values = .ClusterScaleValues(cluster, ".ClusterFillScale"), ...)
}

.ClusterColourScale <- function(cluster, ...) {
  ggplot2::scale_colour_manual(
    values = .ClusterScaleValues(cluster, ".ClusterColourScale"), ...)
}

#' Plot a per-cluster diagnostic value
#'
#' @description Boxplot of any per-participant diagnostic (posterior
#' probability, distance to centroid, outlier score) split by cluster.
#' @param individual Per-participant diagnostic table containing `Cluster`.
#' @param value Diagnostic column name.
#' @param title Optional plot title.
#' @param noise_label Cluster value treated as noise, or `NULL` to disable.
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_MClust(
#'   df_Training, paste0("Var", 1:12), method = "finalize",
#'   final_k = 4, final_model = 1
#' )
#' PlotClusterDiagnostic(model$ProbFit$individual, "PosteriorMax")
#' }
#' @export
PlotClusterDiagnostic <- function(individual, value, title = NULL,
    noise_label = 0L) {
  df_Plot <- dplyr::filter(
    individual, !is.na(.data$Cluster), is.finite(.data[[value]]))
  df_Plot$.ClusterDisplay <- .ClusterFactor(df_Plot$Cluster, noise_label)
  value_label <- .ClusterMetricLabel(value)
  ggplot2::ggplot(
    df_Plot,
    ggplot2::aes(x = .data$.ClusterDisplay, y = .data[[value]])) +
    ggplot2::geom_boxplot() + ggplot2::theme_bw() +
    ggplot2::labs(
      title = if (is.null(title)) paste(value_label, "by cluster") else title,
      x = "Cluster", y = value_label)
}

# A compact distance diagnostic. The dashed line is the frozen upper training
# quantile used to triage projected cases; it is not a cluster boundary.
.ClusterDistancePlots <- function(individual, value, value_label,
    cutoff = NULL, prefix = "Training", noise_label = 0L) {
  df_Plot <- dplyr::filter(individual, is.finite(.data[[value]]))
  if (!nrow(df_Plot)) return(list())
  df_Plot$.ClusterDisplay <- .ClusterFactor(df_Plot$Cluster, noise_label)
  cutoff_layer <- function(plot, vertical) {
    if (is.null(cutoff) || !is.finite(cutoff)) return(plot)
    plot + (if (vertical) {
      ggplot2::geom_vline(
        xintercept = cutoff, linetype = "dashed", color = "#C92A2A")
    } else {
      ggplot2::geom_hline(
        yintercept = cutoff, linetype = "dashed", color = "#C92A2A")
    })
  }
  p_hist <- cutoff_layer(
    ggplot2::ggplot(df_Plot, ggplot2::aes(x = .data[[value]])) +
      ggplot2::geom_histogram(bins = 30, fill = "#3B5BDB") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste(prefix, tolower(value_label)),
        subtitle = if (is.null(cutoff) || !is.finite(cutoff)) NULL else
          "Dashed line marks the frozen training high-distance cutoff",
        x = value_label, y = "Count"),
    vertical = TRUE)
  list(distance_hist = p_hist)
}

.ClusterConfidencePlots <- function(individual, value, value_label,
    threshold = NULL, prefix = "Training", noise_label = 0L) {
  df_Plot <- dplyr::filter(
    individual, !is.na(.data$Cluster), is.finite(.data[[value]]))
  if (!nrow(df_Plot)) return(list())
  df_Plot$.ClusterDisplay <- .ClusterFactor(df_Plot$Cluster, noise_label)
  p_density <- ggplot2::ggplot(
    df_Plot, ggplot2::aes(x = .data[[value]])) +
    ggplot2::geom_density(fill = "darkblue", alpha = 0.7) +
    ggplot2::facet_wrap(~.ClusterDisplay, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Density of", tolower(value_label), "by cluster"),
      x = value_label, y = "Density")
  if (!is.null(threshold) && is.finite(threshold)) {
    p_density <- p_density + ggplot2::geom_vline(
      xintercept = threshold, linetype = "dashed", color = "#C92A2A")
  }
  list(confidence_density = p_density)
}

#' Plot a two-dimensional cluster review map
#'
#' @description Cluster assignments in a two-dimensional review space. The
#' space is frozen at training time and reused for projection so training and
#' projected cases are directly comparable. It is a display space only and does
#' not affect the clustering itself.
#' @param data Data frame containing coordinates and a cluster column.
#' @param x,y Coordinate variable names.
#' @param ClusterVar Cluster column name.
#' @param centroids Optional data frame of centroid coordinates in the same two
#'   columns, overlaid as crosses.
#' @param title,subtitle,xlab,ylab Plot annotations.
#' @param noise_label Cluster value treated as noise, or `NULL` to disable.
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
#' )
#' model$ModelInfo$plots$separation_map
#' PlotClusterMap(model$DataWithClusters, "Var1", "Var2", "Cluster")
#' }
#' @export
PlotClusterMap <- function(data, x, y, ClusterVar = "Cluster",
    centroids = NULL, title = "Cluster review map", subtitle = NULL,
    xlab = NULL, ylab = NULL, noise_label = 0L) {
  df_Plot <- data
  df_Plot$.ClusterDisplay <- .ClusterFactor(df_Plot[[ClusterVar]], noise_label)
  counts <- table(df_Plot$.ClusterDisplay)
  levels(df_Plot$.ClusterDisplay) <- paste0(
    names(counts), " (n = ", as.integer(counts), ", ",
    sprintf("%.0f%%", 100 * as.integer(counts) / sum(counts)), ")")
  plot <- ggplot2::ggplot(
    df_Plot,
    ggplot2::aes(
      x = .data[[x]], y = .data[[y]], color = .data$.ClusterDisplay)) +
    ggplot2::geom_point(alpha = 0.75) +
    .ClusterColourScale(df_Plot$.ClusterDisplay) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title, subtitle = subtitle, color = "Cluster",
      x = if (is.null(xlab)) .ClusterVariableLabel(data, x) else xlab,
      y = if (is.null(ylab)) .ClusterVariableLabel(data, y) else ylab)
  if (!is.null(centroids) && nrow(centroids)) {
    plot <- plot + ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(x = .data[[names(centroids)[[1]]]],
                   y = .data[[names(centroids)[[2]]]]),
      inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.2, color = "black")
  }
  plot
}

#' @rdname PlotClusterMap
#' @export
PlotClusterAssignment <- function(data, x, y, ClusterVar = "Cluster") {
  PlotClusterMap(data, x, y, ClusterVar, title = "Cluster assignment")
}

# A frozen PCA gives every numeric method the same kind of two-dimensional
# review space the SOM grid gives the SOM pipeline, and new cases can be placed
# in it without refitting.
.ClusterReviewSpace <- function(X) {
  if (is.null(X) || ncol(X) < 2L || nrow(X) < 3L) return(NULL)
  fit <- tryCatch(
    stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = 2L),
    error = function(e) NULL)
  if (is.null(fit) || ncol(fit$rotation) < 2L) return(NULL)
  variance <- fit$sdev^2
  list(prcomp = fit,
       variance_explained = variance[1:2] / sum(variance))
}

.ClusterReviewCoordinates <- function(review_space, X) {
  if (is.null(review_space) || is.null(X) || !nrow(X)) return(NULL)
  coordinates <- as.data.frame(
    stats::predict(review_space$prcomp, newdata = X)[, 1:2, drop = FALSE])
  names(coordinates) <- c("ReviewDimension1", "ReviewDimension2")
  coordinates
}

.ClusterReviewAxisLabels <- function(review_space) {
  if (is.null(review_space$variance_explained)) {
    return(c("Review dimension 1", "Review dimension 2"))
  }
  paste0("Review dimension ", 1:2, " (",
         sprintf("%.0f%%", 100 * review_space$variance_explained), " variance)")
}

.ClusterSeparationMap <- function(review_space, X, cluster, centers = NULL,
    prefix = "Training", noise_label = 0L) {
  coordinates <- .ClusterReviewCoordinates(review_space, X)
  if (is.null(coordinates)) return(NULL)
  coordinates$Cluster <- cluster
  centroid_coordinates <- if (is.null(centers)) NULL else
    .ClusterReviewCoordinates(review_space, centers)
  axis_labels <- .ClusterReviewAxisLabels(review_space)
  PlotClusterMap(
    coordinates, "ReviewDimension1", "ReviewDimension2", "Cluster",
    centroids = centroid_coordinates,
    title = paste(prefix, "cluster review map"),
    subtitle = paste0(
      "Training-frozen principal-component review space",
      if (is.null(centroid_coordinates)) "" else "; black crosses mark centres"),
    xlab = axis_labels[[1]], ylab = axis_labels[[2]],
    noise_label = noise_label)
}

#' Plot a per-participant silhouette profile
#'
#' @description The classic silhouette profile: one bar per participant, sorted
#' within cluster, with the average silhouette width marked. Bars near zero or
#' below sit closer to a neighbouring cluster than their own.
#' @param silhouette A `cluster::silhouette()` object or a matrix with
#'   `cluster`, `neighbor`, and `sil_width` columns.
#' @param title Plot title.
#' @return A `ggplot` object, or `NULL` when silhouette widths are unavailable.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
#' )
#' model$ModelInfo$plots$silhouette
#' }
#' @export
PlotClusterSilhouette <- function(silhouette, title = "Silhouette profile") {
  if (is.null(silhouette)) return(NULL)
  df_Silhouette <- as.data.frame(unclass(silhouette))
  if (!all(c("cluster", "sil_width") %in% names(df_Silhouette))) return(NULL)
  df_Silhouette <- df_Silhouette[is.finite(df_Silhouette$sil_width), , drop = FALSE]
  if (!nrow(df_Silhouette)) return(NULL)
  df_Silhouette$Cluster <- .ClusterFactor(df_Silhouette$cluster, NULL)
  df_Silhouette <- df_Silhouette[
    order(df_Silhouette$Cluster, -df_Silhouette$sil_width), , drop = FALSE]
  df_Silhouette$Order <- seq_len(nrow(df_Silhouette))
  average_width <- mean(df_Silhouette$sil_width)
  cluster_means <- df_Silhouette %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      MeanWidth = mean(.data$sil_width), n = dplyr::n(), .groups = "drop")
  ggplot2::ggplot(
    df_Silhouette,
    ggplot2::aes(x = .data$Order, y = .data$sil_width, fill = .data$Cluster)) +
    ggplot2::geom_col(width = 1) +
    .ClusterFillScale(df_Silhouette$Cluster) +
    ggplot2::geom_hline(
      yintercept = average_width, linetype = "dashed", color = "#C92A2A") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_reverse() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Average width ", sprintf("%.2f", average_width),
        " (dashed); per cluster ",
        paste(sprintf("%s: %.2f", cluster_means$Cluster, cluster_means$MeanWidth),
              collapse = ", ")),
      x = "Participant", y = "Silhouette width", fill = "Cluster")
}

.ClusterCenterTable <- function(centers, variable_labels = NULL) {
  if (is.null(centers) || !nrow(centers)) return(NULL)
  df_Centers <- as.data.frame(centers)
  variable_names <- names(df_Centers)
  if (is.null(variable_labels) || length(variable_labels) != length(variable_names)) {
    variable_labels <- variable_names
  }
  df_Centers$Cluster <- .ClusterFactor(
    if (is.null(rownames(centers))) seq_len(nrow(centers)) else rownames(centers),
    NULL)
  df_Long <- tidyr::pivot_longer(
    df_Centers, -dplyr::all_of("Cluster"),
    names_to = "Variable", values_to = "Value")
  df_Long$Variable <- factor(
    df_Long$Variable, levels = variable_names, labels = variable_labels)
  df_Long
}

#' Plot cluster centre profiles as a heatmap
#'
#' @description Cluster centres across clustering variables. Values are the
#' centres in the frozen analysis scale, so a centred and scaled model reads
#' directly as standard deviations from the cohort mean.
#' @param centers Matrix or data frame of cluster centres, one row per cluster.
#' @param variable_labels Optional display labels for the columns.
#' @param title Plot title.
#' @param value_label Legend title describing the centre scale.
#' @return A `ggplot` object, or `NULL` when centres are unavailable.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
#' )
#' model$ModelInfo$plots$centre_heatmap
#' }
#' @export
PlotClusterCentreHeatmap <- function(centers, variable_labels = NULL,
    title = "Cluster centre profiles", value_label = "Centre") {
  df_Long <- .ClusterCenterTable(centers, variable_labels)
  if (is.null(df_Long)) return(NULL)
  limit <- max(abs(df_Long$Value), na.rm = TRUE)
  ggplot2::ggplot(
    df_Long,
    ggplot2::aes(x = .data$Variable, y = .data$Cluster, fill = .data$Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#1864AB", mid = "white", high = "#C92A2A", midpoint = 0,
      limits = c(-limit, limit)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = title, x = "Variable", y = "Cluster", fill = value_label)
}

#' @description `PlotClusterCentreProfile()` shows the same centres as connected
#'   lines, which reads more like the SOM line map when variables have a
#'   meaningful order.
#' @rdname PlotClusterCentreHeatmap
#' @export
PlotClusterCentreProfile <- function(centers, variable_labels = NULL,
    title = "Cluster centre profiles", value_label = "Centre") {
  df_Long <- .ClusterCenterTable(centers, variable_labels)
  if (is.null(df_Long)) return(NULL)
  ggplot2::ggplot(
    df_Long,
    ggplot2::aes(
      x = .data$Variable, y = .data$Value, color = .data$Cluster,
      group = .data$Cluster)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = title, x = "Variable", y = value_label, color = "Cluster")
}

#' Plot labelled numeric profiles by cluster
#' @inheritParams PlotClusterBoxplot
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), method = "finalize", final_k = 4
#' )
#' PlotClusterProfiles(model$DataWithClusters, "Cluster", paste0("Var", 1:12))
#' }
#' @export
PlotClusterProfiles <- function(data, ClusterVar, variables, codebook = NULL) {
  PlotClusterBoxplot(
    data = dplyr::filter(data, !is.na(.data[[ClusterVar]])),
    ClusterVar = ClusterVar, variables = variables, codebook = codebook,
    FillTitle = "Measure") +
    ggplot2::labs(title = "Numeric cluster profiles")
}

#' Plot categorical composition by cluster
#' @param data Data frame containing the categorical variables.
#' @param variables Categorical variable names.
#' @param cluster Cluster assignment vector aligned to `data`.
#' @param facet_by Whether stacked-bar facets represent categorical variables
#'   (the default) or clusters.
#' @param style Either `"stacked"` for composition bars or `"enrichment"` for
#'   a cluster-first heatmap of percentage-point differences from the cohort.
#' @return A `ggplot` object, or `NULL` when no categorical variables are given.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_Gower_PAM(
#'   df_Training, c("Var1", "Var2", "CatVar1", "CatVar2"),
#'   method = "finalize", final_k = 4
#' )
#' model$ModelInfo$plots$categorical_composition
#' }
#' @export
PlotClusterComposition <- function(data, variables, cluster,
    facet_by = c("variable", "cluster"), style = c("stacked", "enrichment")) {
  if (!length(variables)) return(NULL)
  facet_by <- match.arg(facet_by)
  style <- match.arg(style)
  df_Categorical <- dplyr::bind_rows(lapply(variables, function(variable) {
    dplyr::tibble(
      Cluster = .ClusterFactor(cluster, NULL),
      Variable = .ClusterVariableLabel(data, variable),
      Value = sub("^.*: ", "", as.character(data[[variable]])))
  })) %>%
    dplyr::filter(!is.na(.data$Cluster), !is.na(.data$Value)) %>%
    dplyr::count(.data$Cluster, .data$Variable, .data$Value) %>%
    dplyr::group_by(.data$Cluster, .data$Variable) %>%
    dplyr::mutate(Proportion = .data$n / sum(.data$n)) %>%
    dplyr::ungroup()
  if (!nrow(df_Categorical)) return(NULL)
  if (identical(style, "enrichment")) {
    overall <- df_Categorical %>%
      dplyr::group_by(.data$Variable, .data$Value) %>%
      dplyr::summarise(n = sum(.data$n), .groups = "drop") %>%
      dplyr::group_by(.data$Variable) %>%
      dplyr::mutate(OverallProportion = .data$n / sum(.data$n)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"n")
    df_Categorical <- df_Categorical %>%
      dplyr::left_join(overall, by = c("Variable", "Value")) %>%
      dplyr::mutate(
        Enrichment = .data$Proportion - .data$OverallProportion,
        Category = paste0(.data$Variable, ": ", .data$Value)) %>%
      dplyr::group_by(.data$Cluster) %>%
      dplyr::arrange(dplyr::desc(abs(.data$Enrichment)), .by_group = TRUE) %>%
      dplyr::mutate(Category = factor(.data$Category, levels = rev(unique(.data$Category)))) %>%
      dplyr::ungroup()
    limit <- max(abs(df_Categorical$Enrichment), na.rm = TRUE)
    return(ggplot2::ggplot(
      df_Categorical,
      ggplot2::aes(x = 1, y = .data$Category, fill = .data$Enrichment)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = scales::percent(.data$Proportion, accuracy = 1)),
        size = 3) +
      ggplot2::facet_wrap(~Cluster, scales = "free_y") +
      ggplot2::scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030",
        midpoint = 0, limits = c(-limit, limit), labels = scales::label_percent()) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(title = "Categorical features enriched within each cluster",
        subtitle = "Tile colour = percentage-point difference from the full cohort; text = within-cluster percentage",
        x = NULL, y = NULL, fill = "Difference from cohort"))
  }
  facet_formula <- if (facet_by == "variable") ~Variable else ~Cluster
  x_var <- if (facet_by == "variable") "Cluster" else "Variable"

  ggplot2::ggplot(
    df_Categorical,
    ggplot2::aes(x = .data[[x_var]], y = .data$Proportion, fill = .data$Value)) +
    ggplot2::geom_col() +
    .SciDataFillScale() +
    ggplot2::facet_wrap(facet_formula) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Categorical composition by cluster",
      x = if (facet_by == "variable") "Cluster" else "Variable",
      y = "Participants", fill = "Category")
}

#' Plot cluster fit-review metrics
#' @param fit_table Candidate-model fit table.
#' @param x Candidate-count column name.
#' @param metrics Numeric metric columns to display. When `NULL`, a concise
#'   set of decision-relevant raw metrics is selected automatically.
#' @param group Optional candidate-family column, such as `"Model"` for
#'   Mclust or `"Epsilon"` for HDBSCAN. When `NULL`, it is inferred when
#'   possible.
#' @param title Plot title.
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(SimulatedPhenotypeData)
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' model <- CreateClusterModel_KMeans(
#'   df_Training, paste0("Var", 1:12), k_range = 2:4
#' )
#' PlotClusterFitReview(model$ModelInfo$fit_table)
#' }
#' @export
PlotClusterFitReview <- function(fit_table, x = "Classes", metrics = NULL,
    group = NULL, title = "Candidate model review") {
  if (!x %in% names(fit_table)) {
    stop("x was not found in fit_table: ", x)
  }
  if (is.null(metrics)) {
    review_metrics <- c(
      "BIC", "ICL", "AIC", "Entropy", "WSS", "Silhouette",
      "CalinskiHarabasz", "Persistence", "MeanMembershipProbability",
      "NoiseProportion", "AverageSilhouette", "MeanMedoidDistance",
      "MinClusterN", "MinClassN", "SizeBalance", "ReproducibilityScore",
      "MaxUncertainty")
    metrics <- intersect(review_metrics, names(fit_table))
  }
  metrics <- metrics[metrics %in% names(fit_table)]
  metrics <- metrics[vapply(fit_table[metrics], is.numeric, logical(1))]
  metrics <- metrics[vapply(
    fit_table[metrics], function(x) any(is.finite(x)), logical(1))]
  metrics <- metrics[!grepl("_scaled$", metrics)]
  metrics <- setdiff(metrics, c("StabilitySuccessRate", "Classes", "MinPts", "Epsilon"))
  metrics <- setdiff(metrics, x)
  if (!length(metrics)) stop("No numeric fit metrics are available to plot.")
  if (is.null(group)) {
    group <- if ("Model" %in% names(fit_table)) {
      "Model"
    } else if ("Epsilon" %in% names(fit_table) && x != "Epsilon") {
      "Epsilon"
    } else NULL
  }
  if (!is.null(group) && !group %in% names(fit_table)) {
    stop("group was not found in fit_table: ", group)
  }
  df_Long <- tidyr::pivot_longer(
    fit_table, dplyr::all_of(metrics), names_to = "Metric", values_to = "Value")
  if (is.null(group)) {
    plot <- ggplot2::ggplot(
      df_Long, ggplot2::aes(x = .data[[x]], y = .data$Value)) +
      ggplot2::geom_point()
    if (nrow(fit_table) > 1) {
      plot <- plot + ggplot2::geom_line(ggplot2::aes(group = .data$Metric))
    }
  } else {
    plot <- ggplot2::ggplot(
      df_Long,
      ggplot2::aes(
        x = .data[[x]], y = .data$Value, color = factor(.data[[group]]),
        group = interaction(.data$Metric, .data[[group]]))) +
      ggplot2::geom_point()
    if (nrow(fit_table) > 1) plot <- plot + ggplot2::geom_line()
    plot <- plot + .SciDataColourScale() + ggplot2::labs(color = group)
  }
  plot <- plot +
    ggplot2::facet_wrap(
      ~Metric, scales = "free_y", ncol = 1,
      labeller = ggplot2::labeller(Metric = .ClusterMetricLabel)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title, x = .ClusterMetricLabel(x), y = "Fit metric")
  plot
}

# A single-metric raw candidate curve. Used for the elbow, silhouette, BIC,
# and persistence curves each method is judged on.
.ClusterCandidateCurve <- function(fit_table, metric, title, ylab = NULL,
    x = "Classes") {
  if (!metric %in% names(fit_table) ||
      !any(is.finite(fit_table[[metric]]))) return(NULL)
  PlotClusterFitReview(fit_table, x = x, metrics = metric, title = title) +
    ggplot2::facet_null() +
    ggplot2::labs(y = if (is.null(ylab)) .ClusterMetricLabel(metric) else ylab)
}

.ClusterStabilityPlots <- function(Stability) {
  if (is.null(Stability) || is.null(Stability$replicates) ||
      !nrow(Stability$replicates)) return(list())
  identifiers <- intersect(
    c("Model", "Classes", "MinPts", "Epsilon"), names(Stability$replicates))
  CandidateLabel <- function(df_Candidate) {
    if (!length(identifiers)) return(rep("Final model", nrow(df_Candidate)))
    apply(df_Candidate[identifiers], 1, function(row) {
      paste(paste0(identifiers, " = ", trimws(row)), collapse = ", ")
    })
  }
  plots <- list()
  if (!is.null(Stability$cluster_recovery) &&
      nrow(Stability$cluster_recovery)) {
    df_Recovery <- Stability$cluster_recovery
    df_Recovery$Candidate <- CandidateLabel(df_Recovery)
    df_Recovery$ClusterDisplay <- .ClusterFactor(df_Recovery$Cluster)
    plots$cluster_recovery <- ggplot2::ggplot(
      df_Recovery,
      ggplot2::aes(x = .data$ClusterDisplay, y = .data$Jaccard)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~Candidate) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Per-cluster bootstrap recovery",
        x = "Cluster", y = "Jaccard")
  }
  metric_names <- intersect(c("ARI", "VI", "NMI", "FowlkesMallows"),
    names(Stability$replicates))
  if (length(metric_names)) {
    df_Metrics <- Stability$replicates %>%
      dplyr::filter(.data$Status == "success") %>%
      tidyr::pivot_longer(dplyr::all_of(metric_names), names_to = "Metric",
        values_to = "Value")
    if (nrow(df_Metrics)) {
      plots$partition_metrics <- ggplot2::ggplot(
        df_Metrics, ggplot2::aes(x = .data$Value)) +
        ggplot2::geom_histogram(bins = 20, fill = "#3B5BDB") +
        ggplot2::facet_wrap(~.data$Metric, scales = "free") +
        ggplot2::theme_bw() + ggplot2::labs(
          title = "Bootstrap partition agreement metrics", x = NULL, y = "Refits")
    }
  }
  if (!is.null(Stability$cluster_inclusion) && nrow(Stability$cluster_inclusion)) {
    plots$cluster_inclusion <- ggplot2::ggplot(
      Stability$cluster_inclusion,
      ggplot2::aes(x = factor(.data$Cluster), y = .data$MeanInclusion)) +
      ggplot2::geom_col(fill = "#3B5BDB") +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = .data$P05Inclusion, ymax = .data$MeanInclusion), width = .2) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) + ggplot2::theme_bw() +
      ggplot2::labs(title = "Cluster inclusion across bootstrap refits",
        x = "Reference cluster", y = "Inclusion probability")
  }
  if (!is.null(Stability$coassignment) &&
      identical(Stability$coassignment$status, "available")) {
    matrix <- Stability$coassignment$matrix
    if (!is.null(matrix)) {
      df_Matrix <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
      names(df_Matrix) <- c("Participant1", "Participant2", "Probability")
      plots$coassignment <- ggplot2::ggplot(df_Matrix,
        ggplot2::aes(x = .data$Participant1, y = .data$Participant2,
          fill = .data$Probability)) +
        ggplot2::geom_raster() + ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
        ggplot2::theme_bw() + ggplot2::theme(
          axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(title = "Bootstrap co-assignment probability", x = "Participant", y = "Participant")
    }
  }
  plots
}

# Every projection is triaged the same way the SOM pipeline triages it: distance
# beyond the frozen training cutoff, low membership confidence, or both.
.ClusterProjectionFitClass <- function(distance, confidence, distance_cutoff,
    confidence_threshold, higher_confidence_is_better = TRUE) {
  poor_distance <- !is.na(distance) & is.finite(distance_cutoff) &
    distance > distance_cutoff
  uncertain <- if (is.null(confidence) || !is.finite(confidence_threshold)) {
    rep(FALSE, length(distance))
  } else if (higher_confidence_is_better) {
    !is.na(confidence) & confidence < confidence_threshold
  } else {
    !is.na(confidence) & confidence > confidence_threshold
  }
  fit_class <- dplyr::case_when(
    is.na(distance) ~ NA_character_,
    poor_distance & uncertain ~ "Potential novel phenotype",
    poor_distance ~ "Poor fit to training structure",
    uncertain ~ "Uncertain membership",
    TRUE ~ "Good fit")
  factor(fit_class, levels = c(
    "Good fit", "Uncertain membership", "Poor fit to training structure",
    "Potential novel phenotype"))
}

.ClusterFitClassPlot <- function(fit_class) {
  df_Class <- dplyr::tibble(FitClass = fit_class)
  df_Class <- dplyr::filter(df_Class, !is.na(.data$FitClass))
  if (!nrow(df_Class)) return(NULL)
  df_Class <- df_Class %>%
    dplyr::count(.data$FitClass, name = "n", .drop = FALSE) %>%
    dplyr::mutate(
      Proportion = .data$n / sum(.data$n),
      Label = paste0(.data$n, " (", sprintf("%.0f%%", 100 * .data$Proportion), ")"))
  ggplot2::ggplot(
    df_Class,
    ggplot2::aes(x = .data$FitClass, y = .data$n, fill = .data$FitClass)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), vjust = -0.35) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::scale_fill_manual(values = c(
      "Good fit" = "#2F9E44",
      "Uncertain membership" = "#F59F00",
      "Poor fit to training structure" = "#E8590C",
      "Potential novel phenotype" = "#C92A2A"), drop = FALSE) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::labs(
      title = "Projection fit against the frozen training reference",
      x = "Projection fit class", y = "Participants")
}
