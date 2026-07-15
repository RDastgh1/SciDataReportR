#' Make a pairwise referent heatmap
#'
#' Build a heatmap of pairwise group contrasts against a required referent
#' group. Continuous outcomes are always transformed using the referent group
#' before modeling: Z-scores when `Parametric = TRUE`, and M-scores when
#' `Parametric = FALSE`.
#'
#' @param data A data frame.
#' @param group_var Character scalar naming the grouping variable.
#' @param variables Character vector of continuous outcome variables.
#' @param Referent Character scalar naming the referent level of `group_var`.
#' @param covariates Optional character vector of covariates.
#' @param Parametric Logical. If `TRUE`, outcomes are Z-scored before modeling.
#'   If `FALSE`, outcomes are M-scored and HC3 robust covariance is used for
#'   estimated marginal mean contrasts.
#' @param show_caption Logical; add an explanatory caption to the plot.
#' @param x_axis_text_angle Numeric angle for x-axis labels. Defaults to `0`.
#' @param adjust_scope Multiple-comparison correction scope. `"per_group"`
#'   adjusts across variables within each group-vs-referent contrast;
#'   `"per_variable"` adjusts across group contrasts within each variable;
#'   `"matrix"` adjusts across all displayed cells; `"none"` applies no
#'   correction.
#' @param p_adjust_method Method passed to [stats::p.adjust()]. Use `"none"`
#'   for no correction.
#' @param star_p Which p-values should drive cell stars: raw, adjusted, or none.
#' @param adjusted_outline Logical; outline cells significant after adjustment.
#' @param adjusted_significance_threshold Threshold for adjusted-significant
#'   outlines.
#' @param adjusted_outline_color,adjusted_outline_linewidth Appearance of the
#'   adjusted-significant outline.
#' @param low_color,mid_color,high_color Diverging heatmap colors.
#' @param fill_midpoint Numeric midpoint for the fill scale.
#' @param fill_limits Optional numeric vector of length 2. If `NULL`, symmetric
#'   limits are computed from the observed estimated mean differences.
#' @param fill_oob Out-of-bounds handler for the fill scale.
#' @param cluster_rows,cluster_columns Logical; optionally cluster rows or
#'   columns based on estimated mean differences.
#' @param return_models Logical; include fitted model objects in the return.
#' @param star_color,star_size Appearance of p-value stars.
#'
#' @return An object of class `"SciDataReportRPairwiseHeatmap"` with `Plot`,
#'   `Results`, `Models`, `Settings`, `ScalingParameters`, and `Warnings`.
#'   `Results` includes readable audit columns such as `Test`, `Contrast`,
#'   `Adjustment`, and `ModelFormula`.
#'
#' @details
#' Each cell is an estimated marginal mean contrast from the model
#' `referent_scaled_outcome ~ group_var + covariates`, computed as
#' `Group - Referent`. Scaling parameters are estimated in the referent group
#' and projected onto the full dataset before modeling. With the defaults,
#' adjusted p-values use FDR correction within each group-vs-referent contrast
#' across variables.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' res <- MakePairwiseHeatmap(
#'   data = Labelled,
#'   group_var = "Diagnosis",
#'   variables = c("AXL", "Adiponectin", "Cortisol"),
#'   Referent = levels(Labelled$Diagnosis)[1]
#' )
#' res$Plot
#'
#' @export
MakePairwiseHeatmap <- function(data,
                                group_var,
                                variables,
                                Referent,
                                covariates = NULL,
                                Parametric = TRUE,
                                adjust_scope = c("per_group", "per_variable", "matrix", "none"),
                                p_adjust_method = c("fdr", "bonferroni", "holm", "none"),
                                star_p = c("raw", "adjusted", "none"),
                                adjusted_outline = TRUE,
                                adjusted_significance_threshold = 0.05,
                                adjusted_outline_color = "black",
                                adjusted_outline_linewidth = 1.0,
                                low_color = "#52BCA3FF",
                                mid_color = "white",
                                high_color = "#E58606FF",
                                fill_midpoint = 0,
                                fill_limits = NULL,
                                fill_oob = scales::squish,
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                show_caption = FALSE,
                                x_axis_text_angle = 0,
                                return_models = FALSE,
                                star_color = "black",
                                star_size = 4) {
  adjust_scope <- match.arg(adjust_scope)
  p_adjust_method <- match.arg(p_adjust_method)
  star_p <- match.arg(star_p)

  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.character(group_var) || length(group_var) != 1) {
    stop("group_var must be a single character string.")
  }
  if (!group_var %in% names(data)) {
    stop("group_var was not found in data: ", group_var)
  }
  if (!is.character(variables) || length(variables) == 0) {
    stop("variables must be a non-empty character vector.")
  }
  if (!all(variables %in% names(data))) {
    stop(
      "Variable(s) not found: ",
      paste(setdiff(variables, names(data)), collapse = ", ")
    )
  }
  if (missing(Referent) || is.null(Referent) || !is.character(Referent) || length(Referent) != 1) {
    stop("Referent must be supplied as a single character level name.")
  }
  if (!is.logical(Parametric) || length(Parametric) != 1 || is.na(Parametric)) {
    stop("Parametric must be TRUE or FALSE.")
  }
  if (!is.logical(cluster_rows) || length(cluster_rows) != 1 || is.na(cluster_rows)) {
    stop("cluster_rows must be TRUE or FALSE.")
  }
  if (!is.logical(cluster_columns) || length(cluster_columns) != 1 || is.na(cluster_columns)) {
    stop("cluster_columns must be TRUE or FALSE.")
  }
  if (!is.logical(show_caption) || length(show_caption) != 1 || is.na(show_caption)) {
    stop("show_caption must be TRUE or FALSE.")
  }
  if (!is.numeric(x_axis_text_angle) || length(x_axis_text_angle) != 1 || !is.finite(x_axis_text_angle)) {
    stop("x_axis_text_angle must be a single finite numeric value.")
  }
  if (!is.null(covariates) && !all(covariates %in% names(data))) {
    stop(
      "Covariate(s) not found: ",
      paste(setdiff(covariates, names(data)), collapse = ", ")
    )
  }
  if (!is.null(fill_limits)) {
    if (!is.numeric(fill_limits) || length(fill_limits) != 2 || any(!is.finite(fill_limits))) {
      stop("fill_limits must be NULL or a finite numeric vector of length 2.")
    }
    if (fill_limits[[1]] >= fill_limits[[2]]) {
      stop("fill_limits must be ordered from low to high.")
    }
  }

  variables <- unique(variables)
  score_type <- if (isTRUE(Parametric)) "ZScore" else "MScore"
  score_prefix <- if (isTRUE(Parametric)) ".PairwiseHeatmap_Z_" else ".PairwiseHeatmap_M_"

  numeric_ok <- vapply(data[variables], is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop(
      "MakePairwiseHeatmap currently requires continuous numeric variables. Non-numeric: ",
      paste(variables[!numeric_ok], collapse = ", ")
    )
  }

  group_for_scaling <- if (is.factor(data[[group_var]])) {
    droplevels(data[[group_var]])
  } else if (is.logical(data[[group_var]])) {
    factor(data[[group_var]], levels = c(FALSE, TRUE))
  } else {
    factor(data[[group_var]])
  }

  if (!Referent %in% levels(group_for_scaling)) {
    stop("Referent level not found: ", Referent)
  }

  referent_data <- data[group_for_scaling == Referent & !is.na(group_for_scaling), , drop = FALSE]
  if (nrow(referent_data) < 2) {
    stop("At least two referent rows are required to estimate scaling parameters.")
  }

  scale_obj <- if (isTRUE(Parametric)) {
    referent_scale_obj <- CreateZScoreObject(
      data = referent_data,
      variables = variables,
      names_prefix = score_prefix,
      RetainLabels = FALSE,
      RenameLabels = FALSE
    )
    ProjectZScore(
      data = data,
      variables = variables,
      parameters = referent_scale_obj,
      ParameterInputType = "ZScoreObj",
      names_prefix = score_prefix,
      RetainLabels = FALSE,
      RenameLabels = FALSE
    )
  } else {
    referent_scale_obj <- CreateMScoreObject(
      data = referent_data,
      variables = variables,
      names_prefix = score_prefix,
      RetainLabels = FALSE,
      RenameLabels = FALSE
    )
    .ProjectMScore(
      data = data,
      variables = variables,
      parameters = referent_scale_obj$Parameters,
      names_prefix = score_prefix,
      center = referent_scale_obj$Center,
      scale = referent_scale_obj$Scale
    )
  }

  scaled_data <- if (isTRUE(Parametric)) {
    scale_obj$DataWithZ
  } else {
    scale_obj$DataWithM
  }

  transformed_variables <- paste0(score_prefix, variables)

  contrast_res <- .ComputePairwiseReferentContrasts(
    data = scaled_data,
    group_var = group_var,
    variables = variables,
    Referent = Referent,
    covariates = covariates,
    transformed_variables = transformed_variables,
    score_type = score_type,
    Parametric = Parametric,
    adjust_scope = adjust_scope,
    p_adjust_method = p_adjust_method,
    adjusted_significance_threshold = adjusted_significance_threshold,
    star_p = star_p,
    return_models = return_models
  )

  results <- contrast_res$Results

  if (!nrow(results)) {
    empty_plot <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        subtitle = "No pairwise referent contrasts were available."
      )

    out <- list(
      Plot = empty_plot,
      Results = results,
      Models = contrast_res$Models,
      Settings = list(),
      ScalingParameters = scale_obj$Parameters,
      Warnings = contrast_res$Warnings
    )
    class(out) <- c("SciDataReportRPairwiseHeatmap", class(out))
    return(out)
  }

  group_levels <- if (is.factor(scaled_data[[group_var]])) {
    levels(droplevels(scaled_data[[group_var]]))
  } else {
    levels(factor(scaled_data[[group_var]]))
  }
  group_order <- setdiff(group_levels, Referent)
  row_order <- variables

  if (isTRUE(cluster_rows)) {
    wide_rows <- stats::xtabs(
      EstimatedMeanDifference ~ Variable + Group,
      data = results
    )
    if (nrow(wide_rows) > 1 && ncol(wide_rows) > 0 && all(is.finite(wide_rows))) {
      row_order <- rownames(wide_rows)[stats::hclust(stats::dist(wide_rows))$order]
    }
  }

  if (isTRUE(cluster_columns)) {
    wide_cols <- stats::xtabs(
      EstimatedMeanDifference ~ Variable + Group,
      data = results
    )
    if (ncol(wide_cols) > 1 && nrow(wide_cols) > 0 && all(is.finite(wide_cols))) {
      group_order <- colnames(wide_cols)[stats::hclust(stats::dist(t(wide_cols)))$order]
    }
  }

  row_labels <- stats::setNames(
    vapply(row_order, function(v) .SdrLabelOrName(data, v), character(1)),
    row_order
  )

  results$RowLabel <- row_labels[results$Variable]
  results$ColumnLabel <- results$GroupLabel
  results$RowLabel <- factor(results$RowLabel, levels = rev(row_labels))
  results$ColumnLabel <- factor(results$ColumnLabel, levels = group_order)

  fill_limits_used <- if (is.null(fill_limits)) {
    .ComputeSymmetricFillLimits(results$EstimatedMeanDifference)
  } else {
    fill_limits
  }

  caption <- if (isTRUE(show_caption)) {
    .BuildHeatmapCaption(
      score_type = score_type,
      star_p = star_p,
      adjusted_outline = adjusted_outline
    )
  } else {
    NULL
  }

  plot <- ggplot2::ggplot(
    results,
    ggplot2::aes(
      x = .data$ColumnLabel,
      y = .data$RowLabel,
      fill = .data$EstimatedMeanDifference
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    .GetHeatmapColorScale(
      low_color = low_color,
      mid_color = mid_color,
      high_color = high_color,
      fill_midpoint = fill_midpoint,
      fill_limits = fill_limits_used,
      fill_oob = fill_oob,
      name = paste0("Referent-scaled\n", ifelse(score_type == "ZScore", "Z-score", "M-score"))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = x_axis_text_angle,
        hjust = if (x_axis_text_angle == 0) 0.5 else 1,
        vjust = if (x_axis_text_angle == 0) 0.5 else 1
      ),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(caption = caption)

  plot <- .AddHeatmapSignificanceLayers(
    plot = plot,
    data = results,
    x_col = "ColumnLabel",
    y_col = "RowLabel",
    label_col = "SignificanceLabel",
    outline_col = "IsAdjustedSignificant",
    adjusted_outline = adjusted_outline,
    adjusted_outline_color = adjusted_outline_color,
    adjusted_outline_linewidth = adjusted_outline_linewidth,
    star_color = star_color,
    star_size = star_size
  )

  settings <- list(
    group_var = group_var,
    variables = variables,
    Referent = Referent,
    covariates = covariates,
    Parametric = Parametric,
    ScoreType = score_type,
    ScalingReference = "Referent",
    adjust_scope = adjust_scope,
    p_adjust_method = p_adjust_method,
    star_p = star_p,
    adjusted_outline = adjusted_outline,
    adjusted_significance_threshold = adjusted_significance_threshold,
    fill_limits = fill_limits,
    fill_limits_used = fill_limits_used,
    fill_midpoint = fill_midpoint,
    colors = c(
      low_color = low_color,
      mid_color = mid_color,
      high_color = high_color
    ),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_caption = show_caption,
    x_axis_text_angle = x_axis_text_angle
  )

  out <- list(
    Plot = plot,
    Results = results,
    Models = contrast_res$Models,
    Settings = settings,
    ScalingParameters = scale_obj$Parameters,
    Warnings = contrast_res$Warnings
  )

  class(out) <- c("SciDataReportRPairwiseHeatmap", class(out))
  out
}

#' @export
print.SciDataReportRPairwiseHeatmap <- function(x, ...) {
  print(x$Plot)
  invisible(x)
}
