#' Plot a diagnostic likelihood-ratio heatmap
#'
#' Visualizes results created by [DiagnosticLikelihoodRatioTable()] without
#' recomputing diagnostic statistics. Tile fill is log2(LR), making reciprocal
#' evidence (for example, LR 0.25 and LR 4) equally distant from LR 1.
#'
#' @param x An object returned by [DiagnosticLikelihoodRatioTable()], or a tidy
#'   data frame compatible with its `Results` element.
#' @param result Diagnostic result levels to display: `"all"`, `"positive"`, or
#'   `"negative"`. Positive and negative selections require the result object.
#' @param predictor_order Predictor row ordering: `"original"`, `"alphabetical"`,
#'   `"strength"`, or `"cluster"`.
#' @param outcome_order Outcome column ordering with the same choices.
#' @param orientation Tile orientation: `"auto"`, `"predictors_rows"`, or
#'   `"outcomes_rows"`.
#' @param show_values `"auto"`, `TRUE`, or `FALSE`; auto labels at most 150 tiles.
#' @param show_ci_marker Logical; append `*` when the unadjusted LR CI excludes 1.
#' @param facet_strata Logical; facet separate diagnostic strata when present.
#' @param cap Optional positive maximum absolute log2(LR) for color scaling.
#' @param na_color Fill color for unavailable LR values.
#'
#' @return A static `ggplot`. Its `DiagnosticLRData` attribute and tile `text`
#'   aesthetic contain complete hover-ready information for optional Plotly use.
#'
#' @seealso [DiagnosticLikelihoodRatioTable()] to calculate the displayed LRs.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' df_Labelled$DiagnosisBinary <- factor(df_Labelled$Diagnosis,
#'   levels = c("Control", "Impaired"))
#' lr <- DiagnosticLikelihoodRatioTable(df_Labelled, "DiagnosisBinary",
#'   c("sex", "Genotype"))
#' PlotDiagnosticLRHeatmap(lr, result = "all")
#' @export
PlotDiagnosticLRHeatmap <- function(x, result = c("all", "positive", "negative"),
    predictor_order = c("original", "alphabetical", "strength", "cluster"),
    outcome_order = c("original", "alphabetical", "strength", "cluster"),
    orientation = c("auto", "predictors_rows", "outcomes_rows"),
    show_values = "auto", show_ci_marker = TRUE, facet_strata = TRUE,
    cap = NULL, na_color = "grey90") {
  result <- match.arg(result)
  predictor_order <- match.arg(predictor_order)
  outcome_order <- match.arg(outcome_order)
  orientation <- match.arg(orientation)
  if (!is.logical(show_ci_marker) || length(show_ci_marker) != 1 || is.na(show_ci_marker)) stop("show_ci_marker must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(facet_strata) || length(facet_strata) != 1 || is.na(facet_strata)) stop("facet_strata must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(cap) && (!is.numeric(cap) || length(cap) != 1 || is.na(cap) || cap <= 0)) stop("cap must be NULL or a single positive number.", call. = FALSE)

  is_lr_object <- is.list(x) && "Results" %in% names(x)
  Results <- if (is_lr_object) x$Results else x
  BinarySummary <- if (is_lr_object) x$BinarySummary else NULL
  required_columns <- c("Outcome", "OutcomeLabel", "Predictor", "PredictorLabel", "PredictorLevels", "ResultLevel", "Stratum", "N", "NCase", "NControl", "CaseWithResult", "ControlWithResult", "CaseWithoutResult", "ControlWithoutResult", "CaseProbability", "ControlProbability", "LikelihoodRatio", "LRLowerCI", "LRUpperCI", "Log2LR")
  if (!is.data.frame(Results) || length(setdiff(required_columns, names(Results)))) stop("x must be a DiagnosticLikelihoodRatioTable() result or a compatible Results data frame.", call. = FALSE)
  if (!is.character(na_color) || length(na_color) != 1 || is.na(na_color)) stop("na_color must be one color string.", call. = FALSE)

  # Prepare data

  if (result == "all") {
    plot_data <- Results %>% dplyr::mutate(RowId = paste(.data$Predictor, .data$ResultLevel, sep = "\r"), RowLabel = paste0(.data$PredictorLabel, ": ", .data$ResultLevel))
  } else {
    if (is.null(BinarySummary) || !nrow(BinarySummary)) stop("result = '", result, "' requires binary diagnostic predictors in a DiagnosticLikelihoodRatioTable() result.", call. = FALSE)
    selected_levels <- BinarySummary %>%
      dplyr::transmute(Outcome, Predictor, Stratum,
        SelectedResultLevel = if (result == "positive") .data$PositiveTestLevel else .data$NegativeTestLevel)
    plot_data <- Results %>%
      dplyr::inner_join(selected_levels, by = c("Outcome", "Predictor", "Stratum")) %>%
      dplyr::filter(.data$ResultLevel == .data$SelectedResultLevel) %>%
      dplyr::select(-tidyselect::all_of("SelectedResultLevel")) %>%
      dplyr::mutate(RowId = .data$Predictor, RowLabel = .data$PredictorLabel)
  }
  if (!nrow(plot_data)) stop("No diagnostic likelihood-ratio results are available to plot.", call. = FALSE)
  plot_data <- plot_data %>% dplyr::mutate(CIExcludesOne = !is.na(.data$LRLowerCI) & !is.na(.data$LRUpperCI) & (.data$LRLowerCI > 1 | .data$LRUpperCI < 1), ValueLabel = vapply(seq_len(n()), function(index) ScidrDiagnosticFormatLR(.data$LikelihoodRatio[[index]], .data$LRLowerCI[[index]], .data$LRUpperCI[[index]]), character(1)), ValueLabel = sub(" \\([^)]*\\)$", "", .data$ValueLabel), ValueLabel = if (show_ci_marker) ifelse(.data$CIExcludesOne, paste0(.data$ValueLabel, "*"), .data$ValueLabel) else .data$ValueLabel)
  plot_data$HoverText <- paste0("Predictor: ", plot_data$PredictorLabel, "<br>Result: ", plot_data$ResultLevel, "<br>Outcome: ", plot_data$OutcomeLabel, "<br>Stratum: ", plot_data$Stratum, "<br>LR: ", vapply(seq_len(nrow(plot_data)), function(index) ScidrDiagnosticFormatLR(plot_data$LikelihoodRatio[[index]], plot_data$LRLowerCI[[index]], plot_data$LRUpperCI[[index]]), character(1)), "<br>N: ", plot_data$N, " (cases: ", plot_data$NCase, "; controls: ", plot_data$NControl, ")<br>P(result | case): ", scales::percent(plot_data$CaseProbability, accuracy = 0.1), "<br>P(result | control): ", scales::percent(plot_data$ControlProbability, accuracy = 0.1), "<br>Cells [case result, control result, case other, control other]: ", paste(plot_data$CaseWithResult, plot_data$ControlWithResult, plot_data$CaseWithoutResult, plot_data$ControlWithoutResult, sep = ", "))

  # Order axes

  original_rows <- unique(plot_data$RowId)
  original_outcomes <- unique(plot_data$Outcome)
  row_order <- if (predictor_order == "alphabetical") plot_data %>% dplyr::distinct(.data$RowId, .data$RowLabel) %>% dplyr::arrange(.data$RowLabel) %>% dplyr::pull(.data$RowId) else if (predictor_order == "strength") plot_data %>% dplyr::group_by(.data$RowId) %>% dplyr::summarise(Strength = max(abs(.data$Log2LR[is.finite(.data$Log2LR)]), na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(Strength = dplyr::if_else(is.finite(.data$Strength), .data$Strength, -Inf)) %>% dplyr::arrange(dplyr::desc(.data$Strength)) %>% dplyr::pull(.data$RowId) else original_rows
  column_order <- if (outcome_order == "alphabetical") plot_data %>% dplyr::distinct(.data$Outcome, .data$OutcomeLabel) %>% dplyr::arrange(.data$OutcomeLabel) %>% dplyr::pull(.data$Outcome) else if (outcome_order == "strength") plot_data %>% dplyr::group_by(.data$Outcome) %>% dplyr::summarise(Strength = max(abs(.data$Log2LR[is.finite(.data$Log2LR)]), na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(Strength = dplyr::if_else(is.finite(.data$Strength), .data$Strength, -Inf)) %>% dplyr::arrange(dplyr::desc(.data$Strength)) %>% dplyr::pull(.data$Outcome) else original_outcomes
  if (predictor_order == "cluster" || outcome_order == "cluster") {
    AxisOrder <- .OrderHeatmapAxes(plot_data, "RowId", "Outcome", "Log2LR", row_order = row_order, column_order = column_order, cluster_rows = predictor_order == "cluster", cluster_columns = outcome_order == "cluster")
    row_order <- AxisOrder$rows; column_order <- AxisOrder$columns
  }
  row_labels <- plot_data %>% dplyr::distinct(.data$RowId, .data$RowLabel) %>% tibble::deframe()
  outcome_labels <- plot_data %>% dplyr::distinct(.data$Outcome, .data$OutcomeLabel) %>% tibble::deframe()
  finite_values <- abs(plot_data$Log2LR[is.finite(plot_data$Log2LR)])
  if (is.null(cap)) cap <- if (length(finite_values)) max(1, as.numeric(stats::quantile(finite_values, 0.95, na.rm = TRUE))) else 1
  plot_data <- plot_data %>% dplyr::mutate(PlotLog2LR = dplyr::case_when(is.infinite(.data$Log2LR) & .data$Log2LR > 0 ~ cap, is.infinite(.data$Log2LR) & .data$Log2LR < 0 ~ -cap, TRUE ~ pmax(-cap, pmin(cap, .data$Log2LR))), RowId = factor(.data$RowId, levels = rev(row_order)), Outcome = factor(.data$Outcome, levels = column_order))
  if (identical(show_values, "auto")) show_values <- nrow(plot_data) <= 150
  if (!is.logical(show_values) || length(show_values) != 1 || is.na(show_values)) stop("show_values must be 'auto', TRUE, or FALSE.", call. = FALSE)
  if (orientation == "auto") orientation <- if (length(row_order) >= length(column_order)) "predictors_rows" else "outcomes_rows"

  # Build plot

  x_axis <- if (orientation == "predictors_rows") "Outcome" else "RowId"
  y_axis <- if (orientation == "predictors_rows") "RowId" else "Outcome"
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]], fill = .data$PlotLog2LR, text = .data$HoverText)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    .GetHeatmapColorScale(low_color = SciDataPalette()[["Navy"]], mid_color = "white", high_color = SciDataPalette()[["Orange"]], fill_midpoint = 0, fill_limits = c(-cap, cap), fill_oob = scales::squish, name = "Diagnostic\nlikelihood ratio") +
    ggplot2::scale_x_discrete(labels = if (orientation == "predictors_rows") outcome_labels else row_labels) +
    ggplot2::scale_y_discrete(labels = if (orientation == "predictors_rows") row_labels else outcome_labels) +
    ggplot2::labs(x = NULL, y = NULL, caption = if (show_ci_marker) "* Unadjusted likelihood-ratio confidence interval excludes 1." else NULL) +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), strip.text = ggplot2::element_text(face = "bold"))
  if (show_values) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$ValueLabel), size = 3)
  if (facet_strata && dplyr::n_distinct(plot_data$Stratum) > 1) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$Stratum))
  attr(p, "DiagnosticLRData") <- plot_data
  attr(p, "DiagnosticLRCap") <- cap
  attr(p, "DiagnosticLROrientation") <- orientation
  p
}
