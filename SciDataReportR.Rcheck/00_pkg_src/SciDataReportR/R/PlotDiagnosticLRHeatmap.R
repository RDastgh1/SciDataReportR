#' Plot a diagnostic likelihood-ratio heatmap
#'
#' Visualizes results created by [DiagnosticLikelihoodRatioTable()] without
#' recomputing statistics. Tile fill is log2(LR), so reciprocal evidence is
#' equally distant from LR 1. Multi-outcome input can return a compact screening
#' overview alongside outcome-specific diagnostic panels.
#'
#' @param x An object returned by [DiagnosticLikelihoodRatioTable()], or a tidy
#'   data frame compatible with its `Results` element.
#' @param result Diagnostic result levels to display: `"all"`, `"positive"`, or
#'   `"negative"`. Positive and negative selections require the result object.
#' @param predictor_order Predictor ordering: `"original"`, `"alphabetical"`,
#'   `"strength"`, or `"cluster"`.
#' @param outcome_order Outcome ordering with the same choices.
#' @param orientation Tile orientation: `"auto"`, `"predictors_rows"`, or
#'   `"outcomes_rows"`.
#' @param show_values `"auto"`, `TRUE`, or `FALSE`; auto labels at most 150 tiles.
#' @param show_ci_marker Logical; append `*` when the unadjusted LR CI excludes 1.
#' @param facet_strata Logical; facet separate diagnostic strata when present.
#' @param cap Optional positive maximum absolute log2(LR) for color scaling.
#' @param na_color Fill color for unavailable LR values.
#' @param multi_outcome Multi-outcome display: `"auto"` returns a split
#'   overview/panel list for multiple outcomes, `"combined"` returns one
#'   all-results matrix, and `"split"` always returns the linked plot list.
#'
#' @return A named list. Single-outcome and combined displays contain
#'   `DiagnosticLR`, `DiagnosticMatrices`, `DiagnosticLRData`, and
#'   `DiagnosticMatrixData`. Split multi-outcome displays contain `Overview`,
#'   `ByOutcome`, `DiagnosticMatrices`, and their corresponding tidy data.
#'
#' @seealso [DiagnosticLikelihoodRatioTable()] to calculate displayed LRs.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' df_Labelled$DiagnosisBinary <- factor(df_Labelled$Diagnosis,
#'   levels = c("Control", "Impaired"))
#' lr <- DiagnosticLikelihoodRatioTable(df_Labelled, "DiagnosisBinary",
#'   c("sex", "Genotype"))
#' plots <- PlotDiagnosticLRHeatmap(lr)
#' plots$DiagnosticLR
#' plots$DiagnosticMatrices
#' @export
PlotDiagnosticLRHeatmap <- function(x, result = c("all", "positive", "negative"),
    predictor_order = c("original", "alphabetical", "strength", "cluster"),
    outcome_order = c("original", "alphabetical", "strength", "cluster"),
    orientation = c("auto", "predictors_rows", "outcomes_rows"),
    show_values = "auto", show_ci_marker = TRUE, facet_strata = TRUE,
    cap = NULL, na_color = "grey90",
    multi_outcome = c("auto", "combined", "split")) {
  result <- match.arg(result)
  predictor_order <- match.arg(predictor_order)
  outcome_order <- match.arg(outcome_order)
  orientation <- match.arg(orientation)
  multi_outcome <- match.arg(multi_outcome)
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
    selected_levels <- BinarySummary %>% dplyr::transmute(Outcome, Predictor, Stratum, SelectedResultLevel = if (result == "positive") .data$PositiveTestLevel else .data$NegativeTestLevel)
    plot_data <- Results %>% dplyr::inner_join(selected_levels, by = c("Outcome", "Predictor", "Stratum")) %>% dplyr::filter(.data$ResultLevel == .data$SelectedResultLevel) %>% dplyr::select(-tidyselect::all_of("SelectedResultLevel")) %>% dplyr::mutate(RowId = .data$Predictor, RowLabel = .data$PredictorLabel)
  }
  if (!nrow(plot_data)) stop("No diagnostic likelihood-ratio results are available to plot.", call. = FALSE)
  plot_data <- ScidrDiagnosticPrepareHeatmapData(plot_data, show_ci_marker)
  matrix_data <- ScidrDiagnosticMatrixData(Results)
  finite_values <- abs(plot_data$Log2LR[is.finite(plot_data$Log2LR)])
  if (is.null(cap)) cap <- if (length(finite_values)) max(1, as.numeric(stats::quantile(finite_values, 0.95, na.rm = TRUE))) else 1
  plot_data <- ScidrDiagnosticCapHeatmapData(plot_data, cap)
  if (identical(show_values, "auto")) show_values <- nrow(plot_data) <= 150
  if (!is.logical(show_values) || length(show_values) != 1 || is.na(show_values)) stop("show_values must be 'auto', TRUE, or FALSE.", call. = FALSE)

  n_outcomes <- dplyr::n_distinct(plot_data$Outcome)
  split_display <- multi_outcome == "split" || (multi_outcome == "auto" && n_outcomes > 1)
  if (!split_display) {
    diagnostic_lr <- ScidrDiagnosticBuildHeatmap(
      data = plot_data, predictor_order = predictor_order,
      outcome_order = outcome_order, orientation = orientation,
      show_values = show_values, show_ci_marker = show_ci_marker,
      facet_strata = facet_strata, cap = cap, na_color = na_color
    )
    diagnostic_matrices <- ScidrDiagnosticBuildMatrixPlot(matrix_data)
    return(list(
      DiagnosticLR = diagnostic_lr,
      DiagnosticMatrices = diagnostic_matrices,
      DiagnosticLRData = plot_data,
      DiagnosticMatrixData = matrix_data
    ))
  }

  # Build split multi-outcome outputs

  overview_data <- ScidrDiagnosticOverviewData(plot_data)
  overview <- ScidrDiagnosticBuildHeatmap(
    data = overview_data, predictor_order = predictor_order,
    outcome_order = outcome_order, orientation = orientation,
    show_values = FALSE, show_ci_marker = show_ci_marker,
    facet_strata = facet_strata, cap = cap, na_color = na_color,
    overview = TRUE
  )
  by_outcome <- ScidrDiagnosticBuildHeatmap(
    data = plot_data, predictor_order = predictor_order,
    outcome_order = outcome_order, orientation = "predictors_rows",
    show_values = show_values, show_ci_marker = show_ci_marker,
    facet_strata = facet_strata, cap = cap, na_color = na_color,
    facet_outcome = TRUE
  )
  diagnostic_matrices <- ScidrDiagnosticBuildMatrixPlot(matrix_data)
  out <- list(Overview = overview, ByOutcome = by_outcome,
    DiagnosticMatrices = diagnostic_matrices, OverviewData = overview_data,
    ByOutcomeData = plot_data, DiagnosticMatrixData = matrix_data)
  attr(out, "DiagnosticLRCap") <- cap
  attr(out, "DiagnosticLROrientation") <- orientation
  out
}

ScidrDiagnosticPrepareHeatmapData <- function(data, show_ci_marker) {
  data <- data %>% dplyr::mutate(
    .OriginalOrder = dplyr::row_number(),
    CIExcludesOne = !is.na(.data$LRLowerCI) & !is.na(.data$LRUpperCI) &
      (.data$LRLowerCI > 1 | .data$LRUpperCI < 1),
    ValueLabel = vapply(seq_len(dplyr::n()), function(index) {
      ScidrDiagnosticFormatLR(.data$LikelihoodRatio[[index]], .data$LRLowerCI[[index]], .data$LRUpperCI[[index]])
    }, character(1))
  ) %>% dplyr::mutate(
    ValueLabel = sub(" \\([^)]*\\)$", "", .data$ValueLabel),
    ValueLabel = if (show_ci_marker) ifelse(.data$CIExcludesOne,
      paste0(.data$ValueLabel, "*"), .data$ValueLabel) else .data$ValueLabel
  )
  data$HoverText <- paste0(
    "Predictor: ", data$PredictorLabel, "<br>Result: ", data$ResultLevel,
    "<br>Outcome: ", data$OutcomeLabel, "<br>Stratum: ", data$Stratum,
    "<br>LR: ", vapply(seq_len(nrow(data)), function(index) {
      ScidrDiagnosticFormatLR(data$LikelihoodRatio[[index]], data$LRLowerCI[[index]], data$LRUpperCI[[index]])
    }, character(1)), "<br>N: ", data$N, " (cases: ", data$NCase,
    "; controls: ", data$NControl, ")<br>P(result | case): ",
    scales::percent(data$CaseProbability, accuracy = 0.1),
    "<br>P(result | control): ", scales::percent(data$ControlProbability, accuracy = 0.1),
    "<br>Cells [case result, control result, case other, control other]: ",
    paste(data$CaseWithResult, data$ControlWithResult, data$CaseWithoutResult,
      data$ControlWithoutResult, sep = ", ")
  )
  data
}

ScidrDiagnosticCapHeatmapData <- function(data, cap) {
  data %>% dplyr::mutate(PlotLog2LR = dplyr::case_when(
    is.infinite(.data$Log2LR) & .data$Log2LR > 0 ~ cap,
    is.infinite(.data$Log2LR) & .data$Log2LR < 0 ~ -cap,
    TRUE ~ pmax(-cap, pmin(cap, .data$Log2LR))
  ))
}

ScidrDiagnosticMatrixData <- function(results) {
  control_data <- results %>% dplyr::transmute(
    Stratum, Predictor, PredictorLabel, Outcome, OutcomeLabel,
    ResultLevel, OutcomeCondition = "Outcome negative",
    Count = .data$ControlWithResult, Denominator = .data$NControl,
    Probability = .data$ControlProbability
  )
  case_data <- results %>% dplyr::transmute(
    Stratum, Predictor, PredictorLabel, Outcome, OutcomeLabel,
    ResultLevel, OutcomeCondition = "Outcome positive",
    Count = .data$CaseWithResult, Denominator = .data$NCase,
    Probability = .data$CaseProbability
  )
  dplyr::bind_rows(control_data, case_data) %>% dplyr::mutate(
    CellLabel = paste0(.data$Count, "/", .data$Denominator, " (",
      scales::percent(.data$Probability, accuracy = 0.1), ")"),
    MatrixRow = paste0(.data$Predictor, "\r", .data$ResultLevel),
    MatrixRowLabel = .data$ResultLevel,
    HoverText = paste0(
      "Predictor: ", .data$PredictorLabel, "<br>Result: ", .data$ResultLevel,
      "<br>Outcome: ", .data$OutcomeLabel, "<br>Stratum: ", .data$Stratum,
      "<br>", .data$OutcomeCondition, ": ", .data$CellLabel
    )
  )
}

ScidrDiagnosticBuildMatrixPlot <- function(data) {
  row_labels <- data %>%
    dplyr::distinct(.data$MatrixRow, .data$MatrixRowLabel) %>%
    tibble::deframe()
  row_order <- data %>%
    dplyr::distinct(.data$Predictor, .data$MatrixRow) %>%
    dplyr::pull(.data$MatrixRow)
  plot_data <- data %>% dplyr::mutate(
    MatrixRow = factor(.data$MatrixRow, levels = rev(row_order)),
    OutcomeCondition = factor(.data$OutcomeCondition,
      levels = c("Outcome negative", "Outcome positive"))
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$OutcomeCondition, y = .data$MatrixRow,
    fill = .data$Probability, text = .data$HoverText
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = .data$CellLabel), size = 3) +
    ggplot2::scale_fill_gradient(
      low = "white", high = SciDataPalette()[["Navy"]],
      limits = c(0, 1), oob = scales::squish, na.value = "grey90",
      labels = scales::percent_format(accuracy = 1),
      name = "Within-outcome\nresult probability"
    ) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$Stratum, .data$PredictorLabel),
      cols = ggplot2::vars(.data$OutcomeLabel), scales = "free_y", space = "free_y"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(face = "bold")
    )
  attr(p, "DiagnosticMatrixData") <- plot_data
  p
}

ScidrDiagnosticOverviewData <- function(data) {
  data %>% dplyr::group_by(.data$Stratum, .data$Outcome, .data$Predictor) %>%
    dplyr::arrange(dplyr::desc(is.finite(.data$Log2LR)),
      dplyr::desc(abs(dplyr::if_else(is.finite(.data$Log2LR), .data$Log2LR, -Inf))),
      .data$.OriginalOrder, .by_group = TRUE) %>%
    dplyr::slice(1) %>% dplyr::ungroup() %>%
    dplyr::mutate(RowId = .data$Predictor, RowLabel = .data$PredictorLabel,
      HoverText = paste0("Selected strongest finite result<br>", .data$HoverText))
}

ScidrDiagnosticBuildHeatmap <- function(data, predictor_order, outcome_order,
    orientation, show_values, show_ci_marker, facet_strata, cap, na_color,
    overview = FALSE, facet_outcome = FALSE) {
  original_rows <- unique(data$RowId)
  original_outcomes <- unique(data$Outcome)
  row_order <- ScidrDiagnosticAxisOrder(data, "RowId", "RowLabel", predictor_order, original_rows)
  column_order <- ScidrDiagnosticAxisOrder(data, "Outcome", "OutcomeLabel", outcome_order, original_outcomes)
  if (predictor_order == "cluster" || outcome_order == "cluster") {
    AxisOrder <- .OrderHeatmapAxes(data, "RowId", "Outcome", "Log2LR",
      row_order = row_order, column_order = column_order,
      cluster_rows = predictor_order == "cluster",
      cluster_columns = outcome_order == "cluster")
    row_order <- AxisOrder$rows
    column_order <- AxisOrder$columns
  }
  row_labels <- data %>% dplyr::distinct(.data$RowId, .data$RowLabel) %>% tibble::deframe()
  outcome_labels <- data %>% dplyr::distinct(.data$Outcome, .data$OutcomeLabel) %>% tibble::deframe()
  data <- data %>% dplyr::mutate(
    RowId = factor(.data$RowId, levels = rev(row_order)),
    Outcome = factor(.data$Outcome, levels = column_order)
  )
  if (orientation == "auto") orientation <- if (length(row_order) >= length(column_order)) "predictors_rows" else "outcomes_rows"
  x_axis <- if (orientation == "predictors_rows") "Outcome" else "RowId"
  y_axis <- if (orientation == "predictors_rows") "RowId" else "Outcome"
  legend_breaks <- c(-cap, 0, cap)
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]],
      fill = .data$PlotLog2LR, text = .data$HoverText)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(low = SciDataPalette()[["Navy"]], mid = "white",
      high = SciDataPalette()[["Orange"]], midpoint = 0, limits = c(-cap, cap),
      oob = scales::squish, na.value = na_color, breaks = legend_breaks,
      labels = paste0("LR ", trimws(formatC(2^legend_breaks, digits = 2, format = "fg"))),
      name = "Diagnostic\nlikelihood ratio") +
    ggplot2::scale_x_discrete(labels = if (facet_outcome) NULL else if (orientation == "predictors_rows") outcome_labels else row_labels) +
    ggplot2::scale_y_discrete(labels = if (orientation == "predictors_rows") row_labels else outcome_labels) +
    ggplot2::labs(x = NULL, y = NULL,
      caption = if (show_ci_marker) "* Unadjusted likelihood-ratio confidence interval excludes 1." else NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = if (facet_outcome) 0 else 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold"))
  if (show_values) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$ValueLabel), size = 3)
  if (facet_outcome && facet_strata && dplyr::n_distinct(data$Stratum) > 1) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$Stratum), cols = ggplot2::vars(.data$OutcomeLabel), scales = "free_y")
  } else if (facet_outcome) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$OutcomeLabel), scales = "free_y")
  } else if (facet_strata && dplyr::n_distinct(data$Stratum) > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$Stratum))
  }
  attr(p, "DiagnosticLRData") <- data
  attr(p, "DiagnosticLRCap") <- cap
  attr(p, "DiagnosticLROrientation") <- orientation
  attr(p, "DiagnosticLROverview") <- overview
  p
}

ScidrDiagnosticAxisOrder <- function(data, id_col, label_col, order, original_order) {
  if (order == "alphabetical") return(data %>% dplyr::distinct(.data[[id_col]], .data[[label_col]]) %>% dplyr::arrange(.data[[label_col]]) %>% dplyr::pull(.data[[id_col]]))
  if (order == "strength") return(data %>% dplyr::group_by(.data[[id_col]]) %>% dplyr::summarise(Strength = max(abs(.data$Log2LR[is.finite(.data$Log2LR)]), na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(Strength = dplyr::if_else(is.finite(.data$Strength), .data$Strength, -Inf)) %>% dplyr::arrange(dplyr::desc(.data$Strength)) %>% dplyr::pull(.data[[id_col]]))
  original_order
}
