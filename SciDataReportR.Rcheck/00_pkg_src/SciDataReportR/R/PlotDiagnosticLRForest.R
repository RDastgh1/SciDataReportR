#' Plot diagnostic likelihood ratios as a forest plot
#'
#' Visualizes likelihood ratios calculated by [DiagnosticLikelihoodRatioTable()]
#' without recalculating diagnostic statistics. Each point is a diagnostic test
#' result and its horizontal interval is the likelihood-ratio confidence
#' interval. The logarithmic scale makes reciprocal likelihood ratios equally
#' distant from the neutral value of one.
#'
#' @param x An object returned by [DiagnosticLikelihoodRatioTable()], or a tidy
#'   data frame compatible with its `Results` element.
#' @param result Diagnostic result levels to display: `"all"`, `"positive"`, or
#'   `"negative"`. Positive and negative selections require the result object.
#' @param predictor_order Predictor ordering: `"original"`, `"alphabetical"`,
#'   or `"strength"`. Strength orders predictors by their largest finite
#'   absolute log2 likelihood ratio.
#' @param outcome_order Outcome ordering with the same choices.
#' @param facet_by Plot panels by `"outcome"` (the default) or `"predictor"`.
#' @param facet_strata Logical; add strata as facet row groups when present.
#' @param limits Optional numeric vector of two positive, increasing likelihood
#'   ratio limits. By default, limits are chosen from finite estimates and
#'   confidence intervals while always including one.
#' @param p_size Numeric point size.
#'
#' @return A ggplot object. Its `DiagnosticLRForestData`,
#'   `DiagnosticLRForestLimits`, and `DiagnosticLRForestFacetBy` attributes
#'   retain the prepared data and resolved plotting settings.
#'
#' @section Reading the plot:
#' The dashed vertical line marks LR = 1. Dark black estimates have an
#' unadjusted likelihood-ratio confidence interval that excludes one; gray
#' estimates do not. This visual cue is not an FDR-adjusted significance test.
#' Uncorrected zero and infinite likelihood ratios are shown as boundary arrows
#' labelled `0` and `Inf`, because a finite confidence interval is unavailable.
#'
#' @seealso [DiagnosticLikelihoodRatioTable()] to calculate diagnostic LRs, and
#'   [PlotDiagnosticLRHeatmap()] for matrix and count-matrix views.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' df_Labelled$DiagnosisBinary <- factor(df_Labelled$Diagnosis,
#'   levels = c("Control", "Impaired"))
#' lr <- DiagnosticLikelihoodRatioTable(df_Labelled, "DiagnosisBinary",
#'   c("sex", "Genotype"))
#' PlotDiagnosticLRForest(lr)
#' PlotDiagnosticLRForest(lr, result = "positive")
#' @export
PlotDiagnosticLRForest <- function(x, result = c("all", "positive", "negative"),
    predictor_order = c("original", "alphabetical", "strength"),
    outcome_order = c("original", "alphabetical", "strength"),
    facet_by = c("outcome", "predictor"), facet_strata = TRUE,
    limits = NULL, p_size = 2) {

  result <- match.arg(result)
  predictor_order <- match.arg(predictor_order)
  outcome_order <- match.arg(outcome_order)
  facet_by <- match.arg(facet_by)
  if (!is.logical(facet_strata) || length(facet_strata) != 1 || is.na(facet_strata)) stop("facet_strata must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(p_size) || length(p_size) != 1 || is.na(p_size) || p_size <= 0) stop("p_size must be one positive number.", call. = FALSE)
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2 || anyNA(limits) || any(!is.finite(limits)) || any(limits <= 0) || limits[[1]] >= limits[[2]])) stop("limits must be NULL or two increasing positive finite likelihood-ratio values.", call. = FALSE)

  # Prepare data

  is_lr_object <- is.list(x) && "Results" %in% names(x)
  Results <- if (is_lr_object) x$Results else x
  BinarySummary <- if (is_lr_object) x$BinarySummary else NULL
  required_columns <- c("Outcome", "OutcomeLabel", "Predictor", "PredictorLabel", "ResultLevel", "Stratum", "LikelihoodRatio", "LRLowerCI", "LRUpperCI", "Log2LR")
  if (!is.data.frame(Results) || length(setdiff(required_columns, names(Results)))) stop("x must be a DiagnosticLikelihoodRatioTable() result or a compatible Results data frame.", call. = FALSE)

  if (result == "all") {
    df_Plot <- Results
  } else {
    if (is.null(BinarySummary) || !nrow(BinarySummary)) stop("result = '", result, "' requires binary diagnostic predictors in a DiagnosticLikelihoodRatioTable() result.", call. = FALSE)
    df_SelectedLevels <- BinarySummary %>%
      dplyr::transmute(
        Outcome, Predictor, Stratum,
        SelectedResultLevel = if (result == "positive") .data$PositiveTestLevel else .data$NegativeTestLevel
      )
    df_Plot <- Results %>%
      dplyr::inner_join(df_SelectedLevels, by = c("Outcome", "Predictor", "Stratum")) %>%
      dplyr::filter(.data$ResultLevel == .data$SelectedResultLevel) %>%
      dplyr::select(-tidyselect::all_of("SelectedResultLevel"))
  }
  df_Plot <- df_Plot %>%
    dplyr::filter(!is.na(.data$LikelihoodRatio)) %>%
    dplyr::mutate(
      EstimateStatus = dplyr::case_when(
        .data$LikelihoodRatio == 0 ~ "Zero",
        is.infinite(.data$LikelihoodRatio) ~ "Infinite",
        is.finite(.data$LikelihoodRatio) & .data$LikelihoodRatio > 0 ~ "Finite",
        TRUE ~ "Unavailable"
      ),
      CIExcludesOne = is.finite(.data$LRLowerCI) & is.finite(.data$LRUpperCI) &
        (.data$LRUpperCI < 1 | .data$LRLowerCI > 1),
      CIHighlight = ifelse(.data$CIExcludesOne, "CI excludes 1", "CI includes 1 or unavailable"),
      HoverText = paste0(
        "Predictor: ", .data$PredictorLabel,
        "<br>Result: ", .data$ResultLevel,
        "<br>Outcome: ", .data$OutcomeLabel,
        "<br>Stratum: ", .data$Stratum,
        "<br>LR: ", mapply(ScidrDiagnosticFormatLR, .data$LikelihoodRatio,
          .data$LRLowerCI, .data$LRUpperCI)
      )
    )
  if (!nrow(df_Plot)) stop("No diagnostic likelihood-ratio estimates are available to plot.", call. = FALSE)

  # Resolve limits and plotting coordinates

  plot_limits <- ScidrDiagnosticForestLimits(df_Plot, limits)
  lower_limit <- plot_limits[[1]]
  upper_limit <- plot_limits[[2]]
  df_Plot <- df_Plot %>%
    dplyr::mutate(
      PlotEstimate = dplyr::case_when(
        .data$EstimateStatus == "Zero" ~ lower_limit,
        .data$EstimateStatus == "Infinite" ~ upper_limit,
        TRUE ~ pmax(lower_limit, pmin(upper_limit, .data$LikelihoodRatio))
      ),
      PlotLowerCI = dplyr::if_else(
        is.finite(.data$LRLowerCI) & .data$LRLowerCI > 0,
        pmax(lower_limit, pmin(upper_limit, .data$LRLowerCI)), NA_real_
      ),
      PlotUpperCI = dplyr::if_else(
        is.finite(.data$LRUpperCI) & .data$LRUpperCI > 0,
        pmax(lower_limit, pmin(upper_limit, .data$LRUpperCI)), NA_real_
      ),
      ArrowLeft = .data$EstimateStatus == "Zero" |
        (is.finite(.data$LRLowerCI) & .data$LRLowerCI < lower_limit),
      ArrowRight = .data$EstimateStatus == "Infinite" |
        (is.finite(.data$LRUpperCI) & .data$LRUpperCI > upper_limit),
      BoundaryLabel = dplyr::case_when(
        .data$EstimateStatus == "Zero" ~ "0",
        .data$EstimateStatus == "Infinite" ~ "Inf",
        TRUE ~ NA_character_
      )
    )

  # Order rows and facets

  predictor_levels <- ScidrDiagnosticForestOrder(df_Plot, "Predictor", "PredictorLabel", predictor_order)
  outcome_levels <- ScidrDiagnosticForestOrder(df_Plot, "Outcome", "OutcomeLabel", outcome_order)
  if (facet_by == "outcome") {
    df_RowLevels <- df_Plot %>%
      dplyr::mutate(Predictor = factor(.data$Predictor, levels = predictor_levels)) %>%
      dplyr::arrange(.data$Predictor) %>%
      dplyr::distinct(.data$Predictor, .data$ResultLevel, .keep_all = TRUE) %>%
      dplyr::mutate(PlotRowId = paste(.data$Predictor, .data$ResultLevel, sep = "\r"), PlotRowLabel = paste0(.data$PredictorLabel, ": ", .data$ResultLevel))
    row_levels <- df_RowLevels$PlotRowId
    row_labels <- stats::setNames(df_RowLevels$PlotRowLabel, df_RowLevels$PlotRowId)
    df_Plot <- df_Plot %>%
      dplyr::mutate(
        PlotRowId = paste(.data$Predictor, .data$ResultLevel, sep = "\r"),
        PlotFacet = factor(.data$OutcomeLabel, levels = (df_Plot %>% dplyr::distinct(.data$Outcome, .data$OutcomeLabel) %>% dplyr::mutate(Outcome = factor(.data$Outcome, levels = outcome_levels)) %>% dplyr::arrange(.data$Outcome) %>% dplyr::pull(.data$OutcomeLabel)))
      )
  } else {
    df_RowLevels <- df_Plot %>%
      dplyr::mutate(Outcome = factor(.data$Outcome, levels = outcome_levels)) %>%
      dplyr::arrange(.data$Outcome) %>%
      dplyr::distinct(.data$Outcome, .data$ResultLevel, .keep_all = TRUE) %>%
      dplyr::mutate(PlotRowId = paste(.data$Outcome, .data$ResultLevel, sep = "\r"), PlotRowLabel = paste0(.data$OutcomeLabel, ": ", .data$ResultLevel))
    row_levels <- df_RowLevels$PlotRowId
    row_labels <- stats::setNames(df_RowLevels$PlotRowLabel, df_RowLevels$PlotRowId)
    df_Plot <- df_Plot %>%
      dplyr::mutate(
        PlotRowId = paste(.data$Outcome, .data$ResultLevel, sep = "\r"),
        PlotFacet = factor(.data$PredictorLabel, levels = (df_Plot %>% dplyr::distinct(.data$Predictor, .data$PredictorLabel) %>% dplyr::mutate(Predictor = factor(.data$Predictor, levels = predictor_levels)) %>% dplyr::arrange(.data$Predictor) %>% dplyr::pull(.data$PredictorLabel)))
      )
  }
  df_Plot <- df_Plot %>% dplyr::mutate(PlotRow = factor(.data$PlotRowId, levels = rev(row_levels)))

  # Build plot

  p <- ggplot2::ggplot(df_Plot, ggplot2::aes(
    x = .data$PlotEstimate, y = .data$PlotRow, color = .data$CIHighlight,
    text = .data$HoverText
  )) +
    ggplot2::geom_errorbar(
      data = df_Plot %>% dplyr::filter(.data$EstimateStatus == "Finite", !is.na(.data$PlotLowerCI), !is.na(.data$PlotUpperCI)),
      ggplot2::aes(xmin = .data$PlotLowerCI, xmax = .data$PlotUpperCI),
      orientation = "y", width = 0.2, linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = df_Plot %>% dplyr::filter(.data$EstimateStatus == "Finite"),
      size = p_size
    ) +
    ggplot2::geom_segment(
      data = df_Plot %>% dplyr::filter(.data$ArrowLeft),
      ggplot2::aes(x = pmin(upper_limit, lower_limit * 1.12), xend = lower_limit, y = .data$PlotRow, yend = .data$PlotRow),
      arrow = grid::arrow(type = "closed", length = grid::unit(0.10, "inches")), linewidth = 0.45
    ) +
    ggplot2::geom_segment(
      data = df_Plot %>% dplyr::filter(.data$ArrowRight),
      ggplot2::aes(x = pmax(lower_limit, upper_limit / 1.12), xend = upper_limit, y = .data$PlotRow, yend = .data$PlotRow),
      arrow = grid::arrow(type = "closed", length = grid::unit(0.10, "inches")), linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = df_Plot %>% dplyr::filter(!is.na(.data$BoundaryLabel)),
      ggplot2::aes(label = .data$BoundaryLabel), color = "black", size = 3,
      vjust = -0.8, show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35") +
    ggplot2::scale_x_log10(limits = plot_limits, labels = scales::label_number()) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::scale_color_manual(values = c("CI excludes 1" = "black", "CI includes 1 or unavailable" = "grey55")) +
    ggplot2::labs(
      x = "Diagnostic likelihood ratio", y = NULL,
      caption = "Dark black estimates have an unadjusted likelihood-ratio confidence interval that excludes 1; this is not FDR-adjusted."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold")
    )
  if (facet_strata && dplyr::n_distinct(df_Plot$Stratum) > 1) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$Stratum), cols = ggplot2::vars(.data$PlotFacet), scales = "free_y", space = "free_y")
  } else {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$PlotFacet), scales = "free_y")
  }
  attr(p, "DiagnosticLRForestData") <- df_Plot
  attr(p, "DiagnosticLRForestLimits") <- plot_limits
  attr(p, "DiagnosticLRForestFacetBy") <- facet_by
  p
}

ScidrDiagnosticForestLimits <- function(data, limits = NULL) {
  if (!is.null(limits)) return(limits)
  values <- c(data$LikelihoodRatio, data$LRLowerCI, data$LRUpperCI)
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) return(c(0.25, 4))
  lower <- min(c(values, 1))
  upper <- max(c(values, 1))
  c(10^floor(log10(lower / 1.15)), 10^ceiling(log10(upper * 1.15)))
}

ScidrDiagnosticForestOrder <- function(data, id_column, label_column, order) {
  df_Original <- data %>% dplyr::distinct(.data[[id_column]], .data[[label_column]])
  if (order == "original") return(df_Original %>% dplyr::pull(.data[[id_column]]))
  if (order == "alphabetical") return(df_Original %>% dplyr::arrange(.data[[label_column]]) %>% dplyr::pull(.data[[id_column]]))
  data %>%
    dplyr::group_by(.data[[id_column]]) %>%
    dplyr::summarise(Strength = if (any(is.finite(.data$Log2LR))) max(abs(.data$Log2LR[is.finite(.data$Log2LR)])) else -Inf, .groups = "drop") %>%
    dplyr::left_join(df_Original, by = id_column) %>%
    dplyr::arrange(dplyr::desc(.data$Strength), .data[[label_column]]) %>%
    dplyr::pull(.data[[id_column]])
}
