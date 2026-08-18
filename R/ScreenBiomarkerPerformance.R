#' Screen biomarker performance across outcomes
#'
#' Applies [EvaluateBiomarkerPerformance()] across many candidate biomarkers and
#' one or more binary or continuous outcomes. Each biomarker-outcome pair uses
#' all complete observations available for that pair and the requested
#' covariates. The function returns comparison tables and screening plots,
#' including an interactive-ready heatmap whose cells contain hover text.
#'
#' @param data A data frame.
#' @param outcome_vars Character vector of outcome variable names.
#' @param biomarker_vars Character vector of biomarker variable names.
#' @param covariates Optional character vector of covariate variable names.
#' @param PositiveLevel Optional positive level. Supply one value for all binary
#'   outcomes or a named character vector with names matching `outcome_vars`.
#' @param OutcomeType One of `"auto"`, `"binary"`, or `"continuous"`, applied
#'   to all outcomes unless auto-detection is used.
#' @param Validation One of `"none"`, `"bootstrap"`, or
#'   `"cross_validation"`.
#' @param BootstrapR Number of bootstrap resamples for internal validation.
#' @param CVFolds Number of cross-validation folds.
#' @param CIBootstrapR Number of bootstrap resamples for performance confidence
#'   intervals.
#' @param CILevel Confidence level. Default is `0.95`.
#' @param HeatmapMetric Performance metric shown by the heatmap fill. Default is
#'   `"AdjustedAUC"`. Common alternatives are `"AUC"`, `"DeltaAUC"`,
#'   `"AdjustedR2"`, and `"DeltaR2"`.
#' @param Seed Random seed used for resampling.
#' @param Relabel Logical indicating whether labels should be used when available.
#' @param codebook Optional data frame with columns `Variable` and `Label`.
#'
#' @return A named list with `PerformanceTable`, `RegressionTable`,
#'   `ThresholdTable`, `FailureTable`, `Evaluations`, `Plots`, and `Metadata`.
#'   For binary outcomes, `Plots` also includes `BiomarkerPanels` (raw
#'   outcome-stratified distributions for continuous biomarkers) and `ROCFacets`
#'   (adjusted ROC curves), each annotated with pair-specific performance.
#'
#' @examples
#' \dontrun{
#' screen <- ScreenBiomarkerPerformance(
#'   data = df,
#'   outcome_vars = c("DiseaseCohort", "Progression"),
#'   biomarker_vars = c("NfL", "GFAP", "Acrocyanosis"),
#'   covariates = c("Age", "Sex")
#' )
#'
#' screen$PerformanceTable
#' screen$Plots$Heatmap
#' screen$Plots$Heatmap %>% add_biomarker_values()
#' }
#' @export
ScreenBiomarkerPerformance <- function(
    data,
    outcome_vars,
    biomarker_vars,
    covariates = NULL,
    PositiveLevel = NULL,
    OutcomeType = c("auto", "binary", "continuous"),
    Validation = c("none", "bootstrap", "cross_validation"),
    BootstrapR = 500,
    CVFolds = 10,
    CIBootstrapR = 200,
    CILevel = 0.95,
    HeatmapMetric = "AdjustedAUC",
    Seed = 123,
    Relabel = TRUE,
    codebook = NULL) {

  OutcomeType <- match.arg(OutcomeType)
  Validation <- match.arg(Validation)

  # Validate inputs

  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (!is.character(outcome_vars) || length(outcome_vars) == 0) {
    stop("outcome_vars must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.character(biomarker_vars) || length(biomarker_vars) == 0) {
    stop("biomarker_vars must be a non-empty character vector.", call. = FALSE)
  }

  missing_outcomes <- setdiff(outcome_vars, names(data))
  missing_biomarkers <- setdiff(biomarker_vars, names(data))
  missing_covars <- setdiff(covariates, names(data))

  if (length(missing_outcomes) > 0) {
    stop(
      "The following outcome_vars were not found in data: ",
      paste(missing_outcomes, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_covars) > 0) {
    stop(
      "The following covariates were not found in data: ",
      paste(missing_covars, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_biomarkers) > 0) {
    warning(
      "The following biomarker_vars were not found in data and were removed: ",
      paste(missing_biomarkers, collapse = ", "),
      call. = FALSE
    )
    biomarker_vars <- intersect(biomarker_vars, names(data))
  }

  if (length(biomarker_vars) == 0) {
    stop("No valid biomarker_vars remain.", call. = FALSE)
  }

  positive_for_outcome <- function(outcome) {
    if (is.null(PositiveLevel)) return(NULL)
    if (length(PositiveLevel) == 1 && is.null(names(PositiveLevel))) {
      return(as.character(PositiveLevel))
    }
    if (!is.null(names(PositiveLevel)) && outcome %in% names(PositiveLevel)) {
      return(as.character(PositiveLevel[[outcome]]))
    }
    NULL
  }

  # Run evaluations

  evaluation_list <- list()
  performance_rows <- list()
  regression_rows <- list()
  threshold_rows <- list()
  inferred_positive <- list()
  failure_rows <- list()

  for (outcome in outcome_vars) {
    for (biomarker in biomarker_vars) {
      key <- paste(outcome, biomarker, sep = " :: ")

      evaluation <- tryCatch(
        EvaluateBiomarkerPerformance(
          data = data,
          outcome_var = outcome,
          biomarker_var = biomarker,
          covariates = covariates,
          PositiveLevel = positive_for_outcome(outcome),
          OutcomeType = OutcomeType,
          Validation = Validation,
          BootstrapR = BootstrapR,
          CVFolds = CVFolds,
          CIBootstrapR = CIBootstrapR,
          CILevel = CILevel,
          Seed = Seed,
          Relabel = Relabel,
          codebook = codebook,
          Verbose = FALSE
        ),
        error = function(e) {
          failure_rows[[length(failure_rows) + 1]] <<- tibble::tibble(
            Outcome = outcome,
            Biomarker = biomarker,
            Status = "failed",
            Note = conditionMessage(e)
          )
          NULL
        }
      )

      if (is.null(evaluation)) next
      evaluation_list[[key]] <- evaluation

      metadata <- evaluation$Metadata
      inferred_positive[[outcome]] <- metadata$PositiveLevel

      perf <- evaluation$PerformanceTable

      if (metadata$OutcomeType == "binary") {
        adjusted_row <- perf %>% dplyr::filter(.data$Model == "Adjusted")

        screen_row <- tibble::tibble(
          Outcome = outcome,
          OutcomeLabel = metadata$OutcomeLabel,
          OutcomeType = metadata$OutcomeType,
          PositiveLevel = metadata$PositiveLevel,
          Biomarker = biomarker,
          BiomarkerLabel = metadata$BiomarkerLabel,
          BiomarkerType = metadata$BiomarkerType,
          N = adjusted_row$N[1],
          NPositive = adjusted_row$NPositive[1],
          NNegative = adjusted_row$NNegative[1],
          NMissing = adjusted_row$NMissing[1],
          Prevalence = adjusted_row$Prevalence[1],
          AUC = perf$AUC[perf$Model == "Biomarker"][1],
          AUC_Lower = perf$AUC_Lower[perf$Model == "Biomarker"][1],
          AUC_Upper = perf$AUC_Upper[perf$Model == "Biomarker"][1],
          CovariateAUC = perf$AUC[perf$Model == "Covariates"][1],
          AdjustedAUC = perf$AUC[perf$Model == "Adjusted"][1],
          AdjustedAUC_Lower = perf$AUC_Lower[perf$Model == "Adjusted"][1],
          AdjustedAUC_Upper = perf$AUC_Upper[perf$Model == "Adjusted"][1],
          DeltaAUC = adjusted_row$DeltaAUC_vs_Covariates[1],
          Brier = adjusted_row$Brier[1],
          Brier_Lower = adjusted_row$Brier_Lower[1],
          Brier_Upper = adjusted_row$Brier_Upper[1],
          CalibrationIntercept = adjusted_row$CalibrationIntercept[1],
          CalibrationSlope = adjusted_row$CalibrationSlope[1],
          ObservedExpectedRatio = adjusted_row$ObservedExpectedRatio[1]
        )
      } else {
        adjusted_row <- perf %>% dplyr::filter(.data$Model == "Adjusted")

        screen_row <- tibble::tibble(
          Outcome = outcome,
          OutcomeLabel = metadata$OutcomeLabel,
          OutcomeType = metadata$OutcomeType,
          PositiveLevel = NA_character_,
          Biomarker = biomarker,
          BiomarkerLabel = metadata$BiomarkerLabel,
          BiomarkerType = metadata$BiomarkerType,
          N = adjusted_row$N[1],
          NPositive = NA_integer_,
          NNegative = NA_integer_,
          NMissing = adjusted_row$NMissing[1],
          Prevalence = NA_real_,
          R2 = perf$R2[perf$Model == "Biomarker"][1],
          R2_Lower = perf$R2_Lower[perf$Model == "Biomarker"][1],
          R2_Upper = perf$R2_Upper[perf$Model == "Biomarker"][1],
          CovariateR2 = perf$R2[perf$Model == "Covariates"][1],
          AdjustedR2 = perf$R2[perf$Model == "Adjusted"][1],
          AdjustedR2_Lower = perf$R2_Lower[perf$Model == "Adjusted"][1],
          AdjustedR2_Upper = perf$R2_Upper[perf$Model == "Adjusted"][1],
          DeltaR2 = adjusted_row$DeltaR2_vs_Covariates[1],
          RMSE = adjusted_row$RMSE[1],
          RMSE_Lower = adjusted_row$RMSE_Lower[1],
          RMSE_Upper = adjusted_row$RMSE_Upper[1],
          MAE = adjusted_row$MAE[1],
          MAE_Lower = adjusted_row$MAE_Lower[1],
          MAE_Upper = adjusted_row$MAE_Upper[1]
        )
      }

      performance_rows[[length(performance_rows) + 1]] <- screen_row

      if (nrow(evaluation$RegressionTable) > 0) {
        regression_rows[[length(regression_rows) + 1]] <- evaluation$RegressionTable %>%
          dplyr::mutate(
            Outcome = outcome,
            OutcomeLabel = metadata$OutcomeLabel,
            .before = 1
          )
      }

      if (nrow(evaluation$ThresholdTable) > 0) {
        threshold_rows[[length(threshold_rows) + 1]] <- evaluation$ThresholdTable %>%
          dplyr::mutate(
            Outcome = outcome,
            OutcomeLabel = metadata$OutcomeLabel,
            Biomarker = biomarker,
            BiomarkerLabel = metadata$BiomarkerLabel,
            .before = 1
          )
      }
    }
  }

  if (length(performance_rows) == 0) {
    stop("No biomarker-outcome analyses completed successfully.", call. = FALSE)
  }

  PerformanceTable <- dplyr::bind_rows(performance_rows)
  RegressionTable <- dplyr::bind_rows(regression_rows)
  ThresholdTable <- dplyr::bind_rows(threshold_rows)
  FailureTable <- dplyr::bind_rows(failure_rows)

  if (nrow(FailureTable) > 0) {
    message(
      nrow(FailureTable),
      " biomarker-outcome analysis", if (nrow(FailureTable) == 1) " was" else "es were",
      " unavailable; see $FailureTable."
    )
  }

  binary_outcomes <- unique(PerformanceTable$Outcome[PerformanceTable$OutcomeType == "binary"])
  if (length(binary_outcomes) > 0) {
    positive_messages <- vapply(binary_outcomes, function(outcome) {
      level <- inferred_positive[[outcome]]
      if (is.null(level)) return(NA_character_)
      paste0(outcome, ': treating "', level, '" as the positive outcome level.')
    }, character(1))
    positive_messages <- unique(stats::na.omit(positive_messages))
    if (length(positive_messages) > 0) message(paste(positive_messages, collapse = "\n"))
  }

  # Build screening plots

  # Keep a stable screening schema so binary-only and continuous-only screens
  # can use the same plotting code without depending on absent columns.
  stable_numeric_columns <- c(
    "AUC", "AUC_Lower", "AUC_Upper", "CovariateAUC", "AdjustedAUC",
    "AdjustedAUC_Lower", "AdjustedAUC_Upper", "DeltaAUC", "Brier",
    "Brier_Lower", "Brier_Upper", "CalibrationIntercept",
    "CalibrationSlope", "ObservedExpectedRatio", "R2", "R2_Lower",
    "R2_Upper", "CovariateR2", "AdjustedR2", "AdjustedR2_Lower",
    "AdjustedR2_Upper", "DeltaR2", "RMSE", "RMSE_Lower", "RMSE_Upper",
    "MAE", "MAE_Lower", "MAE_Upper"
  )

  for (nm in stable_numeric_columns) {
    if (!nm %in% names(PerformanceTable)) PerformanceTable[[nm]] <- NA_real_
  }

  if (identical(HeatmapMetric, "AdjustedAUC") &&
      all(PerformanceTable$OutcomeType == "continuous")) {
    HeatmapMetric <- "AdjustedR2"
    message("All screened outcomes are continuous; using AdjustedR2 for the heatmap.")
  }

  Heatmap <- NULL
  AUCPlot <- NULL
  IncrementalAUCPlot <- NULL
  R2Plot <- NULL
  IncrementalR2Plot <- NULL
  BiomarkerPanels <- NULL
  ROCFacets <- NULL

  metric_column <- switch(
    HeatmapMetric,
    AdjustedAUC = "AdjustedAUC",
    AUC = "AUC",
    DeltaAUC = "DeltaAUC",
    AdjustedR2 = "AdjustedR2",
    R2 = "R2",
    DeltaR2 = "DeltaR2",
    HeatmapMetric
  )

  if (!metric_column %in% names(PerformanceTable)) {
    stop(
      "HeatmapMetric was not found in the screening results. Requested: ",
      HeatmapMetric,
      ". Available numeric metrics include: ",
      paste(names(PerformanceTable)[vapply(PerformanceTable, is.numeric, logical(1))], collapse = ", "),
      call. = FALSE
    )
  }

  heatmap_data <- PerformanceTable %>%
    dplyr::mutate(
      HeatmapValue = .data[[metric_column]],
      Tooltip = dplyr::case_when(
        .data$OutcomeType == "binary" ~ paste0(
          .data$BiomarkerLabel, "<br>",
          .data$OutcomeLabel, "<br><br>",
          "Adjusted AUC: ", sprintf("%.3f", .data$AdjustedAUC), "<br>",
          "95% CI: ", sprintf("%.3f", .data$AdjustedAUC_Lower), " to ", sprintf("%.3f", .data$AdjustedAUC_Upper), "<br>",
          "Biomarker AUC: ", sprintf("%.3f", .data$AUC), "<br>",
          "Covariate AUC: ", sprintf("%.3f", .data$CovariateAUC), "<br>",
          "Delta AUC: ", sprintf("%+.3f", .data$DeltaAUC), "<br>",
          "Brier: ", sprintf("%.3f", .data$Brier), "<br>",
          "N: ", .data$N
        ),
        TRUE ~ paste0(
          .data$BiomarkerLabel, "<br>",
          .data$OutcomeLabel, "<br><br>",
          "Adjusted R2: ", sprintf("%.3f", .data$AdjustedR2), "<br>",
          "Biomarker R2: ", sprintf("%.3f", .data$R2), "<br>",
          "Covariate R2: ", sprintf("%.3f", .data$CovariateR2), "<br>",
          "Delta R2: ", sprintf("%+.3f", .data$DeltaR2), "<br>",
          "RMSE: ", sprintf("%.3f", .data$RMSE), "<br>",
          "MAE: ", sprintf("%.3f", .data$MAE), "<br>",
          "N: ", .data$N
        )
      )
    )

  Heatmap <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(
      x = .data$OutcomeLabel,
      y = .data$BiomarkerLabel,
      fill = .data$HeatmapValue,
      text = .data$Tooltip
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(na.value = "grey90") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = HeatmapMetric,
      title = "Biomarker performance heatmap"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  binary_data <- PerformanceTable %>% dplyr::filter(.data$OutcomeType == "binary")
  if (nrow(binary_data) > 0) {
    AUCPlot <- ggplot2::ggplot(
      binary_data,
      ggplot2::aes(
        x = .data$AdjustedAUC,
        y = stats::reorder(.data$BiomarkerLabel, .data$AdjustedAUC),
        xmin = .data$AdjustedAUC_Lower,
        xmax = .data$AdjustedAUC_Upper
      )
    ) +
      ggplot2::geom_vline(xintercept = 0.5, linetype = 2) +
      ggplot2::geom_errorbar(height = 0.15, orientation = "y") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~OutcomeLabel, scales = "free_y") +
      ggplot2::labs(x = "Adjusted AUC", y = NULL, title = "Adjusted biomarker discrimination") +
      ggplot2::theme_minimal()

    IncrementalAUCPlot <- ggplot2::ggplot(
      binary_data,
      ggplot2::aes(
        x = .data$DeltaAUC,
        y = stats::reorder(.data$BiomarkerLabel, .data$DeltaAUC)
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$DeltaAUC, yend = .data$BiomarkerLabel),
        linewidth = 0.6
      ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~OutcomeLabel, scales = "free_y") +
      ggplot2::labs(
        x = "Delta AUC after adding biomarker to covariates",
        y = NULL,
        title = "Incremental biomarker discrimination"
      ) +
      ggplot2::theme_minimal()
  }

  continuous_data <- PerformanceTable %>% dplyr::filter(.data$OutcomeType == "continuous")
  if (nrow(continuous_data) > 0) {
    R2Plot <- ggplot2::ggplot(
      continuous_data,
      ggplot2::aes(
        x = .data$AdjustedR2,
        y = stats::reorder(.data$BiomarkerLabel, .data$AdjustedR2),
        xmin = .data$AdjustedR2_Lower,
        xmax = .data$AdjustedR2_Upper
      )
    ) +
      ggplot2::geom_errorbar(height = 0.15, orientation = "y") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~OutcomeLabel, scales = "free_y") +
      ggplot2::labs(x = "Adjusted model R2", y = NULL, title = "Adjusted biomarker model fit") +
      ggplot2::theme_minimal()

    IncrementalR2Plot <- ggplot2::ggplot(
      continuous_data,
      ggplot2::aes(
        x = .data$DeltaR2,
        y = stats::reorder(.data$BiomarkerLabel, .data$DeltaR2)
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$DeltaR2, yend = .data$BiomarkerLabel),
        linewidth = 0.6
      ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~OutcomeLabel, scales = "free_y") +
      ggplot2::labs(
        x = "Delta R2 after adding biomarker to covariates",
        y = NULL,
        title = "Incremental biomarker explanatory value"
      ) +
      ggplot2::theme_minimal()
  }

  # Visual QC for continuous biomarkers screened against binary outcomes.
  # These panels deliberately use each evaluation's own complete-case cohort.
  panel_rows <- list()
  roc_facet_rows <- list()

  for (key in names(evaluation_list)) {
    evaluation <- evaluation_list[[key]]
    metadata <- evaluation$Metadata

    if (metadata$OutcomeType != "binary" ||
        metadata$BiomarkerType != "continuous") {
      next
    }

    adjusted_row <- evaluation$PerformanceTable %>%
      dplyr::filter(.data$Model == "Adjusted")

    if (nrow(adjusted_row) != 1 || nrow(evaluation$Predictions) == 0) {
      next
    }

    panel_rows[[length(panel_rows) + 1]] <- evaluation$Predictions %>%
      dplyr::transmute(
        OutcomeLabel = metadata$OutcomeLabel,
        BiomarkerLabel = metadata$BiomarkerLabel,
        OutcomeLevel = .data$ObservedLevel,
        BiomarkerValue = .data$Biomarker,
        AUC = adjusted_row$AUC,
        AUC_Lower = adjusted_row$AUC_Lower,
        AUC_Upper = adjusted_row$AUC_Upper,
        N = adjusted_row$N
      )

    if (nrow(evaluation$ROCData) > 0) {
      roc_facet_rows[[length(roc_facet_rows) + 1]] <- evaluation$ROCData %>%
        dplyr::filter(.data$Model == "Adjusted") %>%
        dplyr::mutate(
          OutcomeLabel = metadata$OutcomeLabel,
          BiomarkerLabel = metadata$BiomarkerLabel,
          AUC = adjusted_row$AUC,
          AUC_Lower = adjusted_row$AUC_Lower,
          AUC_Upper = adjusted_row$AUC_Upper,
          N = adjusted_row$N
        )
    }
  }

  if (length(panel_rows) > 0) {
    panel_data <- dplyr::bind_rows(panel_rows)
    panel_labels <- panel_data %>%
      dplyr::distinct(
        .data$OutcomeLabel, .data$BiomarkerLabel, .data$AUC,
        .data$AUC_Lower, .data$AUC_Upper, .data$N
      ) %>%
      dplyr::mutate(
        MetricLabel = paste0(
          "Adjusted AUC = ", sprintf("%.2f", .data$AUC),
          " (", sprintf("%.2f", .data$AUC_Lower),
          "–", sprintf("%.2f", .data$AUC_Upper), ")\n",
          "N = ", .data$N
        )
      )

    BiomarkerPanels <- ggplot2::ggplot(
      panel_data,
      ggplot2::aes(x = .data$OutcomeLevel, y = .data$BiomarkerValue)
    ) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.25) +
      ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
      ggplot2::geom_text(
        data = panel_labels,
        ggplot2::aes(x = Inf, y = Inf, label = .data$MetricLabel),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.15,
        size = 3
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$OutcomeLabel),
        cols = ggplot2::vars(.data$BiomarkerLabel),
        scales = "free_y"
      ) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        title = "Continuous biomarker distributions by outcome"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }

  if (length(roc_facet_rows) > 0) {
    roc_facet_data <- dplyr::bind_rows(roc_facet_rows)
    roc_labels <- roc_facet_data %>%
      dplyr::distinct(
        .data$OutcomeLabel, .data$BiomarkerLabel, .data$AUC,
        .data$AUC_Lower, .data$AUC_Upper, .data$N
      ) %>%
      dplyr::mutate(
        FalsePositiveRate = 0.98,
        Sensitivity = 0.02,
        MetricLabel = paste0(
          "AUC = ", sprintf("%.2f", .data$AUC),
          " (", sprintf("%.2f", .data$AUC_Lower),
          "–", sprintf("%.2f", .data$AUC_Upper), ")\n",
          "N = ", .data$N
        )
      )

    ROCFacets <- ggplot2::ggplot(
      roc_facet_data,
      ggplot2::aes(x = .data$FalsePositiveRate, y = .data$Sensitivity)
    ) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_text(
        data = roc_labels,
        ggplot2::aes(
          x = .data$FalsePositiveRate,
          y = .data$Sensitivity,
          label = .data$MetricLabel
        ),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 0,
        size = 3
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$OutcomeLabel),
        cols = ggplot2::vars(.data$BiomarkerLabel)
      ) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        title = "Adjusted ROC curves by biomarker"
      ) +
      ggplot2::theme_minimal()
  }

  list(
    PerformanceTable = PerformanceTable,
    RegressionTable = RegressionTable,
    ThresholdTable = ThresholdTable,
    FailureTable = FailureTable,
    Evaluations = evaluation_list,
    Plots = list(
      Heatmap = Heatmap,
      AUC = AUCPlot,
      IncrementalAUC = IncrementalAUCPlot,
      R2 = R2Plot,
      IncrementalR2 = IncrementalR2Plot,
      BiomarkerPanels = BiomarkerPanels,
      ROCFacets = ROCFacets
    ),
    Metadata = list(
      OutcomeVariables = outcome_vars,
      BiomarkerVariables = biomarker_vars,
      Covariates = covariates,
      ValidationMethod = Validation,
      HeatmapMetric = HeatmapMetric,
      CILevel = CILevel,
      Relabeled = Relabel,
      PairwiseCompleteCases = TRUE
    )
  )
}
