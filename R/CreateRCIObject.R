
#' Create a Reliable Change Index (RCI) object
#'
#' Learn regression-based Reliable Change Index (RCI) models relative to a
#' user-defined reference visit and calculate projected RCI values.
#'
#' Supports both wide and long longitudinal data structures.
#'
#' Long format is recommended for datasets with more than two visits.
#'
#' @param Data A data frame.
#' @param Variables Character vector of canonical variable names.
#' @param DataFormat Either "wide" or "long".
#' @param ID ID column.
#' @param Method Currently only "regression" is supported.
#'
#' @param BaselineSpecifier Baseline visit identifier for wide data.
#' @param FollowupSpecifier Follow-up visit identifier for wide data.
#' @param SpecifierPosition Either "suffix" or "prefix".
#'
#' @param VisitColumn Visit column for long data.
#' @param VisitOrder Optional ordering of visits.
#' @param BaselineVisit Reference visit used for RCI calculations.
#'
#' @param Confidence Confidence interval threshold.
#' @param Relabel Logical; use variable labels when available.
#'
#'
#'#'

#' ## Interpretation guide

#'

#' egression-based RCI values are interpreted similarly to z-scores.

#'
#' | RCI cutoff | Approximate confidence interval |
#' |:-----------|:--------------------------------|
#' | +/-0.50 | ~38% |
#' | +/-1.00 | ~68% |
#' | +/-1.645 | ~90% |
#' | +/-1.96 | ~95% |
#' | +/-2.58 | ~99% |
#'

#' Traditional Jacobson-Truax RCI thresholds typically use +/-1.96,
#' corresponding to approximately 95% confidence.

#'
#' @return A SciDataReportR_RCI object.
#'
#' @export
CreateRCIObject <- function(
    Data,
    Variables,
    DataFormat = c("wide", "long"),
    ID,
    Method = "regression",

    BaselineSpecifier = NULL,
    FollowupSpecifier = NULL,
    SpecifierPosition = c("suffix", "prefix"),

    VisitColumn = NULL,
    VisitOrder = NULL,
    BaselineVisit = NULL,

    Confidence = 0.95,
    Relabel = TRUE
) {

  DataFormat <- match.arg(DataFormat)
  SpecifierPosition <- match.arg(SpecifierPosition)

  z_threshold <- stats::qnorm(
    1 - ((1 - Confidence) / 2)
  )

  ClassificationColors <- c(
    "Reliable Improvement" = "#1f4e79",
    "Stable" = "grey70",
    "Reliable Decline" = "#b35806"
  )

  # ============================================================================
  # Prepare longitudinal data
  # ============================================================================

  if (DataFormat == "long") {

    if (is.null(VisitColumn)) {
      stop("VisitColumn must be supplied for long data.")
    }

    if (is.null(BaselineVisit)) {
      stop("BaselineVisit must be supplied for long data.")
    }

    if (!is.null(VisitOrder)) {

      Data[[VisitColumn]] <- factor(
        Data[[VisitColumn]],
        levels = VisitOrder,
        ordered = TRUE
      )

    }

    LongData <- tidyr::pivot_longer(
      Data,
      cols = tidyselect::all_of(Variables),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
      dplyr::rename(
        Visit = tidyselect::all_of(VisitColumn)
      ) %>%
      dplyr::mutate(
        ID_Internal = .data[[ID]]
      )

  }

  if (DataFormat == "wide") {

    long_list <- list()

    for (var in Variables) {

      if (SpecifierPosition == "prefix") {

        baseline_var <- paste0(
          BaselineSpecifier,
          var
        )

        followup_pattern <- paste0(
          "^",
          FollowupSpecifier,
          var
        )

      } else {

        baseline_var <- paste0(
          var,
          BaselineSpecifier
        )

        followup_pattern <- paste0(
          "^",
          var,
          FollowupSpecifier
        )

      }

      if (!baseline_var %in% names(Data)) {
        next
      }

      followup_vars <- names(Data)[
        grepl(
          followup_pattern,
          names(Data)
        )
      ]

      if (length(followup_vars) == 0) {
        next
      }

      baseline_df <- tibble::tibble(
        ID_Internal = Data[[ID]],
        Variable = var,
        Visit = BaselineVisit,
        Value = Data[[baseline_var]]
      )

      followup_df <- purrr::map_dfr(
        followup_vars,
        function(fu_var) {

          if (SpecifierPosition == "prefix") {

            VisitName <- stringr::str_remove(
              fu_var,
              paste0("_?", var, "$")
            )

          } else {

            VisitName <- stringr::str_remove(
              fu_var,
              paste0("^", var, "_?")
            )

          }

          tibble::tibble(
            ID_Internal = Data[[ID]],
            Variable = var,
            Visit = VisitName,
            Value = Data[[fu_var]]
          )

        }
      )

      long_list[[var]] <- dplyr::bind_rows(
        baseline_df,
        followup_df
      )

    }

    LongData <- dplyr::bind_rows(
      long_list
    )

  }

  # ============================================================================
  # Variable labels
  # ============================================================================

  VariableTable <- tibble::tibble(
    Variable = Variables
  )

  if (Relabel) {

    VariableTable$Label <- purrr::map_chr(
      Variables,
      function(v) {

        label_source <- NULL

        if (DataFormat == "long") {

          label_source <- Data[[v]]

        }

        if (DataFormat == "wide") {

          if (SpecifierPosition == "prefix") {

            label_source <- Data[[paste0(
              BaselineSpecifier,
              v
            )]]

          } else {

            label_source <- Data[[paste0(
              v,
              BaselineSpecifier
            )]]

          }

        }

        label <- sjlabelled::get_label(
          label_source
        )

        if (
          length(label) == 0 ||
          is.null(label) ||
          is.na(label)
        ) {

          return(v)

        }

        as.character(label)

      }
    )

  } else {

    VariableTable$Label <- Variables

  }

  # ============================================================================
  # Fit models
  # ============================================================================

  Models <- list()
  ProjectionList <- list()

  for (var in Variables) {

    temp <- LongData %>%
      dplyr::filter(
        Variable == var
      )

    baseline_values <- temp %>%
      dplyr::filter(
        Visit == BaselineVisit
      ) %>%
      dplyr::select(
        ID_Internal,
        Baseline = Value
      )

    followup_values <- temp %>%
      dplyr::filter(
        Visit != BaselineVisit
      ) %>%
      dplyr::select(
        ID_Internal,
        Visit,
        Followup = Value
      )

    merged <- dplyr::left_join(
      followup_values,
      baseline_values,
      by = "ID_Internal"
    ) %>%
      dplyr::filter(
        is.finite(Baseline),
        is.finite(Followup)
      )

    if (nrow(merged) < 5) {
      next
    }

    model <- stats::lm(
      Followup ~ Baseline,
      data = merged
    )

    see <- summary(model)$sigma

    merged$Predicted <- stats::predict(
      model,
      newdata = merged
    )

    merged$Residual <- (
      merged$Followup -
      merged$Predicted
    )

    merged$RCI <- (
      merged$Residual / see
    )

    merged$Classification <- dplyr::case_when(
      merged$RCI >= z_threshold ~ "Reliable Improvement",
      merged$RCI <= -z_threshold ~ "Reliable Decline",
      TRUE ~ "Stable"
    )

    merged$Classification <- factor(
      merged$Classification,
      levels = c(
        "Reliable Decline",
        "Stable",
        "Reliable Improvement"
      )
    )

    merged$Variable <- var

    ProjectionList[[var]] <- merged

    Models[[var]] <- list(
      Model = model,
      SEE = see,
      Threshold = z_threshold,
      N = nrow(merged),
      RSquared = summary(model)$r.squared
    )

  }

  if (length(ProjectionList) == 0) {

    stop(
      "No valid RCI models could be fit."
    )

  }

  RCIValues <- dplyr::bind_rows(
    ProjectionList
  )

  RCIValues[[ID]] <- RCIValues$ID_Internal

  RCIValues <- RCIValues %>%
    dplyr::left_join(
      VariableTable,
      by = "Variable"
    ) %>%
    dplyr::select(
      tidyselect::all_of(ID),
      Visit,
      Variable,
      Label,
      Baseline,
      Followup,
      Predicted,
      Residual,
      RCI,
      Classification,
      dplyr::everything(),
      -ID_Internal
    )

  # ============================================================================
  # CombinedData
  # ============================================================================

  if (DataFormat == "wide") {

    RCIWide <- RCIValues %>%
      dplyr::mutate(

        VisitClean = make.names(
          Visit
        ),

        RCI_Name = paste0(
          Variable,
          "_RCI_",
          VisitClean
        ),

        Class_Name = paste0(
          Variable,
          "_ChangeClassification_",
          VisitClean
        )

      )

    print(names(RCIWide))
    RCIOnly <- RCIWide %>%
      dplyr::select(
        tidyselect::all_of(ID),
        RCI_Name,
        RCI
      ) %>%
      tidyr::pivot_wider(
        names_from = RCI_Name,
        values_from = RCI
      )

    ClassOnly <- RCIWide %>%
      dplyr::select(
        tidyselect::all_of(ID),
        Class_Name,
        Classification
      ) %>%
      tidyr::pivot_wider(
        names_from = Class_Name,
        values_from = Classification
      )

    CombinedData <- Data %>%
      dplyr::left_join(
        RCIOnly,
        by = ID
      ) %>%
      dplyr::left_join(
        ClassOnly,
        by = ID
      )

  }

  if (DataFormat == "long") {

    DataJoin <- Data %>%
      dplyr::rename(
        Visit = tidyselect::all_of(
          VisitColumn
        )
      )

    RCIWide <- RCIValues %>%
      dplyr::mutate(
        RCI_Name = paste0(
          Variable,
          "_RCI"
        ),
        Class_Name = paste0(
          Variable,
          "_ChangeClassification"
        )
      )

    RCIOnly <- RCIWide %>%
      dplyr::select(
        tidyselect::all_of(ID),
        Visit,
        RCI_Name,
        RCI
      ) %>%
      tidyr::pivot_wider(
        names_from = RCI_Name,
        values_from = RCI
      )

    ClassOnly <- RCIWide %>%
      dplyr::select(
        tidyselect::all_of(ID),
        Visit,
        Class_Name,
        Classification
      ) %>%
      tidyr::pivot_wider(
        names_from = Class_Name,
        values_from = Classification
      )

    CombinedData <- DataJoin %>%
      dplyr::left_join(
        RCIOnly,
        by = c(
          ID,
          "Visit"
        )
      ) %>%
      dplyr::left_join(
        ClassOnly,
        by = c(
          ID,
          "Visit"
        )
      ) %>%
      dplyr::rename(
        !!VisitColumn := Visit
      )

  }

  # Add labels back

    # ============================================================================
  # Apply labels to generated columns
  # ============================================================================

  for (i in seq_len(nrow(VariableTable))) {

    var <- VariableTable$Variable[i]

    label <- VariableTable$Label[i]

    rci_cols <- names(CombinedData)[
      grepl(
        paste0("^", var, "_RCI"),
        names(CombinedData)
      )
    ]

    class_cols <- names(CombinedData)[
      grepl(
        paste0("^", var, "_ChangeClassification"),
        names(CombinedData)
      )
    ]

    # --------------------------------------------------------------------------
    # RCI labels
    # --------------------------------------------------------------------------

    for (col in rci_cols) {

      attr(
        CombinedData[[col]],
        "label"
      ) <- paste0(
        label,
        " RCI"
      )

    }

    # --------------------------------------------------------------------------
    # Classification labels
    # --------------------------------------------------------------------------

    for (col in class_cols) {

      attr(
        CombinedData[[col]],
        "label"
      ) <- paste0(
        label,
        " Change Classification"
      )

    }

  }

    # ============================================================================
  # Plots
  # ============================================================================

  SpaghettiPlots <- list()
  WaterfallPlots <- list()
  QuadrantPlots <- list()

  for (var in Variables) {

    plot_label <- VariableTable$Label[
      VariableTable$Variable == var
    ]

    temp_long <- LongData %>%
      dplyr::filter(
        Variable == var
      )

    temp_proj <- RCIValues %>%
      dplyr::filter(
        Variable == var
      )

    spaghetti_data <- temp_long %>%
      dplyr::left_join(
        temp_proj %>%
          dplyr::select(
            tidyselect::all_of(ID),
            Variable,
            Visit,
            Classification
          ),
        by = c(
          ID,
          "Variable",
          "Visit"
        )
      )

    spaghetti_data$Classification[
      is.na(spaghetti_data$Classification)
    ] <- "Stable"

    spaghetti_data$Classification <- factor(
      spaghetti_data$Classification,
      levels = c(
        "Reliable Decline",
        "Stable",
        "Reliable Improvement"
      )
    )

    if (!is.null(VisitOrder)) {

      spaghetti_data$Visit <- factor(
        spaghetti_data$Visit,
        levels = VisitOrder,
        ordered = TRUE
      )

    }

   segment_data <- spaghetti_data %>%
  dplyr::arrange(
    .data[[ID]],
    Visit
  ) %>%
  dplyr::group_by(
    .data[[ID]]
  ) %>%
  dplyr::mutate(
    xend = dplyr::lead(Visit),
    yend = dplyr::lead(Value),
    SegmentClass = dplyr::lead(Classification)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(xend),
    !is.na(yend)
  )

    temp_proj$PlotOrder <- rank(
      temp_proj$RCI,
      ties.method = "first"
    )

    # --------------------------------------------------------------------------
    # Spaghetti plot
    # --------------------------------------------------------------------------

    SpaghettiPlots[[var]] <- ggplot2::ggplot() +

      ggplot2::geom_segment(
        data = segment_data,
        ggplot2::aes(
          x = Visit,
          y = Value,
          xend = xend,
          yend = yend,
          color = SegmentClass,
          group = .data[[ID]]
        ),
        linewidth = 1,
        alpha = 0.8
      ) +

      ggplot2::geom_point(
        data = spaghetti_data,
        ggplot2::aes(
          x = Visit,
          y = Value
        ),
        color = "grey20",
        size = 2
      ) +

      ggplot2::scale_color_manual(
        values = ClassificationColors,
        drop = FALSE
      ) +

      ggplot2::labs(
        title = plot_label,
        subtitle = paste0(
          "Reference visit: ",
          BaselineVisit
        ),
        x = NULL,
        y = "Value",
        color = "Classification"
      ) +

      ggplot2::theme_minimal()

    # --------------------------------------------------------------------------
    # Waterfall plot
    # --------------------------------------------------------------------------

    WaterfallPlots[[var]] <- ggplot2::ggplot(
      temp_proj,
      ggplot2::aes(
        x = reorder(
          PlotOrder,
          RCI
        ),
        y = RCI,
        fill = Classification
      )
    ) +

      ggplot2::geom_col(
        width = 0.9
      ) +

      ggplot2::geom_hline(
        yintercept = c(
          -z_threshold,
          z_threshold
        ),
        linetype = "dashed",
        color = "grey40"
      ) +

      ggplot2::scale_fill_manual(
        values = ClassificationColors,
        drop = FALSE
      ) +

      ggplot2::labs(
        title = plot_label,
        subtitle = paste0(
          "Reference visit: ",
          BaselineVisit
        ),
        x = NULL,
        y = "RCI",
        fill = "Classification"
      ) +

      ggplot2::theme_minimal() +

      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )

    # --------------------------------------------------------------------------
    # Quadrant plot
    # --------------------------------------------------------------------------

    QuadrantPlots[[var]] <- ggplot2::ggplot(
      temp_proj,
      ggplot2::aes(
        x = Baseline,
        y = RCI,
        color = Classification
      )
    ) +

      ggplot2::geom_hline(
        yintercept = c(
          -z_threshold,
          z_threshold
        ),
        linetype = "dashed",
        color = "grey40"
      ) +

      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "grey50"
      ) +

      ggplot2::geom_point(
        size = 2.5,
        alpha = 0.85
      ) +

      ggplot2::scale_color_manual(
        values = ClassificationColors,
        drop = FALSE
      ) +

      ggplot2::labs(
        title = plot_label,
        subtitle = paste0(
          "Reference visit: ",
          BaselineVisit
        )
      ) +

      ggplot2::theme_minimal()

  }

    # ============================================================================
  # Generated variable metadata
  # ============================================================================

  GeneratedVariables <- tibble::tibble()

  for (i in seq_len(nrow(VariableTable))) {

    var <- VariableTable$Variable[i]

    label <- VariableTable$Label[i]

    rci_cols <- names(CombinedData)[
      grepl(
        paste0("^", var, "_RCI"),
        names(CombinedData)
      )
    ]

    class_cols <- names(CombinedData)[
      grepl(
        paste0("^", var, "_ChangeClassification"),
        names(CombinedData)
      )
    ]

    # --------------------------------------------------------------------------
    # RCI variables
    # --------------------------------------------------------------------------

    if (length(rci_cols) > 0) {

      GeneratedVariables <- dplyr::bind_rows(
        GeneratedVariables,
        tibble::tibble(
          Variable = rci_cols,
          Label = paste0(
            label,
            " RCI"
          ),
          Type = "RCI"
        )
      )

    }

    # --------------------------------------------------------------------------
    # Classification variables
    # --------------------------------------------------------------------------

    if (length(class_cols) > 0) {

      GeneratedVariables <- dplyr::bind_rows(
        GeneratedVariables,
        tibble::tibble(
          Variable = class_cols,
          Label = paste0(
            label,
            " Change Classification"
          ),
          Type = "Classification"
        )
      )

    }

  }

  # ============================================================================
  # Return
  # ============================================================================

  Output <- list(

    Metadata = list(
      Method = Method,
      Confidence = Confidence,
      ZThreshold = z_threshold,
      DataFormat = DataFormat,
      BaselineVisit = BaselineVisit,
      VariablesUsed = Variables,
      N = nrow(Data),
      IDColumn = ID,

      VisitColumn = VisitColumn,
      VisitOrder = VisitOrder,

      BaselineSpecifier = BaselineSpecifier,
      FollowupSpecifier = FollowupSpecifier,
      SpecifierPosition = SpecifierPosition,

      Timestamp = Sys.time()
    ),

    Variables = VariableTable,

GeneratedVariables = GeneratedVariables,

    Models = Models,

    RCIValues = RCIValues,

    CombinedData = CombinedData,
    Plots = list(
  Spaghetti = SpaghettiPlots,
  Waterfall = WaterfallPlots,
  Quadrant = QuadrantPlots
)

  )

  class(Output) <- "SciDataReportR_RCI"

  return(Output)

}
