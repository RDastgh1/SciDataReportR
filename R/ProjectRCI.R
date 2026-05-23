
#' Project a trained RCI object onto new data
#'
#' Apply previously learned regression-based Reliable Change Index (RCI)
#' models to a new or expanded longitudinal dataset WITHOUT refitting models.
#'
#' Supports both wide and long data structures and returns:
#' - tidy longitudinal RCI values
#' - original data merged with projected RCI outputs
#' - publication-ready plots
#'
#' @param Data A data frame.
#' @param Object A SciDataReportR_RCI object created with
#'   `CreateRCIObject()`.
#' @param ID Optional ID column override.
#'
#' @return A projected SciDataReportR_RCI object.
#'
#' @export
ProjectRCI <- function(
    Data,
    Object,
    ID = NULL
) {

  # ============================================================================
  # Extract metadata
  # ============================================================================

  Metadata <- Object$Metadata

  Variables <- Metadata$VariablesUsed

  DataFormat <- Metadata$DataFormat

  BaselineVisit <- Metadata$BaselineVisit

  z_threshold <- Metadata$ZThreshold

  if (is.null(ID)) {

    ID <- Metadata$IDColumn

  }

  if (!ID %in% names(Data)) {

    stop(
      paste0(
        "ID column not found: ",
        ID
      )
    )

  }

  ClassificationColors <- c(
    "Reliable Improvement" = "#1f4e79",
    "Stable" = "grey70",
    "Reliable Decline" = "#b35806"
  )

  # ============================================================================
  # Prepare longitudinal data
  # ============================================================================

  if (DataFormat == "long") {

    VisitColumn <- Metadata$VisitColumn

    VisitOrder <- Metadata$VisitOrder

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
        Visit = tidyselect::all_of(
          VisitColumn
        )
      ) %>%
      dplyr::mutate(
  ID_Internal = .data[[ID]],
  !!ID := .data[[ID]]
)

  }

  if (DataFormat == "wide") {

    BaselineSpecifier <- Metadata$BaselineSpecifier

    FollowupSpecifier <- Metadata$FollowupSpecifier

    SpecifierPosition <- Metadata$SpecifierPosition

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
  UserID = Data[[ID]],
  Variable = var,
  Visit = BaselineVisit,
  Value = Data[[baseline_var]]
)

names(baseline_df)[
  names(baseline_df) == "UserID"
] <- ID

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

    out <- tibble::tibble(
      ID_Internal = Data[[ID]],
      UserID = Data[[ID]],
      Variable = var,
      Visit = VisitName,
      Value = Data[[fu_var]]
    )

    names(out)[
      names(out) == "UserID"
    ] <- ID

    out

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
  # Project models WITHOUT refitting
  # ============================================================================

  ProjectionList <- list()

  for (var in Variables) {

    if (!var %in% names(Object$Models)) {
      next
    }

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

    if (nrow(merged) == 0) {
      next
    }

    model <- Object$Models[[var]]$Model

    see <- Object$Models[[var]]$SEE

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

  }

  if (length(ProjectionList) == 0) {

    stop(
      "No valid RCI projections could be generated."
    )

  }

  # ============================================================================
  # Build RCIValues
  # ============================================================================

  RCIValues <- dplyr::bind_rows(
    ProjectionList
  )

  RCIValues$UserID <- RCIValues$ID_Internal

  names(RCIValues)[
    names(RCIValues) == "UserID"
  ] <- ID

  RCIValues <- RCIValues %>%
    dplyr::left_join(
      Object$Variables,
      by = "Variable"
    ) %>%
    dplyr::select(
      -ID_Internal
    ) %>%
    dplyr::relocate(
      tidyselect::all_of(ID),
      Visit,
      Variable,
      Label,
      Baseline,
      Followup,
      Predicted,
      Residual,
      RCI,
      Classification
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

    # ============================================================================
  # Plots
  # ============================================================================

  SpaghettiPlots <- list()
  WaterfallPlots <- list()
  QuadrantPlots <- list()

  for (var in Variables) {

    plot_label <- Object$Variables$Label[
      Object$Variables$Variable == var
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

    if (!is.null(Metadata$VisitOrder)) {

      spaghetti_data$Visit <- factor(
        spaghetti_data$Visit,
        levels = Metadata$VisitOrder,
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
        SegmentClass = Classification
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
  # Return
  # ============================================================================

  Output <- list(

    Metadata = list(
      ProjectionN = nrow(Data),
      ProjectionTimestamp = Sys.time(),
      SourceObjectTimestamp = Metadata$Timestamp
    ),

    Variables = Object$Variables,

    Models = Object$Models,

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
