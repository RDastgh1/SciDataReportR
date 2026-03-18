#' Plot a spider chart across continuous and binary variables
#'
#' This function summarizes a set of variables and displays them on a spider
#' chart. Continuous variables are plotted as mean z-scores by default using
#' `CalcZScore()`, while binary variables are plotted as percentages. It can
#' overlay groups on one spider chart or facet by group, optionally fill the
#' polygons, relabel spokes using variable labels, wrap long labels, reorder
#' variables to visually emphasize between-group differences, and optionally
#' return an interactive radar chart using plotly.
#'
#' @param Data A data frame.
#' @param Variables Character vector of variable names to plot.
#' @param GroupVariable Optional grouping variable name. If NULL, one overall
#'   summary profile is plotted.
#' @param Relabel Logical; if TRUE, use variable labels when available.
#' @param ContinuousSummary Character; one of `"mean"` or `"median"`.
#' @param ContinuousScaling Character; one of `"zscore"`, `"none"`, or `"minmax"`.
#' @param Fill Logical; if TRUE, add transparent polygon fills in the static
#'   ggplot version and filled polygons in the interactive version.
#' @param FillAlpha Numeric transparency for fills.
#' @param Facet Logical; if TRUE and `GroupVariable` is supplied, facet by group
#'   instead of overlaying all groups on one spider chart. Ignored when
#'   `MakeInteractive = TRUE`.
#' @param VariableOrder Character; one of `"input"`, `"discrimination"`,
#'   `"hierarchical"`, `"greedy"`, or `"category_discrimination"`.
#' @param VariableCategories Optional character vector of categories for
#'   `Variables`. Must be the same length as `Variables` when supplied.
#' @param BinaryPositiveValue Optional positive value to use for non-factor
#'   binary variables. Defaults to `1`. For factor variables, the second factor
#'   level is used.
#' @param Palette Character name of an `hcl.colors()` palette. Defaults to
#'   `"Dark 3"`.
#' @param LineSize Numeric line width for the static ggplot version.
#' @param PointSize Numeric point size for the static ggplot version.
#' @param ShowPoints Logical; if TRUE, show points at each spoke in the static
#'   ggplot version.
#' @param LegendTitle Optional legend title. Defaults to `GroupVariable`.
#' @param PlotTitle Optional plot title.
#' @param Subtitle Optional plot subtitle.
#' @param Caption Optional plot caption.
#' @param AxisLabelSize Numeric axis text size for spoke labels in the static
#'   ggplot version.
#' @param AxisTextSize Numeric text size for radial axis labels in the static
#'   ggplot version.
#' @param StripTextSize Numeric facet strip text size in the static ggplot
#'   version.
#' @param WrapLabels Logical; if TRUE, wrap long spoke labels.
#' @param LabelWrapWidth Numeric wrap width passed to `stringr::str_wrap()`.
#' @param LabelRadiusMultiplier Numeric multiplier controlling how far labels
#'   sit outside the spider in the static ggplot version.
#' @param PlotMarginTop Numeric top plot margin for the static ggplot version.
#' @param PlotMarginRight Numeric right plot margin for the static ggplot version.
#' @param PlotMarginBottom Numeric bottom plot margin for the static ggplot version.
#' @param PlotMarginLeft Numeric left plot margin for the static ggplot version.
#' @param MakeInteractive Logical; if TRUE, return an interactive plotly radar
#'   chart instead of a static ggplot.
#' @param InteractiveHeight Numeric height in pixels for the interactive widget.
#' @param InteractiveWidth Optional width passed to plotly layout. Defaults to
#'   NULL.
#' @param InteractiveAxisMin Optional numeric minimum for the interactive radial
#'   axis. If NULL, auto-detected from the summarized values.
#' @param InteractiveAxisMax Optional numeric maximum for the interactive radial
#'   axis. If NULL, auto-detected from the summarized values.
#' @param TooltipDigits Integer number of digits to show in interactive tooltips.
#'
#' @return A ggplot object when `MakeInteractive = FALSE`, otherwise a plotly
#'   htmlwidget.
#' @export
PlotSpiderChart <- function(Data,
                            Variables,
                            GroupVariable = NULL,
                            Relabel = TRUE,
                            ContinuousSummary = "mean",
                            ContinuousScaling = "zscore",
                            Fill = FALSE,
                            FillAlpha = 0.2,
                            Facet = FALSE,
                            VariableOrder = "input",
                            VariableCategories = NULL,
                            BinaryPositiveValue = 1,
                            Palette = "Dark 3",
                            LineSize = 1,
                            PointSize = 2,
                            ShowPoints = FALSE,
                            LegendTitle = NULL,
                            PlotTitle = NULL,
                            Subtitle = NULL,
                            Caption = NULL,
                            AxisLabelSize = 12,
                            AxisTextSize = 10,
                            StripTextSize = 11,
                            WrapLabels = TRUE,
                            LabelWrapWidth = 22,
                            LabelRadiusMultiplier = 1.22,
                            PlotMarginTop = 40,
                            PlotMarginRight = 120,
                            PlotMarginBottom = 40,
                            PlotMarginLeft = 120,
                            MakeInteractive = FALSE,
                            InteractiveHeight = 700,
                            InteractiveWidth = NULL,
                            InteractiveAxisMin = NULL,
                            InteractiveAxisMax = NULL,
                            TooltipDigits = 2) {

  # Validate inputs

  ContinuousSummary <- match.arg(ContinuousSummary, c("mean", "median"))
  ContinuousScaling <- match.arg(ContinuousScaling, c("zscore", "none", "minmax"))
  VariableOrder <- match.arg(
    VariableOrder,
    c("input", "discrimination", "hierarchical", "greedy", "category_discrimination")
  )

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }

  if (missing(Variables) || length(Variables) == 0) {
    stop("Variables must be a non-empty character vector.")
  }

  if (!is.character(Variables)) {
    stop("Variables must be a character vector of column names.")
  }

  if (!all(Variables %in% names(Data))) {
    stop(
      "The following Variables are not in Data: ",
      paste(Variables[!Variables %in% names(Data)], collapse = ", ")
    )
  }

  if (!is.null(GroupVariable) && !GroupVariable %in% names(Data)) {
    stop("GroupVariable must be a column in Data.")
  }

  if (!is.null(VariableCategories) && length(VariableCategories) != length(Variables)) {
    stop("VariableCategories must be NULL or the same length as Variables.")
  }

  if (!is.logical(Relabel) || length(Relabel) != 1) {
    stop("Relabel must be TRUE or FALSE.")
  }

  if (!is.logical(Fill) || length(Fill) != 1) {
    stop("Fill must be TRUE or FALSE.")
  }

  if (!is.logical(Facet) || length(Facet) != 1) {
    stop("Facet must be TRUE or FALSE.")
  }

  if (!is.logical(ShowPoints) || length(ShowPoints) != 1) {
    stop("ShowPoints must be TRUE or FALSE.")
  }

  if (!is.logical(WrapLabels) || length(WrapLabels) != 1) {
    stop("WrapLabels must be TRUE or FALSE.")
  }

  if (!is.logical(MakeInteractive) || length(MakeInteractive) != 1) {
    stop("MakeInteractive must be TRUE or FALSE.")
  }

  if (!is.numeric(FillAlpha) || length(FillAlpha) != 1 || FillAlpha < 0 || FillAlpha > 1) {
    stop("FillAlpha must be a single number between 0 and 1.")
  }

  if (!is.numeric(LineSize) || length(LineSize) != 1 || LineSize <= 0) {
    stop("LineSize must be a single positive number.")
  }

  if (!is.numeric(PointSize) || length(PointSize) != 1 || PointSize <= 0) {
    stop("PointSize must be a single positive number.")
  }

  if (!is.numeric(LabelWrapWidth) || length(LabelWrapWidth) != 1 || LabelWrapWidth <= 0) {
    stop("LabelWrapWidth must be a single positive number.")
  }

  if (!is.numeric(LabelRadiusMultiplier) || length(LabelRadiusMultiplier) != 1 || LabelRadiusMultiplier <= 1) {
    stop("LabelRadiusMultiplier must be a single number greater than 1.")
  }

  if (!is.null(InteractiveAxisMin) &&
      (!is.numeric(InteractiveAxisMin) || length(InteractiveAxisMin) != 1)) {
    stop("InteractiveAxisMin must be NULL or a single number.")
  }

  if (!is.null(InteractiveAxisMax) &&
      (!is.numeric(InteractiveAxisMax) || length(InteractiveAxisMax) != 1)) {
    stop("InteractiveAxisMax must be NULL or a single number.")
  }

  if (!is.numeric(TooltipDigits) || length(TooltipDigits) != 1 || TooltipDigits < 0) {
    stop("TooltipDigits must be a single non-negative number.")
  }

  if (Facet && is.null(GroupVariable)) {
    stop("Facet = TRUE requires a GroupVariable.")
  }

  if (is.null(LegendTitle) && !is.null(GroupVariable)) {
    LegendTitle <- GroupVariable
  }

  if (is.null(LegendTitle) && is.null(GroupVariable)) {
    LegendTitle <- "Group"
  }

  if (MakeInteractive && !requireNamespace("plotly", quietly = TRUE)) {
    stop("MakeInteractive = TRUE requires the plotly package to be installed.")
  }

  if (MakeInteractive && Facet) {
    warning("Facet is ignored when MakeInteractive = TRUE.")
  }

  # Prepare data

  PlotData <- Data %>%
    dplyr::select(dplyr::all_of(unique(c(GroupVariable, Variables)))) %>%
    dplyr::mutate(
      .Group = if (is.null(GroupVariable)) {
        factor("Overall")
      } else {
        as.factor(.data[[GroupVariable]])
      }
    )

  TypeTable <- tibble::tibble(
    Variable = Variables,
    VariableCategory = if (is.null(VariableCategories)) {
      rep(NA_character_, length(Variables))
    } else {
      as.character(VariableCategories)
    },
    Type = NA_character_,
    PositiveValue = NA_character_,
    SpokeLabel = Variables
  )

  for (i in seq_along(Variables)) {
    CurrentVar <- Variables[i]
    x <- PlotData[[CurrentVar]]

    if (inherits(x, "haven_labelled") || inherits(x, "labelled")) {
      x_for_type <- labelled::to_factor(x)

      if (length(unique(stats::na.omit(x_for_type))) > 2) {
        x_for_type <- suppressWarnings(as.numeric(as.character(x)))
      }
    } else {
      x_for_type <- x
    }

    NonMissing <- stats::na.omit(x_for_type)
    UniqueValues <- unique(NonMissing)

    IsBinary <- FALSE

    if (is.logical(x_for_type)) {
      IsBinary <- TRUE
    }

    if (is.factor(x_for_type) && length(levels(x_for_type)) == 2) {
      IsBinary <- TRUE
    }

    if (!is.factor(x_for_type) && length(UniqueValues) == 2) {
      IsBinary <- TRUE
    }

    if (IsBinary) {
      TypeTable$Type[i] <- "binary"

      if (is.factor(x_for_type)) {
        TypeTable$PositiveValue[i] <- as.character(levels(x_for_type)[2])
      } else if (is.logical(x_for_type)) {
        TypeTable$PositiveValue[i] <- "TRUE"
      } else if (BinaryPositiveValue %in% UniqueValues) {
        TypeTable$PositiveValue[i] <- as.character(BinaryPositiveValue)
      } else {
        TypeTable$PositiveValue[i] <- as.character(sort(as.character(UniqueValues))[2])
      }
    } else {
      if (!is.numeric(x) && !is.integer(x)) {
        stop(
          "Variable ", CurrentVar, " is neither binary nor numeric. ",
          "Spider charts currently support numeric continuous variables and binary variables."
        )
      }

      TypeTable$Type[i] <- "continuous"
    }

    if (Relabel) {
      CurrentLabel <- sjlabelled::get_label(Data[[CurrentVar]], def.value = CurrentVar)

      if (length(CurrentLabel) > 1) {
        CurrentLabel <- CurrentLabel[1]
      }

      if (is.null(CurrentLabel) || is.na(CurrentLabel) || CurrentLabel == "") {
        CurrentLabel <- CurrentVar
      }

      TypeTable$SpokeLabel[i] <- as.character(CurrentLabel)
    }
  }

  ContinuousVars <- TypeTable %>%
    dplyr::filter(.data$Type == "continuous") %>%
    dplyr::pull(.data$Variable)

  if (length(ContinuousVars) > 0 && ContinuousScaling == "zscore") {
    ZResult <- tryCatch(
      CalcZScore(
        df = PlotData,
        variables = ContinuousVars,
        names_prefix = "Z_",
        RetainLabels = FALSE,
        RenameLabels = FALSE,
        center = TRUE,
        scale = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(ZResult)) {
      stop(
        "CalcZScore() call failed inside PlotSpiderChart(). ",
        "Please check the current CalcZScore() argument names in SciDataReportR."
      )
    }

    if (is.data.frame(ZResult)) {
      PlotData <- ZResult
    } else if (is.list(ZResult) && "DataWithZ" %in% names(ZResult)) {
      PlotData <- ZResult$DataWithZ
    } else if (is.list(ZResult) && "df" %in% names(ZResult)) {
      PlotData <- ZResult$df
    } else if (is.list(ZResult) && "Data" %in% names(ZResult)) {
      PlotData <- ZResult$Data
    } else {
      stop(
        "CalcZScore() returned an object format that PlotSpiderChart() does not recognize."
      )
    }

    # Reminder for future package work:
    # CalcZScore() currently uses mean centering and SD scaling.
    # If you later add robust scaling there, update this function so the
    # 'median' option can optionally use median centering and robust scaling.
  }

  if (length(ContinuousVars) > 0 && ContinuousScaling == "minmax") {
    for (CurrentVar in ContinuousVars) {
      x <- PlotData[[CurrentVar]]

      xmin <- suppressWarnings(min(x, na.rm = TRUE))
      xmax <- suppressWarnings(max(x, na.rm = TRUE))

      if (!is.finite(xmin) || !is.finite(xmax)) {
        PlotData[[paste0("MM_", CurrentVar)]] <- NA_real_
      } else if (xmax == xmin) {
        PlotData[[paste0("MM_", CurrentVar)]] <- 0
      } else {
        PlotData[[paste0("MM_", CurrentVar)]] <- (x - xmin) / (xmax - xmin)
      }
    }
  }

  # Build outputs

  SummaryList <- vector("list", length = length(Variables))

  for (i in seq_along(Variables)) {
    CurrentVar <- Variables[i]
    CurrentType <- TypeTable$Type[TypeTable$Variable == CurrentVar]
    CurrentPositive <- TypeTable$PositiveValue[TypeTable$Variable == CurrentVar]
    CurrentLabel <- TypeTable$SpokeLabel[TypeTable$Variable == CurrentVar]
    CurrentCategory <- TypeTable$VariableCategory[TypeTable$Variable == CurrentVar]

    if (CurrentType == "continuous") {
      SummaryVar <- CurrentVar

      if (ContinuousScaling == "zscore") {
        SummaryVar <- paste0("Z_", CurrentVar)
      }

      if (ContinuousScaling == "minmax") {
        SummaryVar <- paste0("MM_", CurrentVar)
      }

      TempSummary <- PlotData %>%
        dplyr::group_by(.data$.Group) %>%
        dplyr::summarise(
          Value = if (ContinuousSummary == "mean") {
            mean(.data[[SummaryVar]], na.rm = TRUE)
          } else {
            stats::median(.data[[SummaryVar]], na.rm = TRUE)
          },
          NonMissingN = sum(!is.na(.data[[SummaryVar]])),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Variable = CurrentVar,
          SpokeLabel = CurrentLabel,
          VariableCategory = CurrentCategory,
          Type = CurrentType
        )
    } else {
      TempSummary <- PlotData %>%
        dplyr::group_by(.data$.Group) %>%
        dplyr::summarise(
          Value = {
            x <- .data[[CurrentVar]]

            if (inherits(x, "haven_labelled") || inherits(x, "labelled")) {
              x <- labelled::to_factor(x)
            }

            mean(as.character(x) == as.character(CurrentPositive), na.rm = TRUE) * 100
          },
          NonMissingN = sum(!is.na(.data[[CurrentVar]])),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Variable = CurrentVar,
          SpokeLabel = CurrentLabel,
          VariableCategory = CurrentCategory,
          Type = CurrentType
        )
    }

    TempSummary$Value[TempSummary$NonMissingN == 0] <- NA_real_
    SummaryList[[i]] <- TempSummary
  }

  SummaryData <- dplyr::bind_rows(SummaryList)

  DropVars <- SummaryData %>%
    dplyr::group_by(.data$Variable) %>%
    dplyr::summarise(AllMissingInAnyGroup = any(is.na(.data$Value)), .groups = "drop") %>%
    dplyr::filter(.data$AllMissingInAnyGroup) %>%
    dplyr::pull(.data$Variable)

  if (length(DropVars) > 0) {
    SummaryData <- SummaryData %>%
      dplyr::filter(!.data$Variable %in% DropVars)

    TypeTable <- TypeTable %>%
      dplyr::filter(!.data$Variable %in% DropVars)
  }

  if (nrow(SummaryData) == 0) {
    stop("No variables remained after dropping spokes with all-missing data in at least one group.")
  }

  FinalVariables <- TypeTable$Variable

  if (VariableOrder == "discrimination") {
    OrderTable <- SummaryData %>%
      dplyr::group_by(.data$Variable) %>%
      dplyr::summarise(
        Discrimination = if (all(is.na(.data$Value))) {
          NA_real_
        } else {
          diff(range(.data$Value, na.rm = TRUE))
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(InputOrder = match(.data$Variable, Variables)) %>%
      dplyr::arrange(dplyr::desc(.data$Discrimination), .data$InputOrder)

    FinalVariables <- OrderTable$Variable
  }

  if (VariableOrder == "hierarchical") {
    WideSummary <- SummaryData %>%
      dplyr::select("Variable", ".Group", "Value") %>%
      tidyr::pivot_wider(names_from = ".Group", values_from = "Value")

    WideMat <- WideSummary %>%
      tibble::column_to_rownames("Variable") %>%
      as.matrix()

    if (nrow(WideMat) > 1) {
      WideMat[is.na(WideMat)] <- 0
      HC <- stats::hclust(stats::dist(WideMat))
      FinalVariables <- rownames(WideMat)[HC$order]
    }
  }

  if (VariableOrder == "greedy") {
    WideSummary <- SummaryData %>%
      dplyr::select("Variable", ".Group", "Value") %>%
      tidyr::pivot_wider(names_from = ".Group", values_from = "Value")

    WideMat <- WideSummary %>%
      tibble::column_to_rownames("Variable") %>%
      as.matrix()

    if (nrow(WideMat) > 1) {
      WideMat[is.na(WideMat)] <- 0
      DistMat <- as.matrix(stats::dist(WideMat))
      RangeVec <- apply(WideMat, 1, function(x) diff(range(x, na.rm = TRUE)))
      OrderedVars <- names(sort(RangeVec, decreasing = TRUE))[1]
      RemainingVars <- setdiff(rownames(WideMat), OrderedVars)

      while (length(RemainingVars) > 0) {
        LastVar <- OrderedVars[length(OrderedVars)]
        NextVar <- RemainingVars[which.max(DistMat[LastVar, RemainingVars])]
        OrderedVars <- c(OrderedVars, NextVar)
        RemainingVars <- setdiff(RemainingVars, NextVar)
      }

      FinalVariables <- OrderedVars
    }
  }

  if (VariableOrder == "category_discrimination") {
    if (is.null(VariableCategories)) {
      stop("VariableOrder = 'category_discrimination' requires VariableCategories.")
    }

    OrderTable <- SummaryData %>%
      dplyr::group_by(.data$Variable, .data$VariableCategory) %>%
      dplyr::summarise(
        Discrimination = if (all(is.na(.data$Value))) {
          NA_real_
        } else {
          diff(range(.data$Value, na.rm = TRUE))
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        CategoryIndex = match(.data$VariableCategory, unique(VariableCategories)),
        VariableIndex = match(.data$Variable, Variables)
      ) %>%
      dplyr::arrange(.data$CategoryIndex, dplyr::desc(.data$Discrimination), .data$VariableIndex)

    FinalVariables <- OrderTable$Variable
  }

  SummaryData <- SummaryData %>%
    dplyr::mutate(
      Variable = factor(.data$Variable, levels = FinalVariables),
      SpokeLabel = factor(
        .data$SpokeLabel,
        levels = TypeTable$SpokeLabel[match(FinalVariables, TypeTable$Variable)]
      )
    ) %>%
    dplyr::arrange(.data$.Group, .data$Variable)

  if (!is.null(Subtitle)) {
    FinalSubtitle <- Subtitle
  } else if (length(unique(TypeTable$Type[TypeTable$Variable %in% FinalVariables])) > 1) {
    FinalSubtitle <- "Continuous variables are summarized on the selected scale; binary variables are shown as percent positive."
  } else {
    FinalSubtitle <- NULL
  }

  GroupLevels <- levels(as.factor(SummaryData$.Group))
  PlotColors <- grDevices::hcl.colors(length(GroupLevels), palette = Palette)
  names(PlotColors) <- GroupLevels

  # Return interactive widget

  if (MakeInteractive) {
    InteractiveValues <- SummaryData$Value[is.finite(SummaryData$Value)]

    if (length(InteractiveValues) == 0) {
      stop("No finite summarized values available for the interactive plot.")
    }

    InteractiveMinAuto <- min(InteractiveValues, na.rm = TRUE)
    InteractiveMaxAuto <- max(InteractiveValues, na.rm = TRUE)

    if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "binary")) {
      InteractiveMinAuto <- 0
      InteractiveMaxAuto <- max(100, InteractiveMaxAuto, na.rm = TRUE)
    }

    if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "continuous") &&
        ContinuousScaling == "minmax") {
      InteractiveMinAuto <- 0
      InteractiveMaxAuto <- max(1, InteractiveMaxAuto, na.rm = TRUE)
    }

    if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "continuous") &&
        ContinuousScaling == "none" &&
        InteractiveMinAuto >= 0) {
      InteractiveMinAuto <- 0
    }

    if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "continuous") &&
        ContinuousScaling == "zscore") {
      InteractiveMinAuto <- floor(InteractiveMinAuto)
      InteractiveMaxAuto <- ceiling(InteractiveMaxAuto)
    }

    if (InteractiveMinAuto == InteractiveMaxAuto) {
      InteractiveMinAuto <- InteractiveMinAuto - 0.5
      InteractiveMaxAuto <- InteractiveMaxAuto + 0.5
    }

    InteractiveRange <- InteractiveMaxAuto - InteractiveMinAuto

    InteractiveAxisMinFinal <- if (is.null(InteractiveAxisMin)) {
      InteractiveMinAuto
    } else {
      InteractiveAxisMin
    }

    InteractiveAxisMaxFinal <- if (is.null(InteractiveAxisMax)) {
      InteractiveMaxAuto + 0.05 * InteractiveRange
    } else {
      InteractiveAxisMax
    }

    if (InteractiveAxisMinFinal >= InteractiveAxisMaxFinal) {
      stop("InteractiveAxisMin must be less than InteractiveAxisMax.")
    }

    PlotlyData <- SummaryData %>%
      dplyr::mutate(
        SpokeLabelWrapped = if (WrapLabels) {
          stringr::str_wrap(as.character(.data$SpokeLabel), width = LabelWrapWidth)
        } else {
          as.character(.data$SpokeLabel)
        },
        Variable = factor(.data$Variable, levels = FinalVariables)
      ) %>%
      dplyr::arrange(.data$.Group, .data$Variable)

    ThetaLevels <- unique(if (WrapLabels) {
      stringr::str_wrap(TypeTable$SpokeLabel[match(FinalVariables, TypeTable$Variable)], width = LabelWrapWidth)
    } else {
      TypeTable$SpokeLabel[match(FinalVariables, TypeTable$Variable)]
    })

    p_int <- plotly::plot_ly()

    for (CurrentGroup in GroupLevels) {
      GroupData <- PlotlyData %>%
        dplyr::filter(.data$.Group == CurrentGroup) %>%
        dplyr::arrange(.data$Variable)

      ThetaVals <- GroupData$SpokeLabelWrapped
      RVals <- GroupData$Value
      RawVarNames <- as.character(GroupData$Variable)

      ThetaClosed <- c(as.character(ThetaVals), as.character(ThetaVals[1]))
      RClosed <- c(RVals, RVals[1])
      RawVarClosed <- c(RawVarNames, RawVarNames[1])

      HoverText <- paste0(
        LegendTitle, ": ", CurrentGroup,
        "<br>Variable: ", ThetaClosed,
        "<br>Raw name: ", RawVarClosed,
        "<br>Value: ", format(round(RClosed, TooltipDigits), nsmall = TooltipDigits)
      )

      p_int <- p_int %>%
        plotly::add_trace(
          type = "scatterpolar",
          mode = if (ShowPoints) "lines+markers" else "lines",
          r = RClosed,
          theta = ThetaClosed,
          name = CurrentGroup,
          line = list(
            color = PlotColors[[CurrentGroup]],
            width = LineSize * 2
          ),
          marker = list(
            color = PlotColors[[CurrentGroup]],
            size = PointSize * 3
          ),
          fill = if (Fill) "toself" else "none",
          fillcolor = grDevices::adjustcolor(PlotColors[[CurrentGroup]], alpha.f = FillAlpha),
          text = HoverText,
          hoverinfo = "text"
        )
    }

    InteractiveAnnotations <- NULL

    if (!is.null(Caption) && Caption != "") {
      InteractiveAnnotations <- list(
        list(
          text = Caption,
          x = 0,
          y = -0.12,
          xref = "paper",
          yref = "paper",
          xanchor = "left",
          yanchor = "top",
          showarrow = FALSE,
          font = list(size = 11, color = "gray40")
        )
      )
    }

    p_int <- p_int %>%
      plotly::layout(
        title = list(
          text = paste(
            c(
              PlotTitle,
              if (!is.null(FinalSubtitle) && FinalSubtitle != "") paste0("<br><sup>", FinalSubtitle, "</sup>")
            ),
            collapse = ""
          )
        ),
        polar = list(
          radialaxis = list(
            visible = TRUE,
            range = c(InteractiveAxisMinFinal, InteractiveAxisMaxFinal)
          ),
          angularaxis = list(
            direction = "clockwise",
            rotation = 90,
            categoryorder = "array",
            categoryarray = ThetaLevels,
            tickfont = list(size = AxisLabelSize)
          )
        ),
        legend = list(
          title = list(text = LegendTitle),
          orientation = "v"
        ),
        annotations = InteractiveAnnotations,
        width = InteractiveWidth,
        height = InteractiveHeight
      )

    return(p_int)
  }

  # Build static plot data

  if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "continuous") &&
      ContinuousScaling == "zscore") {
    YBreaks <- c(-2, -1, 0, 1, 2)
  } else if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "binary")) {
    YBreaks <- c(0, 25, 50, 75, 100)
  } else if (all(TypeTable$Type[TypeTable$Variable %in% FinalVariables] == "continuous") &&
             ContinuousScaling == "minmax") {
    YBreaks <- c(0, 0.25, 0.5, 0.75, 1)
  } else {
    YRange <- range(SummaryData$Value, na.rm = TRUE)
    YBreaks <- pretty(YRange, n = 5)
  }

  YBreaks <- YBreaks[is.finite(YBreaks)]

  if (length(YBreaks) < 2) {
    if (all(is.na(SummaryData$Value))) {
      stop("Unable to determine axis breaks because all plotted values are missing.")
    } else {
      YBreaks <- pretty(range(SummaryData$Value, na.rm = TRUE), n = 5)
    }
  }

  MinBreak <- min(YBreaks, na.rm = TRUE)
  RadiusOffset <- ifelse(MinBreak < 0, abs(MinBreak), 0)

  SummaryData <- SummaryData %>%
    dplyr::mutate(
      PlotValue = .data$Value + RadiusOffset
    )

  PlotBreaks <- YBreaks + RadiusOffset
  MaxRadius <- max(PlotBreaks, na.rm = TRUE)

  XLabels <- TypeTable$SpokeLabel[match(FinalVariables, TypeTable$Variable)]
  n_vars <- length(FinalVariables)

  Angles <- seq(
    from = pi / 2,
    to = pi / 2 - 2 * pi,
    length.out = n_vars + 1
  )[1:n_vars]

  AxisMap <- tibble::tibble(
    Variable = FinalVariables,
    SpokeLabel = XLabels,
    Angle = Angles
  )

  SummaryDataPlot <- SummaryData %>%
    dplyr::left_join(AxisMap, by = c("Variable", "SpokeLabel")) %>%
    dplyr::mutate(
      x = .data$PlotValue * cos(.data$Angle),
      y = .data$PlotValue * sin(.data$Angle)
    )

  SummaryDataClosedPlot <- SummaryDataPlot %>%
    dplyr::group_by(.data$.Group) %>%
    dplyr::arrange(match(as.character(.data$Variable), FinalVariables), .by_group = TRUE) %>%
    dplyr::group_modify(~ dplyr::bind_rows(.x, .x[1, ])) %>%
    dplyr::ungroup()

  GridAngles <- seq(0, 2 * pi, length.out = 300)

  GridData <- tidyr::expand_grid(
    Radius = PlotBreaks,
    Angle = GridAngles
  ) %>%
    dplyr::mutate(
      x = .data$Radius * cos(.data$Angle),
      y = .data$Radius * sin(.data$Angle)
    )

  SpokeData <- AxisMap %>%
    dplyr::mutate(
      x = MaxRadius * cos(.data$Angle),
      y = MaxRadius * sin(.data$Angle),
      x0 = 0,
      y0 = 0
    )

  LabelRadius <- MaxRadius * LabelRadiusMultiplier

  LabelData <- AxisMap %>%
    dplyr::mutate(
      SpokeLabel = if (WrapLabels) {
        stringr::str_wrap(.data$SpokeLabel, width = LabelWrapWidth)
      } else {
        .data$SpokeLabel
      },
      x = LabelRadius * cos(.data$Angle),
      y = LabelRadius * sin(.data$Angle),
      hjust = dplyr::case_when(
        cos(.data$Angle) > 0.2 ~ 0,
        cos(.data$Angle) < -0.2 ~ 1,
        TRUE ~ 0.5
      ),
      vjust = dplyr::case_when(
        sin(.data$Angle) > 0.2 ~ 0,
        sin(.data$Angle) < -0.2 ~ 1,
        TRUE ~ 0.5
      )
    )

  AxisLabelData <- tibble::tibble(
    Radius = PlotBreaks,
    Label = YBreaks,
    x = 0,
    y = PlotBreaks
  )

  p <- ggplot2::ggplot()

  p <- p +
    ggplot2::geom_path(
      data = GridData,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$Radius),
      color = "grey85",
      linewidth = 0.5
    ) +
    ggplot2::geom_segment(
      data = SpokeData,
      ggplot2::aes(
        x = .data$x0,
        y = .data$y0,
        xend = .data$x,
        yend = .data$y
      ),
      color = "grey85",
      linewidth = 0.5
    )

  if (Fill) {
    p <- p +
      ggplot2::geom_polygon(
        data = SummaryDataClosedPlot,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$.Group,
          fill = .data$.Group
        ),
        alpha = FillAlpha,
        color = NA
      )
  }

  p <- p +
    ggplot2::geom_path(
      data = SummaryDataClosedPlot,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$.Group,
        color = .data$.Group
      ),
      linewidth = LineSize,
      lineend = "round"
    )

  if (ShowPoints) {
    p <- p +
      ggplot2::geom_point(
        data = SummaryDataPlot,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          color = .data$.Group
        ),
        size = PointSize
      )
  }

  p <- p +
    ggplot2::geom_text(
      data = LabelData,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$SpokeLabel,
        hjust = .data$hjust,
        vjust = .data$vjust
      ),
      size = AxisLabelSize / 3
    ) +
    ggplot2::geom_text(
      data = AxisLabelData,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$Label
      ),
      hjust = -0.2,
      size = AxisTextSize / 3,
      color = "grey40"
    ) +
    ggplot2::scale_color_manual(values = PlotColors) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = PlotTitle,
      subtitle = FinalSubtitle,
      caption = Caption,
      color = LegendTitle
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = if (Facet) "none" else "right",
      strip.text = ggplot2::element_text(size = StripTextSize),
      plot.margin = ggplot2::margin(
        t = PlotMarginTop,
        r = PlotMarginRight,
        b = PlotMarginBottom,
        l = PlotMarginLeft
      )
    )

  if (Fill) {
    p <- p +
      ggplot2::scale_fill_manual(values = PlotColors) +
      ggplot2::labs(fill = LegendTitle)
  }

  if (Facet) {
    p <- p +
      ggplot2::facet_wrap(~.Group)
  }

  # Return result

  p
}
