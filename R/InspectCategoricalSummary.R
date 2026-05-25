#' Inspect categorical variables
#'
#' Provides a SciDataReportR-native categorical inspection summary and plot
#' inspired by `inspectdf::inspect_cat()`, created because inspectdf is archived
#' and no longer available on CRAN.
#'
#' @param Data A data frame.
#' @param Variables Optional character vector of categorical variables to
#'   summarize. If `NULL`, character, factor, logical, labelled, and
#'   haven-labelled columns are detected automatically.
#' @param Codebook Optional data frame with `Variable` and `Label` columns.
#' @param IncludeMissing Logical. If `TRUE`, missing values are included as a
#'   level and in percentages.
#' @param MissingLabel Character label to display for missing values.
#' @param RetainLabels Logical. If `TRUE`, use variable labels from `Codebook`
#'   and value labels from labelled variables when available.
#' @param SortLevelsBy One of `"Frequency"`, `"Value"`, or `"None"`.
#' @param SortVariablesBy One of `"Input"`, `"Label"`, `"MissingPercent"`,
#'   `"UniqueLevels"`, or `"TotalN"`.
#' @param Descending Logical. If `TRUE`, sort selected summaries in descending
#'   order.
#' @param MaxLevels Positive integer giving the maximum number of levels to
#'   display per variable in the plot before collapsing lower-ranked levels to
#'   `"(Other)"`.
#' @param Plot Logical. If `TRUE`, return a ggplot object.
#' @param PlotType One of `"bar"` or `"lollipop"`.
#' @param FacetScales One of `"free_y"` or `"fixed"`.
#' @param UsePercent Logical. If `TRUE`, plot percentages; otherwise plot
#'   counts.
#' @param LabelBars Logical. If `TRUE`, add readable labels to bars or points.
#' @param WrapLabels Integer number of characters used to wrap displayed level
#'   and facet labels.
#' @param BaseSize Base font size passed to `theme_minimal()`.
#'
#' @return A named list with `Summary`, a tibble containing categorical counts
#'   and percentages, and `Plot`, a ggplot object when `Plot = TRUE` or `NULL`
#'   when `Plot = FALSE`.
#' @references Rushworth A. inspectdf: Inspection, comparison and visualisation
#'   of data frames. Formerly available on CRAN; archived package.
#' @examples
#' df <- data.frame(
#'   group = factor(c("A", "B", "A", NA), levels = c("A", "B", "C")),
#'   flag = c(TRUE, FALSE, TRUE, NA),
#'   text = c("low", "high", "low", "mid")
#' )
#'
#' result <- InspectCategoricalSummary(df, Plot = FALSE)
#' result$Summary
#'
#' plot_result <- InspectCategoricalSummary(df, Variables = "group", Plot = TRUE)
#' plot_result$Plot
#' @importFrom magrittr %>%
#' @export
InspectCategoricalSummary <- function(Data,
                                      Variables = NULL,
                                      Codebook = NULL,
                                      IncludeMissing = TRUE,
                                      MissingLabel = "(Missing)",
                                      RetainLabels = TRUE,
                                      SortLevelsBy = c("Frequency", "Value", "None"),
                                      SortVariablesBy = c("Input", "Label", "MissingPercent", "UniqueLevels", "TotalN"),
                                      Descending = TRUE,
                                      MaxLevels = 30,
                                      Plot = TRUE,
                                      PlotType = c("bar", "lollipop"),
                                      FacetScales = c("free_y", "fixed"),
                                      UsePercent = TRUE,
                                      LabelBars = TRUE,
                                      WrapLabels = 35,
                                      BaseSize = 11) {

  if (!is.data.frame(Data)) {
    stop("`Data` must be a data frame.", call. = FALSE)
  }

  SortLevelsBy <- match.arg(SortLevelsBy)
  SortVariablesBy <- match.arg(SortVariablesBy)
  PlotType <- match.arg(PlotType)
  FacetScales <- match.arg(FacetScales)

  if (!is.logical(IncludeMissing) || length(IncludeMissing) != 1 || is.na(IncludeMissing)) {
    stop("`IncludeMissing` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(MissingLabel) || length(MissingLabel) != 1 || is.na(MissingLabel)) {
    stop("`MissingLabel` must be a single character value.", call. = FALSE)
  }
  if (!is.logical(RetainLabels) || length(RetainLabels) != 1 || is.na(RetainLabels)) {
    stop("`RetainLabels` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(Descending) || length(Descending) != 1 || is.na(Descending)) {
    stop("`Descending` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(MaxLevels) != 1 || is.na(MaxLevels) || MaxLevels < 1 || MaxLevels != as.integer(MaxLevels)) {
    stop("`MaxLevels` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(Plot) || length(Plot) != 1 || is.na(Plot)) {
    stop("`Plot` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(UsePercent) || length(UsePercent) != 1 || is.na(UsePercent)) {
    stop("`UsePercent` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(LabelBars) || length(LabelBars) != 1 || is.na(LabelBars)) {
    stop("`LabelBars` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(WrapLabels) != 1 || is.na(WrapLabels) || WrapLabels < 1) {
    stop("`WrapLabels` must be a positive number.", call. = FALSE)
  }
  if (length(BaseSize) != 1 || is.na(BaseSize) || BaseSize <= 0) {
    stop("`BaseSize` must be a positive number.", call. = FALSE)
  }
  if (Plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` must be installed when `Plot = TRUE`.", call. = FALSE)
  }

  if (!is.null(Codebook)) {
    if (!is.data.frame(Codebook)) {
      stop("`Codebook` must be a data frame.", call. = FALSE)
    }
    if (!all(c("Variable", "Label") %in% colnames(Codebook))) {
      stop("`Codebook` must contain `Variable` and `Label` columns.", call. = FALSE)
    }
  }

  if (is.null(Variables)) {
    Variables <- names(Data)[vapply(Data, sdr_is_categorical, logical(1))]
  }

  if (!is.character(Variables)) {
    stop("`Variables` must be a character vector.", call. = FALSE)
  }

  missing_vars <- setdiff(Variables, names(Data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in `Data`: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  variable_labels <- sdr_variable_labels(Data, Variables, Codebook, RetainLabels)
  variable_order <- sdr_variable_order(
    Data = Data,
    Variables = Variables,
    Labels = variable_labels,
    SortVariablesBy = SortVariablesBy,
    Descending = Descending,
    IncludeMissing = IncludeMissing
  )

  summary <- purrr::map_dfr(Variables, function(variable) {
    x <- Data[[variable]]
    missing <- is.na(x)
    value_labels <- sdr_value_labels(x, RetainLabels)
    raw_level <- as.character(x)

    if (is.factor(x)) {
      non_missing_levels <- levels(x)[levels(x) %in% raw_level[!missing]]
    } else {
      non_missing_levels <- unique(raw_level[!missing])
    }

    tibble::tibble(
      Variable = variable,
      Label = unname(variable_labels[[variable]]),
      Level = raw_level,
      Missing = missing
    ) %>%
      dplyr::filter(IncludeMissing | !.data$Missing) %>%
      dplyr::mutate(
        Level = dplyr::if_else(.data$Missing, MissingLabel, .data$Level),
        LevelLabel = dplyr::if_else(
          .data$Missing,
          MissingLabel,
          sdr_apply_value_labels(.data$Level, value_labels)
        )
      ) %>%
      dplyr::count(.data$Variable, .data$Label, .data$Level, .data$LevelLabel, .data$Missing, name = "N") %>%
      dplyr::mutate(
        TotalN = if (IncludeMissing) length(x) else sum(!missing),
        NonMissingN = sum(!missing),
        MissingN = sum(missing),
        MissingPercent = dplyr::if_else(length(x) == 0, NA_real_, sum(missing) / length(x)),
        Percent = dplyr::if_else(.data$TotalN == 0, NA_real_, .data$N / .data$TotalN),
        UniqueLevels = length(non_missing_levels),
        VariableClass = paste(class(x), collapse = "/"),
        VariableOrder = match(variable, variable_order),
        FactorOrder = match(.data$Level, non_missing_levels)
      )
  })

  summary <- summary %>%
    dplyr::group_by(.data$Variable)

  if (SortLevelsBy == "Frequency") {
    summary <- if (Descending) {
      summary %>% dplyr::arrange(.data$Missing, dplyr::desc(.data$N), .data$LevelLabel, .by_group = TRUE)
    } else {
      summary %>% dplyr::arrange(.data$Missing, .data$N, .data$LevelLabel, .by_group = TRUE)
    }
  } else if (SortLevelsBy == "Value") {
    summary <- if (Descending) {
      summary %>% dplyr::arrange(.data$Missing, dplyr::desc(.data$LevelLabel), .by_group = TRUE)
    } else {
      summary %>% dplyr::arrange(.data$Missing, .data$LevelLabel, .by_group = TRUE)
    }
  } else {
    summary <- summary %>% dplyr::arrange(.data$Missing, .data$FactorOrder, .by_group = TRUE)
  }

  summary <- summary %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$VariableOrder, .data$Missing) %>%
    dplyr::select(
      dplyr::all_of(c(
        "Variable", "Label", "Level", "LevelLabel", "N", "Percent",
        "Missing", "TotalN", "NonMissingN", "MissingN", "MissingPercent",
        "UniqueLevels", "VariableClass"
      ))
    ) %>%
    tibble::as_tibble()

  plot <- NULL
  if (Plot) {
    plot <- sdr_plot_categorical_summary(
      summary = summary,
      variable_order = variable_order,
      IncludeMissing = IncludeMissing,
      MaxLevels = MaxLevels,
      PlotType = PlotType,
      FacetScales = FacetScales,
      UsePercent = UsePercent,
      LabelBars = LabelBars,
      WrapLabels = WrapLabels,
      BaseSize = BaseSize
    )
  }

  list(Summary = summary, Plot = plot)
}

# Detects categorical-like columns without requiring haven at runtime.
sdr_is_categorical <- function(x) {
  is.character(x) ||
    is.factor(x) ||
    is.logical(x) ||
    inherits(x, "labelled") ||
    inherits(x, "haven_labelled")
}

sdr_variable_labels <- function(Data, Variables, Codebook, RetainLabels) {
  labels <- stats::setNames(Variables, Variables)

  if (RetainLabels) {
    data_labels <- vapply(Variables, function(variable) {
      label <- attr(Data[[variable]], "label", exact = TRUE)
      if (is.null(label) || is.na(label) || identical(label, "")) {
        variable
      } else {
        as.character(label)
      }
    }, character(1))
    labels <- stats::setNames(data_labels, Variables)
  }

  if (RetainLabels && !is.null(Codebook)) {
    codebook_labels <- Codebook %>%
      dplyr::select(dplyr::all_of(c("Variable", "Label"))) %>%
      dplyr::mutate(
        Variable = as.character(.data$Variable),
        Label = as.character(.data$Label)
      )

    matched <- match(Variables, codebook_labels$Variable)
    has_label <- !is.na(matched) & !is.na(codebook_labels$Label[matched]) & codebook_labels$Label[matched] != ""
    labels[has_label] <- codebook_labels$Label[matched[has_label]]
  }

  labels
}

sdr_variable_order <- function(Data, Variables, Labels, SortVariablesBy, Descending, IncludeMissing) {
  order_df <- tibble::tibble(
    Variable = Variables,
    Label = unname(Labels[Variables]),
    InputOrder = seq_along(Variables),
    MissingPercent = vapply(Variables, function(variable) {
      mean(is.na(Data[[variable]]))
    }, numeric(1)),
    UniqueLevels = vapply(Variables, function(variable) {
      length(unique(as.character(Data[[variable]][!is.na(Data[[variable]])])))
    }, integer(1)),
    TotalN = vapply(Variables, function(variable) {
      if (IncludeMissing) length(Data[[variable]]) else sum(!is.na(Data[[variable]]))
    }, integer(1))
  )

  if (SortVariablesBy == "Input") {
    return(Variables)
  }

  sort_column <- rlang::sym(SortVariablesBy)
  order_df <- if (Descending) {
    order_df %>% dplyr::arrange(dplyr::desc(!!sort_column), .data$InputOrder)
  } else {
    order_df %>% dplyr::arrange(!!sort_column, .data$InputOrder)
  }

  order_df$Variable
}

sdr_value_labels <- function(x, RetainLabels) {
  if (!RetainLabels) {
    return(NULL)
  }

  value_labels <- attr(x, "labels", exact = TRUE)
  if (is.null(value_labels) || length(value_labels) == 0 || is.null(names(value_labels))) {
    return(NULL)
  }

  stats::setNames(names(value_labels), as.character(unname(value_labels)))
}

sdr_apply_value_labels <- function(levels, value_labels) {
  if (is.null(value_labels)) {
    return(as.character(levels))
  }

  labels <- unname(value_labels[as.character(levels)])
  dplyr::if_else(is.na(labels), as.character(levels), labels)
}

sdr_plot_categorical_summary <- function(summary,
                                         variable_order,
                                         IncludeMissing,
                                         MaxLevels,
                                         PlotType,
                                         FacetScales,
                                         UsePercent,
                                         LabelBars,
                                         WrapLabels,
                                         BaseSize) {
  if (nrow(summary) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_minimal(base_size = BaseSize) +
        ggplot2::labs(
          x = NULL,
          y = NULL,
          caption = if (IncludeMissing) "Missing values included." else "Missing values excluded."
        )
    )
  }

  plot_df <- summary %>%
    dplyr::mutate(VariableOrder = match(.data$Variable, variable_order)) %>%
    dplyr::group_by(.data$Variable) %>%
    dplyr::arrange(.data$Missing, dplyr::desc(.data$N), .by_group = TRUE) %>%
    dplyr::mutate(LevelRank = dplyr::row_number()) %>%
    dplyr::mutate(
      LevelLabelPlot = dplyr::if_else(.data$LevelRank > MaxLevels, "(Other)", .data$LevelLabel),
      LevelPlot = dplyr::if_else(.data$LevelRank > MaxLevels, "(Other)", .data$Level)
    ) %>%
    dplyr::group_by(
      .data$Variable, .data$Label, .data$VariableOrder, .data$LevelPlot,
      .data$LevelLabelPlot
    ) %>%
    dplyr::summarise(
      N = sum(.data$N),
      TotalN = dplyr::first(.data$TotalN),
      Missing = any(.data$Missing),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Percent = dplyr::if_else(.data$TotalN == 0, NA_real_, .data$N / .data$TotalN),
      Value = if (UsePercent) .data$Percent else .data$N,
      FacetLabel = stringr::str_wrap(.data$Label, width = WrapLabels),
      AxisLabel = stringr::str_wrap(.data$LevelLabelPlot, width = WrapLabels),
      DisplayLabel = if (UsePercent) {
        scales::percent(.data$Percent, accuracy = 1)
      } else {
        as.character(.data$N)
      }
    ) %>%
    dplyr::group_by(.data$Variable) %>%
    dplyr::arrange(.data$Missing, .data$Value, .by_group = TRUE) %>%
    dplyr::mutate(AxisLabel = factor(.data$AxisLabel, levels = unique(.data$AxisLabel))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$VariableOrder)

  x_lab <- if (UsePercent) "Percent" else "Count"
  caption <- if (IncludeMissing) "Missing values included." else "Missing values excluded."

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Value, y = .data$AxisLabel)
  )

  if (PlotType == "bar") {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(fill = .data$Missing), width = 0.72, show.legend = FALSE)
  } else {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$Value, yend = .data$AxisLabel, color = .data$Missing),
        linewidth = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(ggplot2::aes(color = .data$Missing), size = 2.2, show.legend = FALSE)
  }

  if (LabelBars) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$DisplayLabel),
        hjust = -0.1,
        size = BaseSize / 3.5
      )
  }

  p <- p +
    ggplot2::facet_wrap(stats::as.formula("~ FacetLabel"), scales = FacetScales) +
    ggplot2::labs(
      x = x_lab,
      y = NULL,
      caption = caption
    ) +
    ggplot2::theme_minimal(base_size = BaseSize) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (PlotType == "bar") {
    p <- p + ggplot2::scale_fill_manual(values = c(`FALSE` = "grey40", `TRUE` = "grey70"))
  } else {
    p <- p + ggplot2::scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "grey65"))
  }

  if (UsePercent) {
    p <- p + ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    )
  }

  p
}
