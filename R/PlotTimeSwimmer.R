#' Plot longitudinal swimmer timelines
#'
#' Create swimmer-style longitudinal timeline plots from long-format visit data.
#' Each row in the input data represents a subject visit or observation. The
#' function can visualize longitudinal state trajectories, visit timing, and
#' events over time.
#'
#' Time can be represented using dates or numeric visit values. Timelines can
#' optionally be normalized relative to each participant's first visit or the
#' first occurrence of a specified event.
#'
#' States can either be displayed as continuous intervals extending forward
#' until the next visit (`StateInterval = "forward"`) or shown only as visit
#' points (`StateInterval = "point"`).
#'
#' @param Data A data frame in long format with one row per visit or observation.
#' @param ID Character string naming the participant ID column.
#' @param Time Character string naming the time variable.
#' @param State Optional character string naming a state/group variable used for
#'   coloring timelines or visit points.
#' @param Event Optional character string naming a logical or binary event
#'   indicator variable.
#' @param EventType Optional character string naming an event category variable.
#'   Used for event point shapes/colors.
#' @param TimeScale Character. One of `"from_first"`, `"observed"`, or
#'   `"from_event"`.
#'
#'   `"observed"` uses raw observed time values.
#'
#'   `"from_first"` normalizes each subject relative to their first observed
#'   timepoint.
#'
#'   `"from_event"` normalizes each subject relative to the first occurrence of
#'   `EventReference`.
#' @param EventReference Optional event value used when
#'   `TimeScale = "from_event"`.
#' @param StateInterval Character. One of `"forward"` or `"point"`.
#'
#'   `"forward"` extends the current state forward until the next visit.
#'
#'   `"point"` only colors visit points without extending intervals.
#' @param Format Character. One of `"state_path"`, `"visit_points"`,
#'   `"event_rug"`, or `"minimal"`.
#' @param SortBy Character. One of `"duration"`, `"last_time"`,
#'   `"first_time"`, `"state"`, or `"id"`.
#' @param TimeUnit Character. One of `"auto"`, `"days"`, `"weeks"`,
#'   `"months"`, `"years"`, or `"visits"`.
#' @param Relabel Logical. If `TRUE`, use labels from `Codebook` or variable
#'   attributes when available.
#' @param Codebook Optional codebook data frame with columns `Variable` and
#'   `Label`.
#' @param LineWidth Numeric line width for swimmer segments. Default is `5`.
#' @param PointSize Numeric point size for visit/event points. Default is `2.5`.
#' @param Alpha Numeric alpha transparency for swimmer segments. Default is
#'   `0.9`.
#' @param BaseSize Base font size for the plot theme. Default is `13`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' df <- tibble::tibble(
#'   ID = c(1,1,1,2,2,2),
#'   Visit = c(0,6,12,0,6,12),
#'   Cluster = c("A","A","B","B","B","C")
#' )
#'
#' PlotTimeSwimmer(
#'   Data = df,
#'   ID = "ID",
#'   Time = "Visit",
#'   State = "Cluster"
#' )
#'
#' @export
PlotTimeSwimmer <- function(
    Data,
    ID,
    Time,
    State = NULL,
    Event = NULL,
    EventType = NULL,
    TimeScale = c("from_first", "observed", "from_event"),
    EventReference = NULL,
    StateInterval = c("forward", "point"),
    Format = c(
      "state_path",
      "visit_points",
      "event_rug",
      "minimal"
    ),
    SortBy = c(
      "duration",
      "last_time",
      "first_time",
      "state",
      "id"
    ),
    TimeUnit = c(
      "auto",
      "days",
      "weeks",
      "months",
      "years",
      "visits"
    ),
    Relabel = TRUE,
    Codebook = NULL,
    LineWidth = 5,
    PointSize = 2.5,
    Alpha = 0.9,
    BaseSize = 13
) {

  TimeScale <- match.arg(TimeScale)
  StateInterval <- match.arg(StateInterval)
  Format <- match.arg(Format)
  SortBy <- match.arg(SortBy)
  TimeUnit <- match.arg(TimeUnit)

  # Validate inputs

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }

  required_vars <- c(ID, Time)

  if (!is.null(State)) {
    required_vars <- c(required_vars, State)
  }

  if (!is.null(Event)) {
    required_vars <- c(required_vars, Event)
  }

  if (!is.null(EventType)) {
    required_vars <- c(required_vars, EventType)
  }

  missing_vars <- setdiff(required_vars, names(Data))

  if (length(missing_vars) > 0) {
    stop(
      "The following required variables were not found in Data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (!is.null(Codebook)) {

    if (!is.data.frame(Codebook)) {
      stop("Codebook must be a data frame.")
    }

    required_codebook_cols <- c("Variable", "Label")

    missing_codebook_cols <- setdiff(
      required_codebook_cols,
      names(Codebook)
    )

    if (length(missing_codebook_cols) > 0) {
      stop(
        "Codebook must contain columns Variable and Label. Missing: ",
        paste(missing_codebook_cols, collapse = ", ")
      )
    }
  }

  # Prepare data

  plot_data <- Data %>%
    dplyr::mutate(
      .ID = as.character(.data[[ID]]),
      .TimeRaw = .data[[Time]]
    )

  is_date_time <- inherits(plot_data$.TimeRaw, c("Date", "POSIXct", "POSIXt"))

  plot_data <- plot_data %>%
    dplyr::arrange(.data$.ID, .data$.TimeRaw)

  # Normalize time

  if (TimeScale == "observed") {

    plot_data <- plot_data %>%
      dplyr::mutate(
        .TimePlot = .data$.TimeRaw
      )
  }

  if (TimeScale == "from_first") {

    if (is_date_time) {

      plot_data <- plot_data %>%
        dplyr::group_by(.data$.ID) %>%
        dplyr::mutate(
          .TimePlot = as.numeric(
            difftime(
              .data$.TimeRaw,
              min(.data$.TimeRaw, na.rm = TRUE),
              units = "days"
            )
          )
        ) %>%
        dplyr::ungroup()

    } else {

      plot_data <- plot_data %>%
        dplyr::group_by(.data$.ID) %>%
        dplyr::mutate(
          .TimePlot = .data$.TimeRaw -
            min(.data$.TimeRaw, na.rm = TRUE)
        ) %>%
        dplyr::ungroup()
    }
  }

  if (TimeScale == "from_event") {

    if (is.null(Event)) {
      stop(
        "Event must be supplied when TimeScale = 'from_event'."
      )
    }

    if (is.null(EventReference)) {
      stop(
        "EventReference must be supplied when TimeScale = 'from_event'."
      )
    }

    reference_df <- plot_data %>%
      dplyr::filter(.data[[Event]] == EventReference) %>%
      dplyr::group_by(.data$.ID) %>%
      dplyr::summarise(
        .ReferenceTime = min(.data$.TimeRaw, na.rm = TRUE),
        .groups = "drop"
      )

    plot_data <- plot_data %>%
      dplyr::left_join(reference_df, by = ".ID")

    if (is_date_time) {

      plot_data <- plot_data %>%
        dplyr::mutate(
          .TimePlot = as.numeric(
            difftime(
              .data$.TimeRaw,
              .data$.ReferenceTime,
              units = "days"
            )
          )
        )

    } else {

      plot_data <- plot_data %>%
        dplyr::mutate(
          .TimePlot = .data$.TimeRaw -
            .data$.ReferenceTime
        )
    }
  }

  # Convert time units

  if (TimeUnit == "auto") {

    if (is_date_time) {
      TimeUnit <- "days"
    } else {
      TimeUnit <- "visits"
    }
  }

  if (TimeUnit == "weeks") {
    plot_data$.TimePlot <- plot_data$.TimePlot / 7
  }

  if (TimeUnit == "months") {
    plot_data$.TimePlot <- plot_data$.TimePlot / 30.4375
  }

  if (TimeUnit == "years") {
    plot_data$.TimePlot <- plot_data$.TimePlot / 365.25
  }

  # Create intervals

  interval_data <- plot_data %>%
    dplyr::group_by(.data$.ID) %>%
    dplyr::arrange(.data$.TimePlot, .by_group = TRUE) %>%
    dplyr::mutate(
      .Start = .data$.TimePlot,
      .End = dplyr::lead(.data$.TimePlot)
    ) %>%
    dplyr::ungroup()

  interval_data <- interval_data %>%
    dplyr::filter(!is.na(.data$.End))

  # Sorting

  ordering_df <- plot_data %>%
    dplyr::group_by(.data$.ID) %>%
    dplyr::summarise(
      FirstTime = min(.data$.TimePlot, na.rm = TRUE),
      LastTime = max(.data$.TimePlot, na.rm = TRUE),
      Duration = LastTime - FirstTime,
      FirstState = if (!is.null(State)) {
        as.character(first(.data[[State]]))
      } else {
        NA_character_
      },
      .groups = "drop"
    )

  if (SortBy == "duration") {
    ordering_df <- ordering_df %>%
      dplyr::arrange(dplyr::desc(.data$Duration))
  }

  if (SortBy == "last_time") {
    ordering_df <- ordering_df %>%
      dplyr::arrange(dplyr::desc(.data$LastTime))
  }

  if (SortBy == "first_time") {
    ordering_df <- ordering_df %>%
      dplyr::arrange(.data$FirstTime)
  }

  if (SortBy == "state") {
    ordering_df <- ordering_df %>%
      dplyr::arrange(.data$FirstState)
  }

  if (SortBy == "id") {
    ordering_df <- ordering_df %>%
      dplyr::arrange(.data$.ID)
  }

  ordered_ids <- ordering_df$.ID

  plot_data <- plot_data %>%
    dplyr::mutate(
      .IDPlot = factor(.data$.ID, levels = rev(ordered_ids))
    )

  interval_data <- interval_data %>%
    dplyr::mutate(
      .IDPlot = factor(.data$.ID, levels = rev(ordered_ids))
    )

  # Axis label

  x_label <- dplyr::case_when(
    TimeScale == "observed" & is_date_time ~ "Observed date",
    TimeScale == "observed" ~ "Observed time",
    TimeScale == "from_first" ~ paste0("Time from first visit (", TimeUnit, ")"),
    TimeScale == "from_event" ~ paste0(
      "Time from event (",
      TimeUnit,
      ")"
    ),
    TRUE ~ "Time"
  )

  # Build plot

  p <- ggplot2::ggplot()

  # Add swimmer intervals

  if (StateInterval == "forward" &&
      Format %in% c("state_path", "minimal")) {

    if (!is.null(State)) {

      p <- p +
        ggplot2::geom_segment(
          data = interval_data,
          ggplot2::aes(
            x = .data$.Start,
            xend = .data$.End,
            y = .data$.IDPlot,
            yend = .data$.IDPlot,
            color = .data[[State]]
          ),
          linewidth = LineWidth,
          alpha = Alpha,
          lineend = "round"
        )

    } else {

      p <- p +
        ggplot2::geom_segment(
          data = interval_data,
          ggplot2::aes(
            x = .data$.Start,
            xend = .data$.End,
            y = .data$.IDPlot,
            yend = .data$.IDPlot
          ),
          linewidth = LineWidth,
          alpha = Alpha,
          color = "steelblue",
          lineend = "round"
        )
    }
  }

  # Visit points

  if (Format %in% c(
    "visit_points",
    "state_path",
    "event_rug"
  )) {

    if (!is.null(State)) {

      p <- p +
        ggplot2::geom_point(
          data = plot_data,
          ggplot2::aes(
            x = .data$.TimePlot,
            y = .data$.IDPlot,
            color = .data[[State]]
          ),
          size = PointSize,
          alpha = 0.95
        )

    } else {

      p <- p +
        ggplot2::geom_point(
          data = plot_data,
          ggplot2::aes(
            x = .data$.TimePlot,
            y = .data$.IDPlot
          ),
          size = PointSize,
          alpha = 0.95,
          color = "steelblue"
        )
    }
  }

  # Event points

  if (!is.null(Event)) {

    event_data <- plot_data %>%
      dplyr::filter(!is.na(.data[[Event]]))

    if (!is.null(EventReference)) {
      event_data <- event_data %>%
        dplyr::filter(.data[[Event]] == EventReference)
    }

    if (nrow(event_data) > 0) {

      if (!is.null(EventType)) {

        p <- p +
          ggplot2::geom_point(
            data = event_data,
            ggplot2::aes(
              x = .data$.TimePlot,
              y = .data$.IDPlot,
              shape = .data[[EventType]]
            ),
            size = PointSize + 1,
            stroke = 1.2,
            fill = "white",
            color = "black"
          )

      } else {

        p <- p +
          ggplot2::geom_point(
            data = event_data,
            ggplot2::aes(
              x = .data$.TimePlot,
              y = .data$.IDPlot
            ),
            shape = 21,
            size = PointSize + 1,
            stroke = 1.2,
            fill = "white",
            color = "black"
          )
      }
    }
  }

  # Theme

  p <- p +
    ggplot2::labs(
      x = x_label,
      y = "Participant"
    ) +
    ggplot2::theme_bw(base_size = BaseSize) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      axis.text.y = ggplot2::element_text(size = BaseSize * 0.6)
    )

  return(p)
}
