# Internal helper:
# Recode common binary encodings to integer 0/1 while preserving NA.
.RecodeBinaryCondition <- function(x, var_name = "status_var") {

  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_integer_, ifelse(x, 1L, 0L)))
  }

  if (is.numeric(x)) {
    non_missing <- unique(x[!is.na(x)])

    if (!all(non_missing %in% c(0, 1))) {
      stop(
        "`", var_name, "` must contain only binary values when numeric. ",
        "Found: ", paste(sort(non_missing), collapse = ", ")
      )
    }

    return(as.integer(x))
  }

  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)
  x_chr_lower <- stringr::str_to_lower(x_chr)

  true_values <- c("1", "yes", "y", "true", "t", "present", "positive", "pos")
  false_values <- c("0", "no", "n", "false", "f", "absent", "negative", "neg")

  out <- rep(NA_integer_, length(x_chr_lower))

  out[!is.na(x_chr_lower) & x_chr_lower %in% true_values] <- 1L
  out[!is.na(x_chr_lower) & x_chr_lower %in% false_values] <- 0L

  unknown_values <- unique(x_chr[!is.na(x_chr) & is.na(out)])

  if (length(unknown_values) > 0) {
    stop(
      "`", var_name, "` contains values that could not be interpreted as binary: ",
      paste(unknown_values, collapse = ", ")
    )
  }

  out
}

# Internal helper:
# Prepare visit-level and participant-level transition data used by plotting and summary functions.
.PrepareTransitionData <- function(data,
                                   id_var,
                                   time_var,
                                   status_var,
                                   date_var = NULL,
                                   participant_subset = NULL,
                                   max_participants = NULL,
                                   order_participants_by = c(
                                     "first_positive",
                                     "first_transition",
                                     "ever_positive",
                                     "ever_positive_then_burden",
                                     "input_order",
                                     "n_visits",
                                     "n_positive",
                                     "pct_positive"
                                   ),
                                   x_axis_type = c("visit", "date")) {

  order_participants_by <- match.arg(order_participants_by)
  x_axis_type <- match.arg(x_axis_type)

  id_quo <- rlang::enquo(id_var)
  time_quo <- rlang::enquo(time_var)
  status_quo <- rlang::enquo(status_var)
  date_quo <- rlang::enquo(date_var)

  has_date_var <- !rlang::quo_is_null(date_quo)

  # Validate inputs

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  if (x_axis_type == "date" && !has_date_var) {
    stop("`date_var` must be provided when `x_axis_type = \"date\"`.")
  }

  # Prepare data

  data_prepped <- data %>%
    dplyr::mutate(
      .plot_id = as.character(!!id_quo),
      .plot_time = !!time_quo,
      .plot_status_raw = !!status_quo,
      .input_row_order = dplyr::row_number()
    )

  if (has_date_var) {
    data_prepped <- data_prepped %>%
      dplyr::mutate(.plot_date = !!date_quo)
  } else {
    data_prepped <- data_prepped %>%
      dplyr::mutate(.plot_date = as.Date(NA))
  }

  if (any(is.na(data_prepped$.plot_id))) {
    stop("`id_var` contains missing values. Participant IDs must be non-missing.")
  }

  if (x_axis_type == "visit" && any(is.na(data_prepped$.plot_time))) {
    stop("`time_var` contains missing values. Visit order/time must be non-missing.")
  }

  if (x_axis_type == "date" && any(is.na(data_prepped$.plot_date))) {
    stop("`date_var` contains missing values when `x_axis_type = \"date\"`.")
  }

  if (!is.null(participant_subset)) {
    data_prepped <- data_prepped %>%
      dplyr::filter(.plot_id %in% as.character(participant_subset))
  }

  if (nrow(data_prepped) == 0) {
    stop("No rows remain after filtering.")
  }

  data_prepped <- data_prepped %>%
    dplyr::mutate(
      .plot_status = .RecodeBinaryCondition(
        x = .plot_status_raw,
        var_name = rlang::as_name(status_quo)
      )
    )

  if (has_date_var) {
    data_prepped <- data_prepped %>%
      dplyr::arrange(.plot_id, .plot_date, .plot_time, .input_row_order)
  } else {
    data_prepped <- data_prepped %>%
      dplyr::arrange(.plot_id, .plot_time, .input_row_order)
  }

  plot_data <- data_prepped %>%
    dplyr::group_by(.plot_id) %>%
    dplyr::mutate(
      .visit_index = dplyr::row_number(),
      .lag_status = dplyr::lag(.plot_status),
      .transition = dplyr::case_when(
        is.na(.lag_status) | is.na(.plot_status) ~ NA_character_,
        .lag_status == 0L & .plot_status == 1L ~ "Developed",
        .lag_status == 1L & .plot_status == 0L ~ "Resolved",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup()

  participant_summary <- plot_data %>%
    dplyr::group_by(.plot_id) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_visits = sum(!is.na(.plot_status)),
      n_positive = sum(.plot_status == 1L, na.rm = TRUE),
      pct_positive = dplyr::if_else(n_visits > 0, n_positive / n_visits, NA_real_),
      ever_positive = any(.plot_status == 1L, na.rm = TRUE),
      developed_condition = any(.transition == "Developed", na.rm = TRUE),
      resolved_condition = any(.transition == "Resolved", na.rm = TRUE),
      first_positive_time = {
        vals <- .plot_time[.plot_status == 1L]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else NA_real_
      },
      first_positive_date = {
        vals <- .plot_date[.plot_status == 1L]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else as.Date(NA)
      },
      first_transition_time = {
        vals <- .plot_time[!is.na(.transition)]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else NA_real_
      },
      first_transition_date = {
        vals <- .plot_date[!is.na(.transition)]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else as.Date(NA)
      },
      first_developed_time = {
        vals <- .plot_time[.transition == "Developed"]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else NA_real_
      },
      first_resolved_time = {
        vals <- .plot_time[.transition == "Resolved"]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else NA_real_
      },
      first_developed_date = {
        vals <- .plot_date[.transition == "Developed"]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else as.Date(NA)
      },
      first_resolved_date = {
        vals <- .plot_date[.transition == "Resolved"]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) min(vals) else as.Date(NA)
      },
      input_order = min(.input_row_order, na.rm = TRUE),
      .groups = "drop"
    )

  sustained_after_development <- plot_data %>%
    dplyr::left_join(
      participant_summary %>%
        dplyr::select(.plot_id, developed_condition, first_developed_time),
      by = ".plot_id"
    ) %>%
    dplyr::group_by(.plot_id) %>%
    dplyr::summarise(
      sustained_after_development = {
        has_dev <- dplyr::first(developed_condition)
        dev_time <- dplyr::first(first_developed_time)

        if (!has_dev || is.na(dev_time)) {
          FALSE
        } else {
          later_status <- .plot_status[.plot_time > dev_time]
          later_status <- later_status[!is.na(later_status)]
          length(later_status) > 0 && all(later_status == 1L)
        }
      },
      .groups = "drop"
    )

  sustained_after_resolution <- plot_data %>%
    dplyr::left_join(
      participant_summary %>%
        dplyr::select(.plot_id, resolved_condition, first_resolved_time),
      by = ".plot_id"
    ) %>%
    dplyr::group_by(.plot_id) %>%
    dplyr::summarise(
      sustained_after_resolution = {
        has_res <- dplyr::first(resolved_condition)
        res_time <- dplyr::first(first_resolved_time)

        if (!has_res || is.na(res_time)) {
          FALSE
        } else {
          later_status <- .plot_status[.plot_time > res_time]
          later_status <- later_status[!is.na(later_status)]
          length(later_status) > 0 && all(later_status == 0L)
        }
      },
      .groups = "drop"
    )

  participant_summary <- participant_summary %>%
    dplyr::left_join(sustained_after_development, by = ".plot_id") %>%
    dplyr::left_join(sustained_after_resolution, by = ".plot_id")

  participant_summary <- participant_summary %>%
    dplyr::mutate(
      .order_group = dplyr::case_when(
        order_participants_by == "ever_positive" & ever_positive ~ 0,
        order_participants_by == "ever_positive" & !ever_positive ~ 1,
        order_participants_by == "ever_positive_then_burden" & ever_positive ~ 0,
        order_participants_by == "ever_positive_then_burden" & !ever_positive ~ 1,
        TRUE ~ 0
      ),
      .order_value = dplyr::case_when(
        order_participants_by == "first_positive" & !is.na(first_positive_time) ~ first_positive_time,
        order_participants_by == "first_positive" & is.na(first_positive_time) ~ Inf,
        order_participants_by == "first_transition" & !is.na(first_transition_time) ~ first_transition_time,
        order_participants_by == "first_transition" & is.na(first_transition_time) ~ Inf,
        order_participants_by == "ever_positive" ~ as.numeric(.order_group),
        order_participants_by == "ever_positive_then_burden" ~ -as.numeric(dplyr::coalesce(pct_positive, -Inf)),
        order_participants_by == "input_order" ~ as.numeric(input_order),
        order_participants_by == "n_visits" ~ -as.numeric(n_visits),
        order_participants_by == "n_positive" ~ -as.numeric(n_positive),
        order_participants_by == "pct_positive" ~ -as.numeric(dplyr::coalesce(pct_positive, -Inf)),
        TRUE ~ as.numeric(input_order)
      )
    ) %>%
    dplyr::arrange(.order_group, .order_value, input_order, .plot_id)

  if (!is.null(max_participants)) {
    if (!is.numeric(max_participants) || length(max_participants) != 1 || max_participants < 1) {
      stop("`max_participants` must be a single positive number when provided.")
    }

    participant_summary <- participant_summary %>%
      dplyr::slice_head(n = max_participants)
  }

  keep_ids <- participant_summary$.plot_id

  participant_summary <- participant_summary %>%
    dplyr::mutate(
      .plot_id_factor = factor(.plot_id, levels = rev(keep_ids))
    )

  plot_data <- plot_data %>%
    dplyr::filter(.plot_id %in% keep_ids) %>%
    dplyr::left_join(
      participant_summary %>%
        dplyr::select(
          .plot_id,
          .plot_id_factor,
          n_rows,
          n_visits,
          n_positive,
          pct_positive,
          ever_positive,
          developed_condition,
          resolved_condition,
          sustained_after_development,
          sustained_after_resolution,
          first_positive_time,
          first_positive_date,
          first_transition_time,
          first_transition_date
        ),
      by = ".plot_id"
    ) %>%
    dplyr::mutate(
      .status_label = dplyr::case_when(
        is.na(.plot_status) ~ "Missing",
        .plot_status == 0L ~ "Absent",
        .plot_status == 1L ~ "Present"
      ),
      .status_label = factor(.status_label, levels = c("Absent", "Present", "Missing")),
      .tooltip_label = paste0(
        "ID: ", .plot_id,
        "\nVisit: ", .plot_time,
        ifelse(!is.na(.plot_date), paste0("\nDate: ", .plot_date), ""),
        "\nStatus: ", as.character(.status_label),
        ifelse(!is.na(.transition), paste0("\nTransition: ", .transition), ""),
        "\nN observed visits: ", n_visits,
        "\nN positive visits: ", n_positive,
        "\nPct positive: ",
        ifelse(is.na(pct_positive), "NA", paste0(round(100 * pct_positive, 1), "%"))
      )
    )

  list(
    plot_data = plot_data,
    participant_summary = participant_summary,
    x_axis_type = x_axis_type
  )
}

#' Plot swimmer-style transitions for a binary condition over repeated visits
#'
#' Create a swimmer-style longitudinal plot for a binary condition measured across
#' repeated visits. The plot shows condition status at each visit, highlights
#' transition points where the condition develops or resolves, and supports
#' multiple participant ordering strategies for exploratory quality control,
#' longitudinal debugging, or figure creation.
#'
#' @param data A data frame containing repeated observations per participant.
#' @param id_var Unquoted column name identifying the participant.
#' @param time_var Unquoted column name representing visit order, visit number, or time index.
#' @param status_var Unquoted column name representing the binary condition status.
#'   Accepted encodings include numeric 0/1, logical TRUE/FALSE, factor values
#'   such as `"Yes"` and `"No"`, and character values such as `"1"` and `"0"`.
#' @param date_var Optional unquoted visit date column. This is required when
#'   `x_axis_type = "date"`.
#' @param participant_subset Optional vector of participant IDs to include.
#' @param max_participants Optional maximum number of participants to display
#'   after ordering is applied.
#' @param order_participants_by Character string controlling participant order in
#'   the plot. Options are `"first_positive"`, `"first_transition"`,
#'   `"ever_positive"`, `"ever_positive_then_burden"`, `"input_order"`,
#'   `"n_visits"`, `"n_positive"`, and `"pct_positive"`.
#' @param x_axis_type Character string indicating whether the x-axis should use
#'   aligned visit number (`"visit"`) or actual calendar date (`"date"`).
#' @param show_transition_points Logical. If `TRUE`, transition visits are highlighted.
#' @param show_lines Logical. If `TRUE`, a swimmer line is drawn across visits
#'   within each participant.
#' @param show_y_axis_labels Logical. If `TRUE`, show participant labels on the
#'   y-axis. Defaults to `FALSE`.
#' @param make_interactive Logical. If `TRUE`, return an interactive plotly object.
#'   Defaults to `FALSE`.
#' @param plot_title Optional custom plot title.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label. Defaults to `NULL`.
#' @param return_data Logical. If `TRUE`, return both the plot and the processed data.
#'
#' @details
#' Visits are sorted within participant by `date_var` if provided, otherwise by
#' `time_var`.
#'
#' Transition rules are defined as:
#' - `0 -> 1` = `"Developed"`
#' - `1 -> 0` = `"Resolved"`
#'
#' Missing status values remain `NA` and are not forced to 0.
#'
#' Participant-level summaries are created internally, including:
#' - whether the participant was ever positive
#' - first positive visit
#' - first transition visit
#' - number of observed visits
#' - number and proportion of positive visits
#' - whether positivity is sustained after first onset
#' - whether resolution is sustained after first resolution
#'
#' Sustained after development means that all non-missing observations after the
#' first `0 -> 1` transition remain positive.
#'
#' Sustained after resolution means that all non-missing observations after the
#' first `1 -> 0` transition remain negative.
#'
#' When `x_axis_type = "visit"`, participants are aligned by visit order.
#' When `x_axis_type = "date"`, visits are shown at actual calendar dates and do
#' not need to align across participants.
#'
#' When `make_interactive = TRUE`, the function returns a `plotly` object and
#' uses the internally prepared tooltip text for hover labels.
#'
#' @return A `ggplot` object by default.
#'
#' If `make_interactive = TRUE`, returns a `plotly` object.
#'
#' If `return_data = TRUE`, returns a list with:
#' - `plot`: the `ggplot` or `plotly` object
#' - `plot_data`: the processed visit-level plotting data
#' - `participant_summary`: the participant-level summary table
#'
#' @examples
#' toy_df <- tibble::tibble(
#'   ParticipantID = rep(paste0("P", 1:4), each = 4),
#'   VisitOrder = rep(1:4, times = 4),
#'   VisitDate = rep(seq.Date(as.Date("2024-01-01"), by = "month", length.out = 4), times = 4),
#'   MetSBinary = c(
#'     0, 0, 1, 1,
#'     1, 1, 0, 0,
#'     0, 0, 0, 0,
#'     "No", "Yes", "Yes", "Yes"
#'   )
#' )
#'
#' PlotSwimmerTransitions(
#'   data = toy_df,
#'   id_var = ParticipantID,
#'   time_var = VisitOrder,
#'   status_var = MetSBinary
#' )
#'
#' PlotSwimmerTransitions(
#'   data = toy_df,
#'   id_var = ParticipantID,
#'   time_var = VisitOrder,
#'   status_var = MetSBinary,
#'   date_var = VisitDate,
#'   x_axis_type = "date"
#' )
#' @export
PlotSwimmerTransitions <- function(data,
                                   id_var,
                                   time_var,
                                   status_var,
                                   date_var = NULL,
                                   participant_subset = NULL,
                                   max_participants = NULL,
                                   order_participants_by = c(
                                     "first_positive",
                                     "first_transition",
                                     "ever_positive",
                                     "ever_positive_then_burden",
                                     "input_order",
                                     "n_visits",
                                     "n_positive",
                                     "pct_positive"
                                   ),
                                   x_axis_type = c("visit", "date"),
                                   show_transition_points = TRUE,
                                   show_lines = TRUE,
                                   show_y_axis_labels = FALSE,
                                   make_interactive = FALSE,
                                   plot_title = NULL,
                                   x_label = NULL,
                                   y_label = NULL,
                                   return_data = FALSE) {

  order_participants_by <- match.arg(order_participants_by)
  x_axis_type <- match.arg(x_axis_type)

  prepared <- .PrepareTransitionData(
    data = data,
    id_var = {{ id_var }},
    time_var = {{ time_var }},
    status_var = {{ status_var }},
    date_var = {{ date_var }},
    participant_subset = participant_subset,
    max_participants = max_participants,
    order_participants_by = order_participants_by,
    x_axis_type = x_axis_type
  )

  plot_data <- prepared$plot_data
  participant_summary <- prepared$participant_summary

  # Validate inputs

  if (nrow(plot_data) == 0) {
    stop("No data available to plot after preprocessing.")
  }

  if (make_interactive && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package `plotly` must be installed when `make_interactive = TRUE`.")
  }

  # Build outputs

  if (is.null(plot_title)) {
    plot_title <- paste(
      "Swimmer plot of",
      rlang::as_name(rlang::ensym(status_var)),
      "across visits"
    )
  }

  if (is.null(x_label)) {
    x_label <- if (x_axis_type == "visit") {
      rlang::as_name(rlang::ensym(time_var))
    } else {
      rlang::as_name(rlang::ensym(date_var))
    }
  }

  if (x_axis_type == "visit") {
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .plot_time,
        y = .plot_id_factor
      )
    )
  } else {
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .plot_date,
        y = .plot_id_factor
      )
    )
  }

  if (show_lines) {
    if (make_interactive) {
      p <- p +
        ggplot2::geom_line(
          ggplot2::aes(group = .plot_id, text = .tooltip_label),
          alpha = 0.25,
          linewidth = 0.3,
          na.rm = TRUE
        )
    } else {
      p <- p +
        ggplot2::geom_line(
          ggplot2::aes(group = .plot_id),
          alpha = 0.25,
          linewidth = 0.3,
          na.rm = TRUE
        )
    }
  }

  if (make_interactive) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .status_label, text = .tooltip_label),
        size = 2.4,
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .status_label),
        size = 2.4,
        na.rm = TRUE
      )
  }

  if (show_transition_points) {
    transition_data <- plot_data %>%
      dplyr::filter(!is.na(.transition))

    if (nrow(transition_data) > 0) {
      if (make_interactive) {
        p <- p +
          ggplot2::geom_point(
            data = transition_data,
            ggplot2::aes(shape = .transition, text = .tooltip_label),
            size = 4,
            stroke = 1.1,
            colour = "black",
            fill = "white",
            na.rm = TRUE
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            data = transition_data,
            ggplot2::aes(shape = .transition),
            size = 4,
            stroke = 1.1,
            colour = "black",
            fill = "white",
            na.rm = TRUE
          )
      }
    }
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = c(
        "Absent" = "#E64B35",
        "Present" = "#00A087",
        "Missing" = "grey70"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Developed" = 17,
        "Resolved" = 25
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      color = "Status",
      shape = "Transition"
    ) +
    ggplot2::theme_bw()

  if (!show_y_axis_labels) {
    p <- p +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  if (x_axis_type == "date") {
    p <- p + ggplot2::scale_x_date()
  }

  if (make_interactive) {
    p <- plotly::ggplotly(p, tooltip = "text")
  }

  # Return result

  if (return_data) {
    return(list(
      plot = p,
      plot_data = plot_data,
      participant_summary = participant_summary
    ))
  }

  p
}
