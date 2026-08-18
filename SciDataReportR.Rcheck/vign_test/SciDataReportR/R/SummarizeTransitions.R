
#' Summarize participant transitions for a binary longitudinal condition
#'
#' Create participant-level and condition-level summary tables for a binary
#' condition observed across repeated visits. This function uses the same
#' transition logic and data preparation workflow as `PlotSwimmerTransitions()`
#' so that plotting and summary outputs remain aligned.
#'
#' @param data A data frame containing repeated observations per participant.
#' @param id_var Unquoted column name identifying the participant.
#' @param time_var Unquoted column name representing visit order, visit number, or time index.
#' @param status_var Unquoted column name representing the binary condition status.
#' @param date_var Optional unquoted visit date column. This is required when
#'   `x_axis_type = "date"` or `x_axis_type = "time_from_baseline"`.
#' @param participant_subset Optional vector of participant IDs to include.
#' @param max_participants Optional maximum number of participants to retain after
#'   ordering is applied.
#' @param order_participants_by Character string controlling participant order.
#'   Options are `"first_positive"`, `"first_transition"`, `"ever_positive"`,
#'   `"ever_positive_then_burden"`, `"input_order"`, `"n_visits"`,
#'   `"n_positive"`, and `"pct_positive"`.
#' @param x_axis_type Character string indicating whether longitudinal ordering
#'   should follow aligned visit number (`"visit"`), actual date (`"date"`), or
#'   elapsed time from baseline (`"time_from_baseline"`).
#' @param time_from_baseline_unit Character string specifying the unit for
#'   `x_axis_type = "time_from_baseline"`. Options are `"days"`, `"months"`,
#'   and `"years"`.
#'
#' @details
#' Transition rules are:
#' - `0 -> 1` = developed condition
#' - `1 -> 0` = resolved condition
#'
#' Missing values remain missing and are not recoded to 0.
#'
#' The returned condition-level summary includes:
#' - number of participants
#' - number ever positive
#' - number who developed the condition
#' - number who resolved the condition
#' - number sustained after development
#' - number sustained after resolution
#'
#' @return A list with:
#' - `participant_summary`: participant-level summary table
#' - `condition_summary`: one-row tibble with overall counts
#' - `Plots`: figures for the two summaries, described below
#'
#' @section Figures:
#' `condition_summary` is a single row of six counts, which is exactly the
#' shape a table reads worst: the numbers are related to each other, and the
#' relationship is the point. `Plots` renders them so the relationship is
#' visible:
#'
#' * `Plots$ConditionCascade` shows the counts as a cascade - the whole
#'   cohort, the part of it ever affected, the transitions that occurred
#'   within that part, and how many of those persisted - with each bar
#'   labelled by its count and its share of the cohort.
#' * `Plots$TransitionPatterns` classifies every participant into one
#'   mutually exclusive longitudinal pattern (never positive, positive
#'   throughout, developed only, resolved only, developed and resolved) and
#'   plots the counts. The indicator columns in `participant_summary` overlap,
#'   so this is the view that shows what the cohort is actually made of.
#'
#' The two agree by construction: the pattern counts sum to
#' `n_participants`, and everything except "never positive" sums to
#' `n_ever_positive`.
#'
#' For the per-participant trajectories behind these counts, use
#' [PlotSwimmerTransitions()], which prepares its data the same way.
#'
#' @seealso [PlotSwimmerTransitions()] for the participant-level swimmer plot.
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
#'     TRUE, TRUE, TRUE, TRUE
#'   )
#' )
#'
#' transitions <- SummarizeTransitions(
#'   data = toy_df,
#'   id_var = ParticipantID,
#'   time_var = VisitOrder,
#'   status_var = MetSBinary,
#'   date_var = VisitDate,
#'   x_axis_type = "time_from_baseline",
#'   time_from_baseline_unit = "months"
#' )
#'
#' transitions$condition_summary
#'
#' \donttest{
#' # A larger cohort, where a positive visit tends to be followed by another
#' set.seed(9)
#' n_participants <- 60
#'
#' df_Longitudinal <- do.call(rbind, lapply(seq_len(n_participants), function(i) {
#'   n_visits <- sample(3:5, 1)
#'   status <- numeric(n_visits)
#'   status[1] <- stats::rbinom(1, 1, 0.3)
#'   for (j in seq_len(n_visits)[-1]) {
#'     status[j] <- stats::rbinom(1, 1, if (status[j - 1] == 1) 0.8 else 0.25)
#'   }
#'   data.frame(
#'     ParticipantID = sprintf("P%02d", i),
#'     VisitOrder = seq_len(n_visits),
#'     MetSBinary = status
#'   )
#' }))
#'
#' transitions <- SummarizeTransitions(
#'   data = df_Longitudinal,
#'   id_var = ParticipantID,
#'   time_var = VisitOrder,
#'   status_var = MetSBinary
#' )
#'
#' # Six counts in one row
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(transitions$condition_summary, full_width = TRUE)
#' )))
#'
#' # The same six numbers as a cascade
#' transitions$Plots$ConditionCascade
#'
#' # The cohort split into mutually exclusive patterns
#' transitions$Plots$TransitionPatterns
#'
#' # The participant-level table
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(
#'     dplyr::select(
#'       transitions$participant_summary,
#'       .plot_id, n_visits, n_positive, pct_positive, ever_positive,
#'       developed_condition, resolved_condition
#'     ),
#'     height = "300px", full_width = TRUE
#'   )
#' )))
#' }
#' @export
SummarizeTransitions <- function(data,
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
                                 x_axis_type = c("visit", "date", "time_from_baseline"),
                                 time_from_baseline_unit = c("days", "months", "years")) {

  order_participants_by <- match.arg(order_participants_by)
  x_axis_type <- match.arg(x_axis_type)
  time_from_baseline_unit <- match.arg(time_from_baseline_unit)

  prepared <- .PrepareTransitionData(
    data = data,
    id_var = {{ id_var }},
    time_var = {{ time_var }},
    status_var = {{ status_var }},
    date_var = {{ date_var }},
    participant_subset = participant_subset,
    max_participants = max_participants,
    order_participants_by = order_participants_by,
    x_axis_type = x_axis_type,
    time_from_baseline_unit = time_from_baseline_unit
  )

  participant_summary <- prepared$participant_summary

  condition_summary <- participant_summary %>%
    dplyr::summarise(
      n_participants = dplyr::n(),
      n_ever_positive = sum(ever_positive, na.rm = TRUE),
      n_developed_condition = sum(developed_condition, na.rm = TRUE),
      n_resolved_condition = sum(resolved_condition, na.rm = TRUE),
      n_sustained_after_development = sum(sustained_after_development, na.rm = TRUE),
      n_sustained_after_resolution = sum(sustained_after_resolution, na.rm = TRUE)
    )

  list(
    participant_summary = participant_summary,
    condition_summary = condition_summary,
    Plots = list(
      ConditionCascade = .PlotTransitionCascade(condition_summary),
      TransitionPatterns = .PlotTransitionPatterns(participant_summary)
    )
  )
}

#' Cascade figure for a transition condition summary
#'
#' Renders the one-row `condition_summary` as a horizontal bar chart, which is
#' the shape the counts actually have: a cohort, the part of it ever affected,
#' the transitions observed within that part, and how many of those persisted.
#'
#' @param condition_summary The one-row tibble built by [SummarizeTransitions()].
#' @return A ggplot object.
#' @noRd
.PlotTransitionCascade <- function(condition_summary) {
  n_participants <- condition_summary$n_participants[[1]]

  df_Cascade <- dplyr::tibble(
    Stage = factor(
      c(
        "Participants",
        "Ever positive",
        "Developed (0 to 1)",
        "Resolved (1 to 0)",
        "Sustained after developing",
        "Sustained after resolving"
      ),
      levels = rev(c(
        "Participants",
        "Ever positive",
        "Developed (0 to 1)",
        "Resolved (1 to 0)",
        "Sustained after developing",
        "Sustained after resolving"
      ))
    ),
    Category = factor(
      c("Cohort", "Ever affected", "Transition", "Transition", "Sustained", "Sustained"),
      levels = c("Cohort", "Ever affected", "Transition", "Sustained")
    ),
    Count = c(
      n_participants,
      condition_summary$n_ever_positive[[1]],
      condition_summary$n_developed_condition[[1]],
      condition_summary$n_resolved_condition[[1]],
      condition_summary$n_sustained_after_development[[1]],
      condition_summary$n_sustained_after_resolution[[1]]
    )
  )

  df_Cascade$Percent <- if (isTRUE(n_participants > 0)) {
    100 * df_Cascade$Count / n_participants
  } else {
    rep(NA_real_, nrow(df_Cascade))
  }

  df_Cascade$Label <- ifelse(
    is.na(df_Cascade$Percent),
    as.character(df_Cascade$Count),
    sprintf("%d (%.0f%%)", as.integer(df_Cascade$Count), df_Cascade$Percent)
  )

  axis_max <- max(c(df_Cascade$Count, 1), na.rm = TRUE)

  ggplot2::ggplot(
    df_Cascade,
    ggplot2::aes(x = .data$Count, y = .data$Stage, fill = .data$Category)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Label),
      hjust = -0.12,
      size = 3.2
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.22)),
      limits = c(0, axis_max)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Cohort" = "#B8C4CE",
        "Ever affected" = "#7BA7BC",
        "Transition" = "#E0A458",
        "Sustained" = "#A3B18A"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Condition transitions across the cohort",
      subtitle = "Percentages are of all participants",
      x = "Participants",
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Participant transition-pattern figure
#'
#' Classifies every participant into one mutually exclusive longitudinal
#' pattern and plots the counts, so the cohort's composition is visible rather
#' than having to be reconstructed from overlapping indicator columns.
#'
#' @param participant_summary The participant-level tibble built by
#'   [SummarizeTransitions()].
#' @return A ggplot object.
#' @noRd
.PlotTransitionPatterns <- function(participant_summary) {
  ever_positive <- as.logical(participant_summary$ever_positive)
  developed <- as.logical(participant_summary$developed_condition)
  resolved <- as.logical(participant_summary$resolved_condition)

  ever_positive[is.na(ever_positive)] <- FALSE
  developed[is.na(developed)] <- FALSE
  resolved[is.na(resolved)] <- FALSE

  pattern <- dplyr::case_when(
    !ever_positive ~ "Never positive",
    developed & resolved ~ "Developed and resolved",
    developed ~ "Developed only",
    resolved ~ "Resolved only",
    TRUE ~ "Positive throughout"
  )

  pattern_levels <- c(
    "Never positive",
    "Positive throughout",
    "Developed only",
    "Resolved only",
    "Developed and resolved"
  )

  df_Patterns <- dplyr::tibble(
    Pattern = factor(pattern, levels = pattern_levels)
  ) %>%
    dplyr::count(.data$Pattern, .drop = FALSE)

  n_total <- sum(df_Patterns$n)

  df_Patterns$Label <- if (isTRUE(n_total > 0)) {
    sprintf("%d (%.0f%%)", df_Patterns$n, 100 * df_Patterns$n / n_total)
  } else {
    as.character(df_Patterns$n)
  }

  axis_max <- max(c(df_Patterns$n, 1), na.rm = TRUE)

  ggplot2::ggplot(
    df_Patterns,
    ggplot2::aes(x = .data$n, y = stats::reorder(.data$Pattern, dplyr::row_number(.data$Pattern)))
  ) +
    ggplot2::geom_col(fill = "#7BA7BC", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), hjust = -0.12, size = 3.2) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.22)),
      limits = c(0, axis_max)
    ) +
    ggplot2::labs(
      title = "Participant longitudinal patterns",
      subtitle = "Each participant appears exactly once",
      x = "Participants",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
