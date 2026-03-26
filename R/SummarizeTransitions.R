
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
#'   `x_axis_type = "date"`.
#' @param participant_subset Optional vector of participant IDs to include.
#' @param max_participants Optional maximum number of participants to retain after
#'   ordering is applied.
#' @param order_participants_by Character string controlling participant order.
#'   Options are `"first_positive"`, `"first_transition"`, `"ever_positive"`,
#'   `"ever_positive_then_burden"`, `"input_order"`, `"n_visits"`,
#'   `"n_positive"`, and `"pct_positive"`.
#' @param x_axis_type Character string indicating whether longitudinal ordering
#'   should follow aligned visit number (`"visit"`) or actual date (`"date"`).
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
#'
#' @examples
#' toy_df <- tibble::tibble(
#'   ParticipantID = rep(paste0("P", 1:4), each = 4),
#'   VisitOrder = rep(1:4, times = 4),
#'   MetSBinary = c(
#'     0, 0, 1, 1,
#'     1, 1, 0, 0,
#'     0, 0, 0, 0,
#'     TRUE, TRUE, TRUE, TRUE
#'   )
#' )
#'
#' SummarizeTransitions(
#'   data = toy_df,
#'   id_var = ParticipantID,
#'   time_var = VisitOrder,
#'   status_var = MetSBinary
#' )
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
                                 x_axis_type = c("visit", "date")) {

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
    condition_summary = condition_summary
  )
}
