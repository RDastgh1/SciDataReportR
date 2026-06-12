#' Merge fragmented records into a single observation
#'
#' Merges multiple rows representing fragments of the same observation
#' (e.g., participant visit, assessment session, study encounter, or
#' questionnaire administration) into a single row by selecting the
#' first non-missing value within each variable.
#'
#' This function is useful when data collection is interrupted and
#' restarted, resulting in multiple partial records for the same
#' participant and date. Examples include tablet-based cognitive testing,
#' mobile applications, REDCap surveys, wearable device uploads, and
#' electronic assessments where internet connectivity or software issues
#' may split a single visit across multiple records.
#'
#' Rows are grouped by the combination of `id_var` and `date_var`.
#' Within each group, observations are ordered by `session_var`, and
#' the first non-missing value encountered for each variable is retained.
#'
#' @param df A data frame containing fragmented records.
#' @param id_var Character string specifying the participant identifier
#'   variable. Default is `"subject"`.
#' @param date_var Character string specifying the visit or assessment
#'   date variable. Default is `"date"`.
#' @param session_var Character string specifying the session identifier
#'   variable used to order fragmented records. Default is `"session"`.
#' @param keep_session Logical. If `TRUE`, the first session value
#'   encountered within each group is retained. Default is `TRUE`.
#' @param session_name Character string specifying the name of the
#'   retained session variable. Default is `"first_session"`.
#' @param n_rows_name Character string specifying the name of the
#'   variable recording the number of rows merged. Default is
#'   `"n_rows_collapsed"`.
#' @param arrange_desc_session Logical. If `TRUE`, records are ordered
#'   by descending session number before merging. Default is `FALSE`.
#' @param empty_strings_to_na Logical. If `TRUE`, empty character strings
#'   are converted to missing values prior to merging. Default is `TRUE`.
#'
#' @return
#' A data frame containing one row per unique combination of
#' `id_var` and `date_var`.
#'
#' Additional variables may include:
#' \describe{
#'   \item{n_rows_collapsed}{Number of fragmented rows merged.}
#'   \item{first_session}{First session value retained, if
#'   `keep_session = TRUE`.}
#' }
#'
#' @details
#' For each variable within a participant-date grouping, the function
#' returns the first non-missing value after sorting by session order.
#'
#' For example:
#'
#' \preformatted{
#' subject   date         session   Stroop   Trails
#' 101       2024-01-01   1         50       NA
#' 101       2024-01-01   2         NA       100
#'
#' becomes
#'
#' subject   date         Stroop   Trails
#' 101       2024-01-01   50       100
#' }
#'
#' No attempt is made to resolve conflicting non-missing values across
#' sessions. If multiple non-missing values exist for the same variable,
#' the first value encountered after sorting is retained.
#' @export
MergeFragmentedRecords <- function(
  df,
  id_var = "subject",
  date_var = "date",
  session_var = "session",
  keep_session = TRUE,
  session_name = "first_session",
  n_rows_name = "n_rows_collapsed",
  arrange_desc_session = FALSE,
  empty_strings_to_na = TRUE
) {

  stopifnot(is.data.frame(df))

  required_vars <- c(id_var, date_var, session_var)

  missing_vars <- setdiff(required_vars, names(df))

  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from df: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  if (empty_strings_to_na) {
    df <- df %>%
      dplyr::mutate(
        dplyr::across(
          where(is.character),
          ~ dplyr::na_if(.x, "")
        )
      )
  }

  group_vars <- c(id_var, date_var)

  drop_vars <- c(id_var, date_var, session_var)

  collapse_vars <- setdiff(names(df), drop_vars)

  if (arrange_desc_session) {
    df <- df %>%
      dplyr::arrange(
        dplyr::across(dplyr::all_of(group_vars)),
        dplyr::desc(.data[[session_var]])
      )
  } else {
    df <- df %>%
      dplyr::arrange(
        dplyr::across(dplyr::all_of(group_vars)),
        .data[[session_var]]
      )
  }

  out <- df %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_vars))
    ) %>%
    dplyr::summarise(
      "{n_rows_name}" := dplyr::n(),
      dplyr::across(
        dplyr::all_of(collapse_vars),
        ~ {
          vals <- .x[!is.na(.x)]
          if (length(vals) == 0) NA else vals[[1]]
        }
      ),
      .groups = "drop"
    )

  if (keep_session) {

    session_df <- df %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(group_vars))
      ) %>%
      dplyr::summarise(
        "{session_name}" := {
          vals <- .data[[session_var]][!is.na(.data[[session_var]])]
          if (length(vals) == 0) NA else vals[[1]]
        },
        .groups = "drop"
      )

    out <- out %>%
      dplyr::left_join(
        session_df,
        by = group_vars
      ) %>%
      dplyr::relocate(
        dplyr::all_of(session_name),
        .after = dplyr::all_of(date_var)
      )
  }

  out
}
