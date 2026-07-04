#' Merge two data frames and validate the result in one step
#'
#' Perform a left join or closest-time merge and immediately audit the result
#' with [ValidateMerge()]. `safe_merge()` can also harmonize merge key types
#' before joining, so common numeric-vs-character key mismatches do not break
#' an otherwise valid merge.
#'
#' @param df_before A data frame. The left dataset that rows are added to.
#' @param df_add A data frame. The right dataset being merged in.
#' @param by Character vector of merge keys. Multiple keys are supported.
#' @param name Required. A short label for this merge.
#' @param method Either `"exact"` or `"closest_time"`.
#' @param time_var_before For `method = "closest_time"`, the time variable in
#'   `df_before`.
#' @param time_var_add For `method = "closest_time"`, the time variable in
#'   `df_add`.
#' @param min_match_rate Minimum acceptable key match rate before a merge is
#'   marked as `"WARNING"`.
#' @param harmonize_keys Logical. If `TRUE`, harmonize merge key types before
#'   joining.
#' @param key_parser Optional named list of parser functions for specific merge
#'   keys. Names should match values in `by`.
#' @param stop_on_failed_numeric Logical. If `TRUE`, stop when one key is
#'   numeric-like and the other contains values that cannot be safely converted
#'   to numeric.
#' @param ... Additional arguments passed to [Merge_ByClosestTime()] when
#'   `method = "closest_time"`.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{data}{The merged data frame.}
#'   \item{validation}{The full [ValidateMerge()] result, including
#'   `KeyHarmonization` and `SummaryTable`.}
#'   \item{log}{A one-row tibble with merge-level status and metrics.}
#'   \item{summary}{A compact formatted summary table.}
#' }
#'
#' @export
safe_merge <- function(df_before,
                       df_add,
                       by,
                       name,
                       method = c("exact", "closest_time"),
                       time_var_before = NULL,
                       time_var_add = NULL,
                       min_match_rate = 0.95,
                       harmonize_keys = TRUE,
                       key_parser = NULL,
                       stop_on_failed_numeric = TRUE,
                       ...) {

  method <- match.arg(method)

  if (!is.data.frame(df_before)) {
    stop("df_before must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(df_add)) {
    stop("df_add must be a data.frame.", call. = FALSE)
  }

  if (missing(by) || !is.character(by) || length(by) == 0) {
    stop("by must be a character vector of one or more merge keys.", call. = FALSE)
  }

  if (missing(name)) {
    stop(
      "name is required. Provide a short label for this merge, e.g. name = \"Baseline + labs\".",
      call. = FALSE
    )
  }

  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("name must be a single non-empty character value.", call. = FALSE)
  }

  if (!is.numeric(min_match_rate) ||
      length(min_match_rate) != 1 ||
      is.na(min_match_rate) ||
      min_match_rate < 0 ||
      min_match_rate > 1) {
    stop("min_match_rate must be a single number between 0 and 1.", call. = FALSE)
  }

  missing_before <- setdiff(by, names(df_before))
  missing_add <- setdiff(by, names(df_add))

  if (length(missing_before) > 0) {
    stop(
      "The following merge key(s) are missing from df_before: ",
      paste(missing_before, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_add) > 0) {
    stop(
      "The following merge key(s) are missing from df_add: ",
      paste(missing_add, collapse = ", "),
      call. = FALSE
    )
  }

  key_report <- tibble::tibble(
    Key = character(),
    LeftTypeBefore = character(),
    RightTypeBefore = character(),
    LeftTypeAfter = character(),
    RightTypeAfter = character(),
    ParserUsed = logical(),
    HarmonizedTo = character(),
    Action = character(),
    Status = character(),
    FailedLeftExamples = character(),
    FailedRightExamples = character()
  )

  df_before_join <- df_before
  df_add_join <- df_add

  if (isTRUE(harmonize_keys)) {
    key_harmonized <- HarmonizeMergeKeys(
      df_before = df_before,
      df_add = df_add,
      by = by,
      key_parser = key_parser,
      stop_on_failed_numeric = stop_on_failed_numeric
    )

    df_before_join <- key_harmonized$left
    df_add_join <- key_harmonized$right
    key_report <- key_harmonized$report
  }

  if (method == "exact") {
    df_after <- dplyr::left_join(
      df_before_join,
      df_add_join,
      by = by
    )
  } else {
    if (is.null(time_var_before) || is.null(time_var_add)) {
      stop(
        "method = \"closest_time\" requires both time_var_before and time_var_add.",
        call. = FALSE
      )
    }

    if (!time_var_before %in% names(df_before_join)) {
      stop(
        "time_var_before ('",
        time_var_before,
        "') was not found in df_before.",
        call. = FALSE
      )
    }

    if (!time_var_add %in% names(df_add_join)) {
      stop(
        "time_var_add ('",
        time_var_add,
        "') was not found in df_add.",
        call. = FALSE
      )
    }

    closest_result <- Merge_ByClosestTime(
      DataFrame1 = df_before_join,
      DataFrame2 = df_add_join,
      TimeVar1 = time_var_before,
      TimeVar2 = time_var_add,
      keys = by,
      ...
    )

    df_after <- closest_result$merged_dataframe
  }

  validation <- ValidateMerge(
    LeftData = df_before_join,
    RightData = df_add_join,
    MergedData = df_after,
    keys = by
  )

  summary_row <- validation$Summary

  duplicate_key_groups <- summary_row$DuplicateKeyGroups_Left +
    summary_row$DuplicateKeyGroups_Right +
    summary_row$DuplicateKeyGroups_Merged

  expected_cols_added <- ncol(df_add_join) - length(by)
  actual_cols_added <- ncol(df_after) - ncol(df_before_join)

  rows_changed <- nrow(df_after) != nrow(df_before_join)
  cols_wrong <- actual_cols_added != expected_cols_added
  duplicate_keys <- duplicate_key_groups > 0

  match_rate <- ifelse(
    summary_row$LeftUniqueKeys > 0,
    summary_row$MatchingKeys / summary_row$LeftUniqueKeys,
    NA_real_
  )

  merge_status <- dplyr::case_when(
    rows_changed ~ "FAIL",
    duplicate_keys ~ "FAIL",
    cols_wrong ~ "FAIL",
    !is.na(match_rate) && match_rate < min_match_rate ~ "WARNING",
    TRUE ~ "PASS"
  )

  ready_for_analysis <- merge_status != "FAIL"

  status_label <- dplyr::case_when(
    merge_status == "PASS" ~ "PASS",
    merge_status == "WARNING" ~ "WARNING",
    merge_status == "FAIL" ~ "FAIL",
    TRUE ~ merge_status
  )

  status_background <- dplyr::case_when(
    merge_status == "PASS" ~ "#2E7D32",
    merge_status == "WARNING" ~ "#F9A825",
    merge_status == "FAIL" ~ "#C62828",
    TRUE ~ "#616161"
  )

  key_harmonization_note <- if (nrow(key_report) == 0) {
    "No key harmonization performed."
  } else if (any(key_report$LeftTypeBefore != key_report$RightTypeBefore) ||
             any(key_report$ParserUsed)) {
    paste(
      paste0(
        key_report$Key,
        ": ",
        key_report$LeftTypeBefore,
        " / ",
        key_report$RightTypeBefore,
        " -> ",
        key_report$HarmonizedTo
      ),
      collapse = "; "
    )
  } else {
    "Key types already compatible."
  }

  status_note <- dplyr::case_when(
    rows_changed ~ "Rows changed after merge. Review for row multiplication or filtering.",
    duplicate_keys ~ "Duplicate key groups detected after merge.",
    cols_wrong ~ "Actual columns added did not match expected columns added.",
    !is.na(match_rate) && match_rate < min_match_rate ~
      paste0(
        "Match rate below threshold: ",
        round(100 * match_rate, 1),
        "% matched; threshold is ",
        round(100 * min_match_rate, 1),
        "%."
      ),
    TRUE ~ "Merge structure looks valid."
  )

  log <- tibble::tibble(
    Merge = name,
    Status = merge_status,
    ReadyForAnalysis = ready_for_analysis,
    RowsBefore = nrow(df_before_join),
    RowsAfter = nrow(df_after),
    ColsBefore = ncol(df_before_join),
    ColsAfter = ncol(df_after),
    ExpectedColsAdded = expected_cols_added,
    ActualColsAdded = actual_cols_added,
    MatchedKeys = summary_row$MatchingKeys,
    LeftUniqueKeys = summary_row$LeftUniqueKeys,
    MatchRate = match_rate,
    DuplicateKeyGroups = duplicate_key_groups,
    UnresolvedDupVars = summary_row$UnresolvedDuplicateVariables,
    KeyHarmonization = key_harmonization_note,
    Note = status_note
  )

  summary_df <- tibble::tibble(
    Metric = c(
      "Status",
      "Rows (before -> after)",
      "Columns (before -> after)",
      "Columns added (expected vs actual)",
      "Keys matched",
      "Match rate",
      "Duplicate key groups",
      "Key harmonization",
      "Note"
    ),
    Value = c(
      status_label,
      paste0(nrow(df_before_join), " -> ", nrow(df_after)),
      paste0(ncol(df_before_join), " -> ", ncol(df_after)),
      paste0("expected +", expected_cols_added, "; actual +", actual_cols_added),
      paste0(summary_row$MatchingKeys, " / ", summary_row$LeftUniqueKeys),
      ifelse(
        is.na(match_rate),
        NA_character_,
        paste0(round(100 * match_rate, 1), "%")
      ),
      as.character(duplicate_key_groups),
      key_harmonization_note,
      status_note
    )
  )

  summary_display <- summary_df

  if (requireNamespace("kableExtra", quietly = TRUE)) {
    summary_display$Value[summary_display$Metric == "Status"] <-
      kableExtra::cell_spec(
        summary_display$Value[summary_display$Metric == "Status"],
        bold = TRUE,
        color = "white",
        background = status_background
      )

    summary_kable <- summary_display %>%
      knitr::kable(
        caption = name,
        escape = FALSE,
        align = c("l", "l")
      ) %>%
      kableExtra::kable_styling(
        full_width = FALSE
      ) %>%
      kableExtra::row_spec(
        1,
        bold = TRUE
      )
  } else {
    summary_kable <- knitr::kable(
      summary_df,
      caption = name,
      align = c("l", "l")
    )
  }

  validation$KeyHarmonization <- key_report
  validation$SummaryTable <- summary_df
  validation$MergeStatus <- merge_status
  validation$MergeReadyForAnalysis <- ready_for_analysis
  validation$MergeNote <- status_note

  list(
    data = df_after,
    validation = validation,
    log = log,
    summary = summary_kable
  )
}
