#' Merge two data frames and validate the result in one step
#'
#' Perform a left join (exact-key or closest-time) and immediately audit the
#' result with [ValidateMerge()]. The returned object bundles the merged data,
#' the full validation object, a one-row log tibble suitable for stacking
#' across a pipeline with [merge_summary_table()], and a compact
#' `knitr::kable()` summary. Nothing is printed automatically.
#'
#' @param df_before A data frame. The "left" dataset that rows are added to.
#' @param df_add A data frame. The "right" dataset being merged in.
#' @param by Character vector of merge keys. Multiple keys are supported,
#'   such as `c("subject_id", "visit")` for longitudinal data.
#' @param name Required. A short character label for this merge, used in the
#'   log tibble and as the caption of the returned summary table.
#' @param method Either `"exact"` (default) or `"closest_time"`.
#'   `"exact"` performs `dplyr::left_join(df_before, df_add, by = by)`.
#'   `"closest_time"` calls [Merge_ByClosestTime()] using `by` as the
#'   exact-match key and matches each `df_before` row to the `df_add` row with
#'   the nearest time.
#' @param time_var_before For `method = "closest_time"`, the name of the time
#'   variable in `df_before` (as a string).
#' @param time_var_add For `method = "closest_time"`, the name of the time
#'   variable in `df_add` (as a string).
#' @param ... Additional arguments passed to [Merge_ByClosestTime()] when
#'   `method = "closest_time"` (for example `is_date = TRUE`). Ignored for
#'   `method = "exact"`.
#'
#' @details
#' For `method = "closest_time"`, [ValidateMerge()] is still run with `by` as
#' the keys. In longitudinal data where `df_before` legitimately contains
#' repeated visits per key, the Duplicate Keys check will report those repeats
#' as failures even though they are expected. The validation output is left
#' untouched (it is not suppressed or reinterpreted); review the duplicate-key
#' results with that caveat in mind.
#'
#' @return A list with four elements (nothing is auto-printed):
#' \describe{
#'   \item{data}{The merged data frame.}
#'   \item{validation}{The full [ValidateMerge()] result.}
#'   \item{log}{A one-row tibble with columns `Merge`, `Status` (worst of
#'     `validation$Checks$Status`, ordered FAIL > WARNING > PASS),
#'     `ReadyForAnalysis`, `RowsBefore`, `RowsAfter`, `ColsBefore`,
#'     `ColsAfter`, `MatchedKeys`, `DuplicateKeyGroups`, and
#'     `UnresolvedDupVars`.}
#'   \item{summary}{A `knitr::kable()` table (returned, not printed) in
#'     metric/value layout captioned with `name`.}
#' }
#'
#' @seealso [ValidateMerge()], [merge_summary_table()], [merge_detail()],
#'   [ExploreMergeValidation()]
#'
#' @examples
#' baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
#' labs <- data.frame(id = c(1, 2, 4), glucose = c(90, 110, 100))
#' m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")
#' m$log
#'
#' @export
safe_merge <- function(df_before,
                       df_add,
                       by,
                       name,
                       method = c("exact", "closest_time"),
                       time_var_before = NULL,
                       time_var_add = NULL,
                       ...) {

  method <- match.arg(method)

  # Validate inputs

  if (!is.data.frame(df_before)) {
    stop("df_before must be a data.frame.")
  }

  if (!is.data.frame(df_add)) {
    stop("df_add must be a data.frame.")
  }

  if (missing(by) || !is.character(by) || length(by) == 0) {
    stop("by must be a character vector of one or more merge keys.")
  }

  if (missing(name)) {
    stop("name is required. Provide a short label for this merge, e.g. name = \"Baseline + labs\".")
  }

  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("name must be a single non-empty character value.")
  }

  # Perform the merge

  if (method == "exact") {

    df_after <- dplyr::left_join(df_before, df_add, by = by)

  } else {

    if (is.null(time_var_before) || is.null(time_var_add)) {
      stop(
        "method = \"closest_time\" requires both time_var_before and time_var_add."
      )
    }

    if (!time_var_before %in% names(df_before)) {
      stop("time_var_before ('", time_var_before, "') was not found in df_before.")
    }

    if (!time_var_add %in% names(df_add)) {
      stop("time_var_add ('", time_var_add, "') was not found in df_add.")
    }

    closest_result <- Merge_ByClosestTime(
      DataFrame1 = df_before,
      DataFrame2 = df_add,
      TimeVar1 = time_var_before,
      TimeVar2 = time_var_add,
      keys = by,
      ...
    )

    df_after <- closest_result$merged_dataframe
  }

  # Validate the merge.
  #
  # Caveat for method = "closest_time": ValidateMerge() audits duplicate key
  # combinations using `by` only (not the time variables), so legitimate
  # repeated visits per key in longitudinal data are reported as duplicate-key
  # failures. The validation output is intentionally left as-is rather than
  # suppressed or reinterpreted; the caveat is documented in the roxygen
  # Details section.
  validation <- ValidateMerge(
    LeftData = df_before,
    RightData = df_add,
    MergedData = df_after,
    keys = by
  )

  # Build the one-row log

  status_levels <- c("PASS", "WARNING", "FAIL")

  worst_status <- status_levels[
    max(match(validation$Checks$Status, status_levels), na.rm = TRUE)
  ]

  summary_row <- validation$Summary

  duplicate_key_groups <- summary_row$DuplicateKeyGroups_Left +
    summary_row$DuplicateKeyGroups_Right +
    summary_row$DuplicateKeyGroups_Merged

  log <- tibble::tibble(
    Merge = name,
    Status = worst_status,
    ReadyForAnalysis = validation$ReadyForAnalysis,
    RowsBefore = nrow(df_before),
    RowsAfter = nrow(df_after),
    ColsBefore = ncol(df_before),
    ColsAfter = ncol(df_after),
    MatchedKeys = summary_row$MatchingKeys,
    DuplicateKeyGroups = duplicate_key_groups,
    UnresolvedDupVars = summary_row$UnresolvedDuplicateVariables
  )

  # Build the kable summary (returned, never printed here)

  expected_cols_added <- ncol(df_add) - length(by)
  actual_cols_added <- ncol(df_after) - ncol(df_before)

  summary_df <- data.frame(
    Metric = c(
      "Rows (before → after)",
      "Columns (before → after)",
      "Columns added (expected vs actual)",
      "Keys matched",
      "Duplicate key groups",
      "Status"
    ),
    Value = c(
      paste0(nrow(df_before), " → ", nrow(df_after)),
      paste0(ncol(df_before), " → ", ncol(df_after)),
      paste0("expected +", expected_cols_added, "; actual +", actual_cols_added),
      paste0(summary_row$MatchingKeys, " / ", summary_row$LeftUniqueKeys),
      as.character(duplicate_key_groups),
      worst_status
    ),
    stringsAsFactors = FALSE
  )

  summary_kable <- knitr::kable(
    summary_df,
    caption = name
  )

  list(
    data = df_after,
    validation = validation,
    log = log,
    summary = summary_kable
  )
}
