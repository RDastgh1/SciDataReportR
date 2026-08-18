#' Safely merge two data frames with relationship-aware validation
#'
#' `safe_merge()` performs a merge, validates its structure, logs merge metrics,
#' and returns the merged data plus validation results. It is designed for
#' reproducible database construction pipelines where row count, column count,
#' key coverage, duplicate keys, expected merge relationships, and unresolved
#' duplicate variables need to be audited every time the merge is run.
#'
#' The default `expected_relationship` is `"one-to-one"` to preserve strict
#' historical behavior. For longitudinal databases, use `"many-to-one"` when
#' the left/master data frame has repeated participant IDs because of multiple
#' visits and the right/add-on data frame has one row per participant.
#'
#' Duplicate-variable handling is current-merge aware. If unresolved duplicate
#' variable pairs such as `sy_x` / `sy_y` or `sy.x` / `sy.y` already exist in
#' `df_before`, they are classified as inherited duplicate variables. By default,
#' inherited duplicate variables are reported but do not cause the current merge
#' to fail. New unresolved duplicate variables introduced by the current merge
#' are treated as failures by default.
#'
#' @param df_before Data frame on the left side of the merge.
#' @param df_add Data frame to add to `df_before`.
#' @param by Character vector of merge key column names.
#' @param name Single character label used in the merge log and summary table.
#' @param method Merge method. `"exact"` uses `dplyr::left_join()`.
#'   `"closest_time"` uses `Merge_ByClosestTime()`.
#' @param time_var_before Required when `method = "closest_time"`. Time variable
#'   in `df_before`.
#' @param time_var_add Required when `method = "closest_time"`. Time variable in
#'   `df_add`.
#' @param min_match_rate Numeric between 0 and 1. Merges below this left-key
#'   match rate are marked as `"WARNING"` unless another structural blocker
#'   makes them `"FAIL"`.
#' @param harmonize_keys Logical. If `TRUE`, keys are harmonized using
#'   `HarmonizeMergeKeys()` before merging.
#' @param key_parser Optional parser passed to `HarmonizeMergeKeys()`.
#' @param stop_on_failed_numeric Logical passed to `HarmonizeMergeKeys()`.
#' @param expected_relationship Character. Expected relationship between
#'   `df_before` and `df_add` under `by`. Defaults to `"one-to-one"` to preserve
#'   strict historical behavior. One of:
#'   \itemize{
#'     \item `"one-to-one"`: both sides should be unique by `by`.
#'     \item `"many-to-one"`: left side may repeat keys, right side should be
#'       unique.
#'     \item `"one-to-many"`: left side should be unique, right side may repeat.
#'     \item `"many-to-many"`: both sides may repeat.
#'     \item `"auto"`: infer and report relationship, but do not enforce it.
#'   }
#' @param fail_on_new_duplicate_variables Logical. If `TRUE`, unresolved duplicate
#'   variable pairs introduced by the current merge cause the merge status to be
#'   `"FAIL"`. Default is `TRUE`.
#' @param fail_on_inherited_duplicate_variables Logical. If `TRUE`, unresolved
#'   duplicate variable pairs already present in `df_before` cause the merge
#'   status to be `"FAIL"`. Default is `FALSE`.
#' @param ... Additional arguments passed to `Merge_ByClosestTime()` when
#'   `method = "closest_time"`.
#'
#' @return A list with:
#'   \itemize{
#'     \item `data`: Merged data frame.
#'     \item `validation`: Full validation object from `ValidateMerge()`, with
#'       additional current-merge-aware duplicate-variable diagnostics.
#'     \item `log`: One-row tibble containing merge log metrics.
#'     \item `summary`: A `knitr::kable()` summary table.
#'   }
#'
#' @examples
#' left <- data.frame(
#'   record_id = c(1, 1, 2, 2),
#'   visit_type = c(1, 2, 1, 2),
#'   age = c(40, 40, 55, 55)
#' )
#'
#' right <- data.frame(
#'   record_id = c(1, 2),
#'   imaging_score = c(0.4, 0.8)
#' )
#'
#' m <- safe_merge(
#'   df_before = left,
#'   df_add = right,
#'   by = "record_id",
#'   name = "Example imaging merge",
#'   expected_relationship = "many-to-one"
#' )
#'
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
                       min_match_rate = 0.95,
                       harmonize_keys = TRUE,
                       key_parser = NULL,
                       stop_on_failed_numeric = TRUE,
                       expected_relationship = c(
                         "one-to-one",
                         "many-to-one",
                         "one-to-many",
                         "many-to-many",
                         "auto"
                       ),
                       fail_on_new_duplicate_variables = TRUE,
                       fail_on_inherited_duplicate_variables = FALSE,
                       ...) {

  method <- match.arg(method)
  expected_relationship <- match.arg(expected_relationship)

  ############################################################
  ## Input checks
  ############################################################

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

  if (!is.logical(fail_on_new_duplicate_variables) ||
      length(fail_on_new_duplicate_variables) != 1 ||
      is.na(fail_on_new_duplicate_variables)) {
    stop(
      "fail_on_new_duplicate_variables must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.logical(fail_on_inherited_duplicate_variables) ||
      length(fail_on_inherited_duplicate_variables) != 1 ||
      is.na(fail_on_inherited_duplicate_variables)) {
    stop(
      "fail_on_inherited_duplicate_variables must be TRUE or FALSE.",
      call. = FALSE
    )
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

  ############################################################
  ## Helper: detect unresolved duplicate variable pairs
  ############################################################

  detect_unresolved_duplicate_variables <- function(data) {

    variable_names <- names(data)

    x_dot <- variable_names[grepl("\\.x$", variable_names)]
    y_dot <- variable_names[grepl("\\.y$", variable_names)]

    x_under <- variable_names[grepl("_x$", variable_names)]
    y_under <- variable_names[grepl("_y$", variable_names)]

    dot_pairs <- tibble::tibble(
      Variable = sub("\\.x$", "", x_dot),
      XVariable = x_dot,
      SuffixStyle = ".x/.y"
    ) %>%
      dplyr::mutate(
        YVariable = paste0(Variable, ".y")
      ) %>%
      dplyr::filter(YVariable %in% y_dot)

    under_pairs <- tibble::tibble(
      Variable = sub("_x$", "", x_under),
      XVariable = x_under,
      SuffixStyle = "_x/_y"
    ) %>%
      dplyr::mutate(
        YVariable = paste0(Variable, "_y")
      ) %>%
      dplyr::filter(YVariable %in% y_under)

    dplyr::bind_rows(dot_pairs, under_pairs) %>%
      dplyr::distinct(
        Variable,
        XVariable,
        YVariable,
        SuffixStyle,
        .keep_all = TRUE
      ) %>%
      dplyr::arrange(Variable, XVariable, YVariable)
  }

  ############################################################
  ## Initialize key harmonization report
  ############################################################

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

  ############################################################
  ## Harmonize keys
  ############################################################

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

  ############################################################
  ## Current-merge-aware duplicate-variable baseline
  ############################################################

  unresolved_duplicate_variables_before <- detect_unresolved_duplicate_variables(
    df_before_join
  )

  overlapping_variables_before_merge <- intersect(
    setdiff(names(df_before_join), by),
    setdiff(names(df_add_join), by)
  )

  overlapping_variables_before_merge <- tibble::tibble(
    Variable = overlapping_variables_before_merge
  ) %>%
    dplyr::arrange(Variable)

  ############################################################
  ## Perform merge
  ############################################################

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

  ############################################################
  ## Validate merge
  ############################################################

  validation <- ValidateMerge(
    LeftData = df_before_join,
    RightData = df_add_join,
    MergedData = df_after,
    keys = by,
    expected_relationship = expected_relationship
  )

  summary_row <- validation$Summary

  ############################################################
  ## Current-merge-aware duplicate-variable classification
  ############################################################

  unresolved_duplicate_variables_after <- detect_unresolved_duplicate_variables(
    df_after
  )

  inherited_duplicate_variables <- unresolved_duplicate_variables_after %>%
    dplyr::semi_join(
      unresolved_duplicate_variables_before,
      by = c("Variable", "XVariable", "YVariable", "SuffixStyle")
    )

  new_duplicate_variables <- unresolved_duplicate_variables_after %>%
    dplyr::anti_join(
      unresolved_duplicate_variables_before,
      by = c("Variable", "XVariable", "YVariable", "SuffixStyle")
    )

  n_unresolved_duplicate_variables_before <- nrow(
    unresolved_duplicate_variables_before
  )

  n_unresolved_duplicate_variables_after <- nrow(
    unresolved_duplicate_variables_after
  )

  n_inherited_duplicate_variables <- nrow(inherited_duplicate_variables)
  n_new_duplicate_variables <- nrow(new_duplicate_variables)

  ############################################################
  ## Core merge metrics
  ############################################################

  expected_cols_added <- ncol(df_add_join) - length(by)
  actual_cols_added <- ncol(df_after) - ncol(df_before_join)

  rows_changed <- nrow(df_after) != nrow(df_before_join)

  row_change_allowed <- expected_relationship %in% c(
    "one-to-many",
    "many-to-many",
    "auto"
  )

  rows_changed_blocker <- rows_changed && !row_change_allowed

  cols_wrong <- actual_cols_added != expected_cols_added

  duplicate_key_blockers <- summary_row$DuplicateKeyBlockers

  relationship_matches_expected <- summary_row$RelationshipMatchesExpected

  match_rate <- ifelse(
    summary_row$LeftUniqueKeys > 0,
    summary_row$MatchingKeys / summary_row$LeftUniqueKeys,
    NA_real_
  )

  new_duplicate_blocker <- isTRUE(fail_on_new_duplicate_variables) &&
    n_new_duplicate_variables > 0

  inherited_duplicate_blocker <- isTRUE(fail_on_inherited_duplicate_variables) &&
    n_inherited_duplicate_variables > 0

  ############################################################
  ## Merge status
  ############################################################

  merge_status <- dplyr::case_when(
    rows_changed_blocker ~ "FAIL",
    cols_wrong ~ "FAIL",
    duplicate_key_blockers > 0 ~ "FAIL",
    expected_relationship != "auto" && !relationship_matches_expected ~ "FAIL",
    new_duplicate_blocker ~ "FAIL",
    inherited_duplicate_blocker ~ "FAIL",
    !is.na(match_rate) && match_rate < min_match_rate ~ "WARNING",
    n_inherited_duplicate_variables > 0 ~ "WARNING",
    TRUE ~ "PASS"
  )

  ready_for_analysis <- merge_status != "FAIL"

  ############################################################
  ## Notes
  ############################################################

  key_harmonization_note <- if (nrow(key_report) == 0) {
    "No key harmonization performed."
  } else if (
    any(key_report$LeftTypeBefore != key_report$RightTypeBefore) ||
      any(key_report$ParserUsed)
  ) {
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

  duplicate_variable_note <- dplyr::case_when(
    n_new_duplicate_variables > 0 && n_inherited_duplicate_variables > 0 ~
      paste0(
        "Merged data contains ",
        n_new_duplicate_variables,
        " new unresolved duplicate variable pair(s) and ",
        n_inherited_duplicate_variables,
        " inherited unresolved duplicate variable pair(s)."
      ),
    n_new_duplicate_variables > 0 ~
      paste0(
        "Merged data contains ",
        n_new_duplicate_variables,
        " new unresolved duplicate variable pair(s) introduced by this merge."
      ),
    n_inherited_duplicate_variables > 0 ~
      paste0(
        "Merged data contains ",
        n_inherited_duplicate_variables,
        " inherited unresolved duplicate variable pair(s) already present in df_before."
      ),
    TRUE ~
      "No unresolved duplicate variable pairs detected."
  )

  status_note <- dplyr::case_when(
    rows_changed_blocker ~
      "Rows changed after merge, but expected_relationship should preserve left-side row count.",
    cols_wrong ~
      "Actual columns added did not match expected columns added.",
    duplicate_key_blockers > 0 ~
      paste0(
        "Duplicate key groups violate expected_relationship = '",
        expected_relationship,
        "'."
      ),
    expected_relationship != "auto" && !relationship_matches_expected ~
      paste0(
        "Detected relationship '",
        summary_row$DetectedRelationship,
        "' does not match expected_relationship = '",
        expected_relationship,
        "'."
      ),
    new_duplicate_blocker ~
      paste0(
        "This merge introduced unresolved duplicate variable pair(s): ",
        paste(new_duplicate_variables$Variable, collapse = ", "),
        ". Resolve or rename overlapping variables before merging."
      ),
    inherited_duplicate_blocker ~
      paste0(
        "df_before already contained unresolved duplicate variable pair(s): ",
        paste(inherited_duplicate_variables$Variable, collapse = ", "),
        "."
      ),
    !is.na(match_rate) && match_rate < min_match_rate ~
      paste0(
        "Match rate below threshold: ",
        round(100 * match_rate, 1),
        "% matched; threshold is ",
        round(100 * min_match_rate, 1),
        "%."
      ),
    n_inherited_duplicate_variables > 0 ~
      paste0(
        "Merge structure looks valid, but inherited unresolved duplicate variable pair(s) remain: ",
        paste(inherited_duplicate_variables$Variable, collapse = ", "),
        "."
      ),
    TRUE ~
      "Merge structure looks valid."
  )

  ############################################################
  ## Log
  ############################################################

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
    ExpectedRelationship = expected_relationship,
    DetectedRelationship = summary_row$DetectedRelationship,
    RelationshipMatchesExpected = relationship_matches_expected,
    MatchedKeys = summary_row$MatchingKeys,
    LeftUniqueKeys = summary_row$LeftUniqueKeys,
    MatchRate = match_rate,
    DuplicateKeyBlockers = duplicate_key_blockers,
    DuplicateKeyGroups_Left = summary_row$DuplicateKeyGroups_Left,
    DuplicateKeyGroups_Right = summary_row$DuplicateKeyGroups_Right,
    DuplicateKeyGroups_Merged = summary_row$DuplicateKeyGroups_Merged,
    DuplicateKeyGroups = duplicate_key_blockers,
    UnresolvedDupVars = n_unresolved_duplicate_variables_after,
    NewUnresolvedDupVars = n_new_duplicate_variables,
    InheritedUnresolvedDupVars = n_inherited_duplicate_variables,
    PreExistingUnresolvedDupVars = n_unresolved_duplicate_variables_before,
    KeyHarmonization = key_harmonization_note,
    DuplicateVariableNote = duplicate_variable_note,
    Note = status_note
  )

  ############################################################
  ## Summary table
  ############################################################

  summary_df <- tibble::tibble(
    Metric = c(
      "Status",
      "Rows (before -> after)",
      "Columns (before -> after)",
      "Columns added (expected vs actual)",
      "Expected relationship",
      "Detected relationship",
      "Relationship matches expected",
      "Keys matched",
      "Match rate",
      "Duplicate key blockers",
      "Duplicate key groups left/right/merged",
      "Unresolved duplicate variables before merge",
      "Unresolved duplicate variables after merge",
      "New unresolved duplicate variables",
      "Inherited unresolved duplicate variables",
      "Overlapping variables before merge",
      "Key harmonization",
      "Duplicate variable note",
      "Note"
    ),
    Value = c(
      merge_status,
      paste0(nrow(df_before_join), " -> ", nrow(df_after)),
      paste0(ncol(df_before_join), " -> ", ncol(df_after)),
      paste0("expected +", expected_cols_added, "; actual +", actual_cols_added),
      expected_relationship,
      summary_row$DetectedRelationship,
      as.character(relationship_matches_expected),
      paste0(summary_row$MatchingKeys, " / ", summary_row$LeftUniqueKeys),
      ifelse(
        is.na(match_rate),
        NA_character_,
        paste0(round(100 * match_rate, 1), "%")
      ),
      as.character(duplicate_key_blockers),
      paste0(
        summary_row$DuplicateKeyGroups_Left,
        " / ",
        summary_row$DuplicateKeyGroups_Right,
        " / ",
        summary_row$DuplicateKeyGroups_Merged
      ),
      as.character(n_unresolved_duplicate_variables_before),
      as.character(n_unresolved_duplicate_variables_after),
      as.character(n_new_duplicate_variables),
      as.character(n_inherited_duplicate_variables),
      as.character(nrow(overlapping_variables_before_merge)),
      key_harmonization_note,
      duplicate_variable_note,
      status_note
    )
  )

  summary_kable <- summary_df %>%
    knitr::kable(
      caption = name,
      escape = TRUE,
      align = c("l", "l")
    )

  if (requireNamespace("kableExtra", quietly = TRUE)) {
    summary_kable <- summary_kable %>%
      kableExtra::kable_styling(full_width = FALSE) %>%
      kableExtra::row_spec(1, bold = TRUE)
  }

  ############################################################
  ## Add enhanced diagnostics to validation object
  ############################################################

  validation$KeyHarmonization <- key_report
  validation$SummaryTable <- summary_df
  validation$MergeStatus <- merge_status
  validation$MergeReadyForAnalysis <- ready_for_analysis
  validation$MergeNote <- status_note

  validation$DuplicateVariableAudit <- list(
    Before = unresolved_duplicate_variables_before,
    After = unresolved_duplicate_variables_after,
    New = new_duplicate_variables,
    Inherited = inherited_duplicate_variables,
    OverlappingVariablesBeforeMerge = overlapping_variables_before_merge
  )

  validation$DuplicateVariableSummary <- tibble::tibble(
    Metric = c(
      "Before",
      "After",
      "New",
      "Inherited",
      "Overlapping variables before merge"
    ),
    Count = c(
      n_unresolved_duplicate_variables_before,
      n_unresolved_duplicate_variables_after,
      n_new_duplicate_variables,
      n_inherited_duplicate_variables,
      nrow(overlapping_variables_before_merge)
    )
  )

  ############################################################
  ## Return
  ############################################################

  list(
    data = df_after,
    validation = validation,
    log = log,
    summary = summary_kable
  )
}
