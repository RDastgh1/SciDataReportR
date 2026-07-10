#' Validate a merge between two source data frames and a merged result
#'
#' `ValidateMerge()` audits a completed merge by comparing the left/source data,
#' right/add-on data, and merged result. It checks key coverage, key uniqueness,
#' expected merge relationship, row inflation, overlapping non-key variables,
#' unresolved `.x`/`.y` or `_x`/`_y` columns, and duplicated-variable conflicts.
#'
#' This version supports expected merge relationships. The default is
#' `"one-to-one"` to preserve strict historical behavior. For longitudinal
#' databases, use `"many-to-one"` when the left/master data frame is a
#' participant-visit table and the right/add-on data frame is participant-level.
#'
#' @param LeftData A data frame used as the left side of the merge.
#' @param RightData A data frame used as the right side of the merge.
#' @param MergedData A data frame produced by merging `LeftData` and `RightData`.
#' @param keys Character vector of merge key column names.
#' @param Keys Deprecated. Use `keys`.
#' @param expected_relationship Character. Expected relationship between
#'   `LeftData` and `RightData` under `keys`. Defaults to `"one-to-one"` to
#'   preserve strict historical behavior. One of:
#'   \itemize{
#'     \item `"one-to-one"`: left keys and right keys should both be unique.
#'     \item `"many-to-one"`: left keys may repeat, right keys should be unique.
#'       This is common when merging participant-level data into a longitudinal
#'       participant-visit master table by `record_id`.
#'     \item `"one-to-many"`: left keys should be unique, right keys may repeat.
#'     \item `"many-to-many"`: both sides may repeat. This is allowed only when
#'       explicitly expected.
#'     \item `"auto"`: infer and report relationship, but do not fail solely
#'       because of duplicate keys or relationship type.
#'   }
#'
#' @return A list containing merge validation summaries, checks, relationship
#'   audits, duplicate-key audits, coverage tables, overlap audits, duplicated
#'   variable audits, conflict tables, and suggested actions.
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
#' merged <- dplyr::left_join(left, right, by = "record_id")
#'
#' validation <- ValidateMerge(
#'   LeftData = left,
#'   RightData = right,
#'   MergedData = merged,
#'   keys = "record_id",
#'   expected_relationship = "many-to-one"
#' )
#'
#' validation$Summary
#' @export
ValidateMerge <- function(LeftData,
                          RightData,
                          MergedData,
                          keys,
                          Keys = lifecycle::deprecated(),
                          expected_relationship = c(
                            "one-to-one",
                            "many-to-one",
                            "one-to-many",
                            "many-to-many",
                            "auto"
                          )) {

  expected_relationship <- match.arg(expected_relationship)

  if (lifecycle::is_present(Keys)) {
    lifecycle::deprecate_warn(
      "19.15.0",
      "ValidateMerge(Keys)",
      "ValidateMerge(keys)"
    )
    keys <- Keys
  }

  if (!is.data.frame(LeftData)) {
    stop("LeftData must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(RightData)) {
    stop("RightData must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(MergedData)) {
    stop("MergedData must be a data.frame.", call. = FALSE)
  }

  if (missing(keys) || length(keys) == 0) {
    stop("keys must be supplied as a character vector.", call. = FALSE)
  }

  if (!is.character(keys)) {
    stop("keys must be a character vector.", call. = FALSE)
  }

  Keys <- keys

  missing_left <- setdiff(Keys, names(LeftData))
  missing_right <- setdiff(Keys, names(RightData))
  missing_merged <- setdiff(Keys, names(MergedData))

  if (length(missing_left) > 0) {
    stop(
      "The following key variable(s) are missing from LeftData: ",
      paste(missing_left, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_right) > 0) {
    stop(
      "The following key variable(s) are missing from RightData: ",
      paste(missing_right, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_merged) > 0) {
    missing_merged_audit <- purrr::map_chr(
      missing_merged,
      function(key) {
        x_name_dot <- paste0(key, ".x")
        y_name_dot <- paste0(key, ".y")
        x_name_under <- paste0(key, "_x")
        y_name_under <- paste0(key, "_y")

        if (x_name_dot %in% names(MergedData) &&
            y_name_dot %in% names(MergedData)) {
          paste0(
            key,
            " appears to be present as ",
            x_name_dot,
            " and ",
            y_name_dot,
            ". This often happens when the merge was performed using fewer keys than intended."
          )
        } else if (x_name_under %in% names(MergedData) &&
                   y_name_under %in% names(MergedData)) {
          paste0(
            key,
            " appears to be present as ",
            x_name_under,
            " and ",
            y_name_under,
            ". This often happens when the merge was performed using fewer keys than intended."
          )
        } else {
          paste0(key, " was not found in MergedData.")
        }
      }
    )

    stop(
      "The following intended key variable(s) are missing from MergedData: ",
      paste(missing_merged, collapse = ", "),
      "\n\n",
      paste(missing_merged_audit, collapse = "\n"),
      "\n\nCheck whether the merge was performed with the intended `by = ` variables.",
      call. = FALSE
    )
  }

  empty_duplicate_variables <- function() {
    tibble::tibble(
      Variable = character(),
      XVariable = character(),
      YVariable = character(),
      LeftClass = character(),
      RightClass = character(),
      Agreement = numeric(),
      Conflicts = integer(),
      MissingnessConflicts = integer(),
      BothMissing = integer(),
      TotalRows = integer()
    )
  }

  empty_variable_conflicts <- function(keys) {
    key_cols <- stats::setNames(rep(list(character()), length(keys)), keys)

    tibble::tibble(
      !!!key_cols,
      Variable = character(),
      LeftValue = character(),
      RightValue = character(),
      ConflictType = character()
    )
  }

  empty_key_table <- function(keys) {
    tibble::as_tibble(stats::setNames(rep(list(character()), length(keys)), keys))
  }

  key_complete <- function(df, keys) {
    df %>%
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(keys),
          ~ !is.na(.) & trimws(as.character(.)) != ""
        )
      )
  }

  count_missing_key_rows <- function(df, keys) {
    df %>%
      dplyr::filter(
        dplyr::if_any(
          dplyr::all_of(keys),
          ~ is.na(.) | trimws(as.character(.)) == ""
        )
      ) %>%
      nrow()
  }

  distinct_keys <- function(df, keys) {
    if (nrow(df) == 0) {
      return(empty_key_table(keys))
    }

    df %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(keys)))
  }

  duplicate_key_groups_table <- function(df, keys) {
    if (nrow(df) == 0) {
      out <- empty_key_table(keys)
      out$.n <- integer()
      return(out)
    }

    df %>%
      dplyr::count(dplyr::across(dplyr::all_of(keys)), name = ".n") %>%
      dplyr::filter(.n > 1)
  }

  key_types <- tibble::tibble(
    Key = Keys,
    LeftType = purrr::map_chr(
      Keys,
      ~ paste(class(LeftData[[.x]]), collapse = "/")
    ),
    RightType = purrr::map_chr(
      Keys,
      ~ paste(class(RightData[[.x]]), collapse = "/")
    ),
    MergedType = purrr::map_chr(
      Keys,
      ~ paste(class(MergedData[[.x]]), collapse = "/")
    )
  )

  key_type_count <- key_types %>%
    dplyr::filter(
      LeftType != RightType |
        LeftType != MergedType |
        RightType != MergedType
    ) %>%
    nrow()

  for (key in Keys) {
    LeftData[[key]] <- as.character(LeftData[[key]])
    RightData[[key]] <- as.character(RightData[[key]])
    MergedData[[key]] <- as.character(MergedData[[key]])
  }

  missing_key_rows_left <- count_missing_key_rows(LeftData, Keys)
  missing_key_rows_right <- count_missing_key_rows(RightData, Keys)
  missing_key_rows_merged <- count_missing_key_rows(MergedData, Keys)

  LeftData_keyed <- key_complete(LeftData, Keys)
  RightData_keyed <- key_complete(RightData, Keys)
  MergedData_keyed <- key_complete(MergedData, Keys)

  left_keys <- distinct_keys(LeftData_keyed, Keys)
  right_keys <- distinct_keys(RightData_keyed, Keys)
  merged_keys <- distinct_keys(MergedData_keyed, Keys)

  matching_keys <- dplyr::inner_join(left_keys, right_keys, by = Keys)
  left_only_keys <- dplyr::anti_join(left_keys, right_keys, by = Keys)
  right_only_keys <- dplyr::anti_join(right_keys, left_keys, by = Keys)

  id_coverage <- list(
    Matching = matching_keys,
    LeftOnly = left_only_keys,
    RightOnly = right_only_keys
  )

  fingerprint <- tibble::tibble(
    Metric = c(
      "Rows",
      "Columns",
      "Unique Key Combinations",
      "Rows With Missing Keys"
    ),
    Left = c(
      nrow(LeftData),
      ncol(LeftData),
      nrow(left_keys),
      missing_key_rows_left
    ),
    Right = c(
      nrow(RightData),
      ncol(RightData),
      nrow(right_keys),
      missing_key_rows_right
    ),
    Merged = c(
      nrow(MergedData),
      ncol(MergedData),
      nrow(merged_keys),
      missing_key_rows_merged
    )
  )

  duplicate_keys_left <- duplicate_key_groups_table(LeftData_keyed, Keys)
  duplicate_keys_right <- duplicate_key_groups_table(RightData_keyed, Keys)
  duplicate_keys_merged <- duplicate_key_groups_table(MergedData_keyed, Keys)

  duplicate_key_groups_left <- nrow(duplicate_keys_left)
  duplicate_key_groups_right <- nrow(duplicate_keys_right)
  duplicate_key_groups_merged <- nrow(duplicate_keys_merged)

  duplicate_keys <- list(
    Left = duplicate_keys_left,
    Right = duplicate_keys_right,
    Merged = duplicate_keys_merged
  )

  left_unique <- duplicate_key_groups_left == 0
  right_unique <- duplicate_key_groups_right == 0
  merged_unique <- duplicate_key_groups_merged == 0

  detected_relationship <- dplyr::case_when(
    left_unique & right_unique ~ "one-to-one",
    !left_unique & right_unique ~ "many-to-one",
    left_unique & !right_unique ~ "one-to-many",
    TRUE ~ "many-to-many"
  )

  relationship_matches_expected <- if (expected_relationship == "auto") {
    TRUE
  } else {
    identical(detected_relationship, expected_relationship)
  }

  relationship_duplicate_blockers <- dplyr::case_when(
    expected_relationship == "auto" ~ 0L,

    expected_relationship == "one-to-one" ~
      duplicate_key_groups_left +
      duplicate_key_groups_right +
      duplicate_key_groups_merged,

    expected_relationship == "many-to-one" ~
      duplicate_key_groups_right,

    expected_relationship == "one-to-many" ~
      duplicate_key_groups_left,

    expected_relationship == "many-to-many" ~ 0L,

    TRUE ~
      duplicate_key_groups_left +
      duplicate_key_groups_right +
      duplicate_key_groups_merged
  )

  relationship_status <- dplyr::case_when(
    expected_relationship == "auto" ~ "PASS",
    relationship_matches_expected ~ "PASS",
    relationship_duplicate_blockers > 0 ~ "FAIL",
    TRUE ~ "FAIL"
  )

  relationship_details <- dplyr::case_when(
    expected_relationship == "auto" ~ paste0(
      "Detected relationship is ",
      detected_relationship,
      ". No expected relationship was enforced."
    ),
    relationship_matches_expected ~ paste0(
      "Detected relationship matches expected_relationship = '",
      expected_relationship,
      "'."
    ),
    relationship_duplicate_blockers > 0 ~ paste0(
      "Detected relationship is ",
      detected_relationship,
      ", but expected_relationship = '",
      expected_relationship,
      "'. Duplicate key blockers: ",
      relationship_duplicate_blockers,
      "."
    ),
    TRUE ~ paste0(
      "Detected relationship is ",
      detected_relationship,
      ", but expected_relationship = '",
      expected_relationship,
      "'."
    )
  )

  relationship <- tibble::tibble(
    Dataset = c("Left", "Right", "Merged"),
    UniqueKeys = c(left_unique, right_unique, merged_unique),
    DuplicateKeyGroups = c(
      duplicate_key_groups_left,
      duplicate_key_groups_right,
      duplicate_key_groups_merged
    ),
    ExpectedRelationship = expected_relationship,
    DetectedRelationship = detected_relationship,
    RelationshipMatchesExpected = relationship_matches_expected,
    RelationshipDuplicateBlockers = relationship_duplicate_blockers,
    RelationshipStatus = relationship_status
  )

  source_overlap_vars <- intersect(names(LeftData), names(RightData))
  overlapping_non_key_vars <- setdiff(source_overlap_vars, Keys)

  overlapping_variables <- tibble::tibble(
    Variable = sort(overlapping_non_key_vars)
  )

  potential_merge_risk <- tibble::tibble(
    Variable = sort(overlapping_non_key_vars),
    Risk = "Present in both source datasets but not specified as a merge key."
  )

  join_audit <- tibble::tibble(
    Variable = sort(source_overlap_vars)
  ) %>%
    dplyr::mutate(
      InBoth = TRUE,
      IsKey = Variable %in% Keys,
      JoinRole = dplyr::case_when(
        IsKey ~ "Specified Key",
        TRUE ~ "Overlap Not In Keys"
      )
    )

  all_source_vars <- sort(unique(c(names(LeftData), names(RightData))))

  overlap_audit <- tibble::tibble(
    Variable = all_source_vars
  ) %>%
    dplyr::mutate(
      InLeft = Variable %in% names(LeftData),
      InRight = Variable %in% names(RightData),
      IsKey = Variable %in% Keys,
      Risk = dplyr::case_when(
        InLeft & InRight & !IsKey ~ "Potential merge risk",
        InLeft & InRight & IsKey ~ "Expected key overlap",
        TRUE ~ "None"
      )
    )

  detect_duplicate_pairs <- function(data_names, suffix_x, suffix_y) {
    x_vars <- grep(paste0("\\", suffix_x, "$"), data_names, value = TRUE)

    if (length(x_vars) == 0) {
      return(
        tibble::tibble(
          Variable = character(),
          XVar = character(),
          YVar = character()
        )
      )
    }

    purrr::map_dfr(
      x_vars,
      function(x_var) {
        base_var <- sub(paste0("\\", suffix_x, "$"), "", x_var)
        y_var <- paste0(base_var, suffix_y)

        if (y_var %in% data_names) {
          tibble::tibble(
            Variable = base_var,
            XVar = x_var,
            YVar = y_var
          )
        } else {
          NULL
        }
      }
    )
  }

  duplicate_pairs_dot <- detect_duplicate_pairs(names(MergedData), ".x", ".y")
  duplicate_pairs_under <- detect_duplicate_pairs(names(MergedData), "_x", "_y")

  duplicate_pairs <- dplyr::bind_rows(
    duplicate_pairs_dot,
    duplicate_pairs_under
  ) %>%
    dplyr::distinct(Variable, XVar, YVar) %>%
    dplyr::arrange(Variable)

  unresolved_duplicate_variables <- duplicate_pairs %>%
    dplyr::select(Variable) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Variable)

  if (nrow(unresolved_duplicate_variables) == 0) {
    unresolved_duplicate_variables <- tibble::tibble(Variable = character())
  }

  empty_duplicate_variables <- function() {
    tibble::tibble(
      Variable = character(),
      XVariable = character(),
      YVariable = character(),
      LeftClass = character(),
      RightClass = character(),
      Agreement = numeric(),
      Conflicts = integer(),
      MissingnessConflicts = integer(),
      BothMissing = integer(),
      TotalRows = integer()
    )
  }

  if (nrow(duplicate_pairs) > 0) {
    duplicate_variables <- purrr::map_dfr(
      seq_len(nrow(duplicate_pairs)),
      function(i) {
        x_var <- duplicate_pairs$XVar[i]
        y_var <- duplicate_pairs$YVar[i]

        x <- as.character(MergedData[[x_var]])
        y <- as.character(MergedData[[y_var]])

        agree <- dplyr::case_when(
          is.na(x) & is.na(y) ~ TRUE,
          is.na(x) | is.na(y) ~ FALSE,
          TRUE ~ x == y
        )

        tibble::tibble(
          Variable = duplicate_pairs$Variable[i],
          XVariable = x_var,
          YVariable = y_var,
          LeftClass = paste(class(MergedData[[x_var]]), collapse = "/"),
          RightClass = paste(class(MergedData[[y_var]]), collapse = "/"),
          Agreement = round(mean(agree) * 100, 2),
          Conflicts = sum(!agree),
          MissingnessConflicts = sum(is.na(x) != is.na(y)),
          BothMissing = sum(is.na(x) & is.na(y)),
          TotalRows = length(agree)
        )
      }
    ) %>%
      dplyr::arrange(dplyr::desc(Conflicts), Variable)
  } else {
    duplicate_variables <- empty_duplicate_variables()
  }

  if (nrow(duplicate_variables) > 0) {
    suspicious_conflicts <- duplicate_variables %>%
      dplyr::filter(Agreement < 75 | LeftClass != RightClass) %>%
      dplyr::arrange(Variable)
  } else {
    suspicious_conflicts <- empty_duplicate_variables()
  }

  if (nrow(duplicate_pairs) > 0) {
    variable_conflicts <- purrr::map_dfr(
      seq_len(nrow(duplicate_pairs)),
      function(i) {
        x_var <- duplicate_pairs$XVar[i]
        y_var <- duplicate_pairs$YVar[i]

        tmp <- MergedData %>%
          dplyr::select(dplyr::all_of(Keys), dplyr::all_of(c(x_var, y_var)))

        x <- as.character(tmp[[x_var]])
        y <- as.character(tmp[[y_var]])

        agree <- dplyr::case_when(
          is.na(x) & is.na(y) ~ TRUE,
          is.na(x) | is.na(y) ~ FALSE,
          TRUE ~ x == y
        )

        idx <- !agree

        if (!any(idx)) {
          return(NULL)
        }

        tmp[idx, Keys, drop = FALSE] %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            Variable = duplicate_pairs$Variable[i],
            LeftValue = x[idx],
            RightValue = y[idx],
            ConflictType = dplyr::case_when(
              is.na(x[idx]) & !is.na(y[idx]) ~ "Left missing, Right present",
              !is.na(x[idx]) & is.na(y[idx]) ~ "Left present, Right missing",
              TRUE ~ "Different non-missing values"
            )
          )
      }
    )

    if (nrow(variable_conflicts) == 0) {
      variable_conflicts <- empty_variable_conflicts(Keys)
    } else {
      variable_conflicts <- variable_conflicts %>%
        dplyr::arrange(Variable)
    }
  } else {
    variable_conflicts <- empty_variable_conflicts(Keys)
  }

  max_source_rows <- max(nrow(LeftData), nrow(RightData))

  row_inflation_factor <- ifelse(
    max_source_rows > 0,
    round(nrow(MergedData) / max_source_rows, 3),
    NA_real_
  )

  row_count_expected_to_change <- expected_relationship %in% c(
    "one-to-many",
    "many-to-many",
    "auto"
  )

  row_count_changed <- nrow(MergedData) != nrow(LeftData)
  row_count_blocker <- row_count_changed && !row_count_expected_to_change

  coverage_total <- nrow(left_only_keys) + nrow(right_only_keys)
  overlap_total <- nrow(overlapping_variables)
  unresolved_duplicate_total <- nrow(unresolved_duplicate_variables)
  conflict_total <- nrow(variable_conflicts)
  suspicious_conflict_total <- nrow(suspicious_conflicts)
  missing_key_total <- missing_key_rows_left +
    missing_key_rows_right +
    missing_key_rows_merged

  checks <- tibble::tibble(
    Check = c(
      "Key Types",
      "Missing Keys",
      "Expected Relationship",
      "Duplicate Keys",
      "Coverage",
      "Row Count",
      "Row Inflation",
      "Overlapping Variables",
      "Unresolved Duplicate Variables",
      "Variable Conflicts",
      "Suspicious Conflicts"
    ),
    Count = c(
      key_type_count,
      missing_key_total,
      relationship_duplicate_blockers,
      relationship_duplicate_blockers,
      coverage_total,
      as.integer(row_count_changed),
      row_inflation_factor,
      overlap_total,
      unresolved_duplicate_total,
      conflict_total,
      suspicious_conflict_total
    )
  ) %>%
    dplyr::mutate(
      Status = dplyr::case_when(
        Check == "Key Types" & Count > 0 ~ "WARNING",
        Check == "Missing Keys" & Count > 0 ~ "WARNING",
        Check == "Expected Relationship" &
          expected_relationship != "auto" &
          !relationship_matches_expected ~ "FAIL",
        Check == "Expected Relationship" ~ "PASS",
        Check == "Duplicate Keys" & relationship_duplicate_blockers > 0 ~ "FAIL",
        Check == "Duplicate Keys" ~ "PASS",
        Check == "Coverage" & Count > 0 ~ "WARNING",
        Check == "Row Count" & row_count_blocker ~ "FAIL",
        Check == "Row Count" & row_count_changed ~ "WARNING",
        Check == "Row Inflation" & Count > 2 ~ "FAIL",
        Check == "Row Inflation" & Count > 1.05 ~ "WARNING",
        Check == "Overlapping Variables" & Count > 0 ~ "WARNING",
        Check == "Unresolved Duplicate Variables" & Count > 0 ~ "FAIL",
        Check == "Variable Conflicts" & Count > 0 ~ "WARNING",
        Check == "Suspicious Conflicts" & Count > 0 ~ "WARNING",
        TRUE ~ "PASS"
      ),
      Details = dplyr::case_when(
        Check == "Key Types" & Count > 0 ~
          "At least one key has different storage classes across datasets. Keys were coerced to character for auditing.",
        Check == "Key Types" ~
          "Key storage classes match across datasets.",
        Check == "Missing Keys" & Count > 0 ~
          "At least one source or merged dataset contains rows with missing key values. These rows are excluded from duplicate-key and coverage calculations.",
        Check == "Missing Keys" ~
          "No missing key rows detected.",
        Check == "Expected Relationship" ~ relationship_details,
        Check == "Duplicate Keys" & relationship_duplicate_blockers > 0 ~
          "Duplicate complete key combinations violate the expected relationship.",
        Check == "Duplicate Keys" ~
          "No duplicate key blockers detected for the expected relationship.",
        Check == "Coverage" & Count > 0 ~
          "Some complete key combinations appear only in one source dataset.",
        Check == "Coverage" ~
          "All complete source key combinations overlap.",
        Check == "Row Count" & row_count_blocker ~
          "MergedData row count changed even though the expected relationship should preserve left-side row count.",
        Check == "Row Count" & row_count_changed ~
          "MergedData row count changed. This may be expected for one-to-many, many-to-many, or auto relationship checks.",
        Check == "Row Count" ~
          "MergedData row count matches LeftData row count.",
        Check == "Row Inflation" & Count > 2 ~
          "MergedData has more than twice as many rows as the larger source dataset. This may indicate an accidental many-to-many merge.",
        Check == "Row Inflation" & Count > 1.05 ~
          "MergedData has more rows than expected. Review whether row multiplication was intentional.",
        Check == "Row Inflation" ~
          "No meaningful row inflation detected.",
        Check == "Overlapping Variables" & Count > 0 ~
          "Variables appear in both source datasets but were not specified as keys.",
        Check == "Overlapping Variables" ~
          "No non-key variables overlap across source datasets.",
        Check == "Unresolved Duplicate Variables" & Count > 0 ~
          "MergedData still contains unresolved .x/.y or _x/_y variable pairs.",
        Check == "Unresolved Duplicate Variables" ~
          "No unresolved duplicate variable pairs detected.",
        Check == "Variable Conflicts" & Count > 0 ~
          "At least one duplicated variable pair contains conflicting values.",
        Check == "Variable Conflicts" ~
          "No duplicated-variable value conflicts detected.",
        Check == "Suspicious Conflicts" & Count > 0 ~
          "At least one duplicated variable has low agreement or mismatched classes.",
        Check == "Suspicious Conflicts" ~
          "No low-agreement or class-mismatched duplicated variables detected.",
        TRUE ~ ""
      )
    )

  ready_for_analysis <- !any(checks$Status == "FAIL", na.rm = TRUE)

  checks <- dplyr::bind_rows(
    checks,
    tibble::tibble(
      Check = "Merge Readiness",
      Count = ifelse(ready_for_analysis, 0, 1),
      Status = ifelse(ready_for_analysis, "PASS", "FAIL"),
      Details = ifelse(
        ready_for_analysis,
        "No major merge-integrity blockers detected.",
        "Major merge-integrity blockers detected. Review failed checks."
      )
    )
  )

  suggested_actions <- tibble::tibble(
    Priority = character(),
    Action = character()
  )

  if (key_type_count > 0) {
    key_type_actions <- key_types %>%
      dplyr::filter(
        LeftType != RightType |
          LeftType != MergedType |
          RightType != MergedType
      ) %>%
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Review key type mismatch for ",
          Key,
          " (Left: ",
          LeftType,
          "; Right: ",
          RightType,
          "; Merged: ",
          MergedType,
          ")."
        )
      ) %>%
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(suggested_actions, key_type_actions)
  }

  if (missing_key_total > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "Medium",
        Action = paste0(
          "Review missing key rows. Left: ",
          missing_key_rows_left,
          "; Right: ",
          missing_key_rows_right,
          "; Merged: ",
          missing_key_rows_merged,
          "."
        )
      )
    )
  }

  if (expected_relationship != "auto" && !relationship_matches_expected) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Review expected relationship. Expected ",
          expected_relationship,
          " but detected ",
          detected_relationship,
          "."
        )
      )
    )
  }

  if (relationship_duplicate_blockers > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Investigate duplicate key blockers for expected relationship: ",
          relationship_duplicate_blockers,
          " duplicated key group(s)."
        )
      )
    )
  }

  if (row_count_blocker) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = "Review row count change. This expected relationship should preserve the left-side row count."
      )
    )
  }

  if (nrow(potential_merge_risk) > 0) {
    overlap_actions <- potential_merge_risk %>%
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Review overlapping non-key variable: ",
          Variable,
          ". If dplyr join was run without `by =`, this variable may have been used as an unintended join key."
        )
      ) %>%
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(suggested_actions, overlap_actions)
  }

  if (nrow(unresolved_duplicate_variables) > 0) {
    dup_var_actions <- unresolved_duplicate_variables %>%
      dplyr::mutate(
        Priority = "High",
        Action = paste0("Resolve duplicated variable pair: ", Variable, ".")
      ) %>%
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(suggested_actions, dup_var_actions)
  }

  if (nrow(left_only_keys) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "Medium",
        Action = paste0(
          "Review ",
          nrow(left_only_keys),
          " complete key combination(s) present only in LeftData."
        )
      )
    )
  }

  if (nrow(right_only_keys) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "Medium",
        Action = paste0(
          "Review ",
          nrow(right_only_keys),
          " complete key combination(s) present only in RightData."
        )
      )
    )
  }

  if (nrow(variable_conflicts) > 0) {
    conflict_actions <- variable_conflicts %>%
      dplyr::count(Variable, name = "Conflicts") %>%
      dplyr::arrange(dplyr::desc(Conflicts), Variable) %>%
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Review ",
          Conflicts,
          " conflict(s) for duplicated variable ",
          Variable,
          "."
        )
      ) %>%
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(suggested_actions, conflict_actions)
  }

  if (nrow(suspicious_conflicts) > 0) {
    suspicious_actions <- suspicious_conflicts %>%
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Inspect suspicious duplicated variable ",
          Variable,
          " (Agreement: ",
          Agreement,
          "%; Left class: ",
          LeftClass,
          "; Right class: ",
          RightClass,
          ")."
        )
      ) %>%
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(suggested_actions, suspicious_actions)
  }

  left_match_rate <- ifelse(
    nrow(left_keys) > 0,
    round(100 * nrow(matching_keys) / nrow(left_keys), 1),
    NA_real_
  )

  right_match_rate <- ifelse(
    nrow(right_keys) > 0,
    round(100 * nrow(matching_keys) / nrow(right_keys), 1),
    NA_real_
  )

  summary_tbl <- tibble::tibble(
    LeftRows = nrow(LeftData),
    RightRows = nrow(RightData),
    MergedRows = nrow(MergedData),
    RowInflationFactor = row_inflation_factor,
    LeftColumns = ncol(LeftData),
    RightColumns = ncol(RightData),
    MergedColumns = ncol(MergedData),
    LeftUniqueKeys = nrow(left_keys),
    RightUniqueKeys = nrow(right_keys),
    MergedUniqueKeys = nrow(merged_keys),
    MatchingKeys = nrow(matching_keys),
    LeftOnlyKeys = nrow(left_only_keys),
    RightOnlyKeys = nrow(right_only_keys),
    LeftMatchRate = left_match_rate,
    RightMatchRate = right_match_rate,
    MissingKeyRows_Left = missing_key_rows_left,
    MissingKeyRows_Right = missing_key_rows_right,
    MissingKeyRows_Merged = missing_key_rows_merged,
    DuplicateKeyGroups_Left = duplicate_key_groups_left,
    DuplicateKeyGroups_Right = duplicate_key_groups_right,
    DuplicateKeyGroups_Merged = duplicate_key_groups_merged,
    DuplicateKeyBlockers = relationship_duplicate_blockers,
    ExpectedRelationship = expected_relationship,
    DetectedRelationship = detected_relationship,
    RelationshipMatchesExpected = relationship_matches_expected,
    KeyTypeMismatches = key_type_count,
    OverlappingVariables = nrow(overlapping_variables),
    UnresolvedDuplicateVariables = nrow(unresolved_duplicate_variables),
    VariableConflictCount = nrow(variable_conflicts),
    SuspiciousConflictCount = nrow(suspicious_conflicts),
    RowCountChanged = row_count_changed,
    RowCountBlocker = row_count_blocker,
    ReadyForAnalysis = ready_for_analysis
  )

  coverage_summary <- tibble::tibble(
    Metric = c(
      "Matching",
      "LeftOnly",
      "RightOnly",
      "LeftMatchRate",
      "RightMatchRate"
    ),
    Value = c(
      nrow(matching_keys),
      nrow(left_only_keys),
      nrow(right_only_keys),
      left_match_rate,
      right_match_rate
    )
  )

  overlap_text <- if (nrow(overlapping_variables) > 0) {
    paste(overlapping_variables$Variable, collapse = ", ")
  } else {
    "None"
  }

  unresolved_text <- if (nrow(unresolved_duplicate_variables) > 0) {
    paste(unresolved_duplicate_variables$Variable, collapse = ", ")
  } else {
    "None"
  }

  key_type_text <- if (key_type_count > 0) {
    key_types %>%
      dplyr::filter(
        LeftType != RightType |
          LeftType != MergedType |
          RightType != MergedType
      ) %>%
      dplyr::mutate(
        Text = paste0(
          Key,
          " (Left: ",
          LeftType,
          "; Right: ",
          RightType,
          "; Merged: ",
          MergedType,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "; ")
  } else {
    "None"
  }

  suspicious_text <- if (nrow(suspicious_conflicts) > 0) {
    paste(suspicious_conflicts$Variable, collapse = ", ")
  } else {
    "None"
  }

  summary_text <- paste0(
    "MERGE VALIDATION SUMMARY\n\n",
    "Keys:\n",
    paste(Keys, collapse = ", "),
    "\n\nRows:\n",
    "Left: ",
    nrow(LeftData),
    "\nRight: ",
    nrow(RightData),
    "\nMerged: ",
    nrow(MergedData),
    "\nRow Inflation Factor: ",
    row_inflation_factor,
    "\n\nUnique Complete Key Combinations:\n",
    "Left: ",
    nrow(left_keys),
    "\nRight: ",
    nrow(right_keys),
    "\nMerged: ",
    nrow(merged_keys),
    "\n\nRows With Missing Keys:\n",
    "Left: ",
    missing_key_rows_left,
    "\nRight: ",
    missing_key_rows_right,
    "\nMerged: ",
    missing_key_rows_merged,
    "\n\nCoverage:\n",
    "Matching: ",
    nrow(matching_keys),
    "\nLeft Only: ",
    nrow(left_only_keys),
    "\nRight Only: ",
    nrow(right_only_keys),
    "\nLeft Match Rate: ",
    left_match_rate,
    "%\nRight Match Rate: ",
    right_match_rate,
    "%\n\nDuplicate Complete Key Groups:\n",
    "Left: ",
    duplicate_key_groups_left,
    "\nRight: ",
    duplicate_key_groups_right,
    "\nMerged: ",
    duplicate_key_groups_merged,
    "\nBlockers For Expected Relationship: ",
    relationship_duplicate_blockers,
    "\n\nExpected Relationship:\n",
    expected_relationship,
    "\nDetected Relationship:\n",
    detected_relationship,
    "\nRelationship Matches Expected:\n",
    relationship_matches_expected,
    "\n\nKey Type Mismatches:\n",
    key_type_text,
    "\n\nOverlapping Variables Not In Keys:\n",
    overlap_text,
    "\n\nUnresolved Duplicate Variables:\n",
    unresolved_text,
    "\n\nVariable Conflicts:\n",
    nrow(variable_conflicts),
    "\n\nSuspicious Conflicts:\n",
    suspicious_text,
    "\n\nReady For Analysis:\n",
    ready_for_analysis
  )

  result <- list(
    SummaryText = summary_text,
    ReadyForAnalysis = ready_for_analysis,
    Summary = summary_tbl,
    Fingerprint = fingerprint,
    KeyTypes = key_types,
    Checks = checks,
    SuggestedActions = suggested_actions,
    Relationship = relationship,
    IDCoverage = id_coverage,
    DuplicateKeys = duplicate_keys,
    OverlappingVariables = overlapping_variables,
    PotentialMergeRisk = potential_merge_risk,
    JoinAudit = join_audit,
    OverlapAudit = overlap_audit,
    UnresolvedDuplicateVariables = unresolved_duplicate_variables,
    DuplicateVariables = duplicate_variables,
    SuspiciousConflicts = suspicious_conflicts,
    VariableConflicts = variable_conflicts,
    CoverageSummary = coverage_summary
  )

  return(result)
}

