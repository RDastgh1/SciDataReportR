#' Validate a merged dataset against its source datasets
#'
#' Audit a merged dataset against the two source datasets used to create it.
#' This function is designed to catch common merge problems before analysis,
#' including incompatible key types, duplicate key combinations, incomplete
#' key coverage, unexpected overlapping variables, unresolved `.x` / `.y`
#' variables, and value conflicts between duplicated variables.
#'
#' This function does not perform a merge. It assumes the user has already
#' created a merged dataset and wants to check whether the merge can be trusted.
#'
#' @param LeftData A data frame used as one source for the merge.
#' @param RightData A data frame used as the other source for the merge.
#' @param MergedData The merged data frame to audit.
#' @param keys Character vector of key variables intended to define the merge.
#'   Multiple keys are supported, such as `c("study_id", "visit")`.
#' @param Keys Deprecated. Use `keys` instead.
#'
#' @return A list with merge validation results.
#'
#' @export
ValidateMerge <- function(LeftData,
                          RightData,
                          MergedData,
                          keys,
                          Keys = lifecycle::deprecated()) {

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

        if (x_name_dot %in% names(MergedData) && y_name_dot %in% names(MergedData)) {
          paste0(
            key,
            " appears to be present as ",
            x_name_dot,
            " and ",
            y_name_dot,
            ". This often happens when the merge was performed using fewer keys than intended."
          )
        } else if (x_name_under %in% names(MergedData) && y_name_under %in% names(MergedData)) {
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
      "\n\n",
      "Check whether the merge was performed with the intended `by = ` variables.",
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
    key_cols <- stats::setNames(
      rep(list(character()), length(keys)),
      keys
    )

    tibble::tibble(
      !!!key_cols,
      Variable = character(),
      LeftValue = character(),
      RightValue = character(),
      ConflictType = character()
    )
  }

  empty_key_table <- function(keys) {
    tibble::as_tibble(
      stats::setNames(
        rep(list(character()), length(keys)),
        keys
      )
    )
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
      dplyr::distinct(
        dplyr::across(dplyr::all_of(keys))
      )
  }

  duplicated_key_rows <- function(df, keys) {
    if (nrow(df) == 0) {
      out <- empty_key_table(keys)
      out$.n <- integer()
      return(out)
    }

    df %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(keys))
      ) %>%
      dplyr::mutate(.n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.n > 1)
  }

  count_duplicate_key_groups <- function(df, keys) {
    if (nrow(df) == 0) {
      return(0L)
    }

    df %>%
      dplyr::distinct(
        dplyr::across(dplyr::all_of(keys))
      ) %>%
      nrow()
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

  matching_keys <- dplyr::inner_join(
    left_keys,
    right_keys,
    by = Keys
  )

  left_only_keys <- dplyr::anti_join(
    left_keys,
    right_keys,
    by = Keys
  )

  right_only_keys <- dplyr::anti_join(
    right_keys,
    left_keys,
    by = Keys
  )

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

  duplicate_keys_left <- duplicated_key_rows(LeftData_keyed, Keys)
  duplicate_keys_right <- duplicated_key_rows(RightData_keyed, Keys)
  duplicate_keys_merged <- duplicated_key_rows(MergedData_keyed, Keys)

  duplicate_key_groups_left <- count_duplicate_key_groups(duplicate_keys_left, Keys)
  duplicate_key_groups_right <- count_duplicate_key_groups(duplicate_keys_right, Keys)
  duplicate_key_groups_merged <- count_duplicate_key_groups(duplicate_keys_merged, Keys)

  duplicate_keys <- list(
    Left = duplicate_keys_left,
    Right = duplicate_keys_right,
    Merged = duplicate_keys_merged
  )

  source_overlap_vars <- intersect(
    names(LeftData),
    names(RightData)
  )

  overlapping_non_key_vars <- setdiff(
    source_overlap_vars,
    Keys
  )

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

  all_source_vars <- sort(unique(c(
    names(LeftData),
    names(RightData)
  )))

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
      return(tibble::tibble(
        Variable = character(),
        XVar = character(),
        YVar = character()
      ))
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
    unresolved_duplicate_variables <- tibble::tibble(
      Variable = character()
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
      dplyr::filter(
        Agreement < 75 |
          LeftClass != RightClass
      ) %>%
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
          dplyr::select(
            dplyr::all_of(Keys),
            dplyr::all_of(c(x_var, y_var))
          )

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

  left_unique <- duplicate_key_groups_left == 0
  right_unique <- duplicate_key_groups_right == 0
  merged_unique <- duplicate_key_groups_merged == 0

  relationship_type <- dplyr::case_when(
    left_unique & right_unique ~ "one-to-one",
    left_unique & !right_unique ~ "one-to-many",
    !left_unique & right_unique ~ "many-to-one",
    TRUE ~ "many-to-many"
  )

  relationship <- tibble::tibble(
    Dataset = c("Left", "Right", "Merged"),
    UniqueKeys = c(
      left_unique,
      right_unique,
      merged_unique
    ),
    DuplicateKeyGroups = c(
      duplicate_key_groups_left,
      duplicate_key_groups_right,
      duplicate_key_groups_merged
    ),
    Relationship = c(
      relationship_type,
      relationship_type,
      paste0(
        "Merged result: ",
        ifelse(merged_unique, "unique keys", "duplicate keys")
      )
    )
  )

  duplicate_key_total <- duplicate_key_groups_left +
    duplicate_key_groups_right +
    duplicate_key_groups_merged

  coverage_total <- nrow(left_only_keys) + nrow(right_only_keys)

  overlap_total <- nrow(overlapping_variables)

  unresolved_duplicate_total <- nrow(unresolved_duplicate_variables)

  conflict_total <- nrow(variable_conflicts)

  suspicious_conflict_total <- nrow(suspicious_conflicts)

  missing_key_total <- missing_key_rows_left +
    missing_key_rows_right +
    missing_key_rows_merged

  row_inflation_factor <- ifelse(
    max(nrow(LeftData), nrow(RightData)) > 0,
    round(nrow(MergedData) / max(nrow(LeftData), nrow(RightData)), 3),
    NA_real_
  )

  checks <- tibble::tibble(
    Check = c(
      "Key Types",
      "Missing Keys",
      "Duplicate Keys",
      "Coverage",
      "Row Inflation",
      "Overlapping Variables",
      "Unresolved Duplicate Variables",
      "Variable Conflicts",
      "Suspicious Conflicts"
    ),
    Count = c(
      key_type_count,
      missing_key_total,
      duplicate_key_total,
      coverage_total,
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
        Check == "Duplicate Keys" & Count > 0 ~ "FAIL",
        Check == "Coverage" & Count > 0 ~ "WARNING",
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

        Check == "Duplicate Keys" & Count > 0 ~
          "Duplicate complete key combinations were detected.",
        Check == "Duplicate Keys" ~
          "No duplicate complete key combinations detected.",

        Check == "Coverage" & Count > 0 ~
          "Some complete key combinations appear only in one source dataset.",
        Check == "Coverage" ~
          "All complete source key combinations overlap.",

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

  ready_for_analysis <- duplicate_key_total == 0 &&
    unresolved_duplicate_total == 0

  checks <- dplyr::bind_rows(
    checks,
    tibble::tibble(
      Check = "Merge Readiness",
      Count = ifelse(ready_for_analysis, 0, 1),
      Status = ifelse(ready_for_analysis, "PASS", "FAIL"),
      Details = ifelse(
        ready_for_analysis,
        "No major merge-integrity blockers detected.",
        "Major merge-integrity blockers detected. Review duplicate keys and unresolved duplicate variables."
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

    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      key_type_actions
    )
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

  if (duplicate_key_groups_left > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Investigate duplicate complete key combinations in LeftData: ",
          duplicate_key_groups_left,
          " duplicated key group(s)."
        )
      )
    )
  }

  if (duplicate_key_groups_right > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Investigate duplicate complete key combinations in RightData: ",
          duplicate_key_groups_right,
          " duplicated key group(s)."
        )
      )
    )
  }

  if (duplicate_key_groups_merged > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Investigate duplicate complete key combinations in MergedData: ",
          duplicate_key_groups_merged,
          " duplicated key group(s)."
        )
      )
    )
  }

  if (nrow(potential_merge_risk) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      potential_merge_risk %>%
        dplyr::mutate(
          Priority = "Medium",
          Action = paste0(
            "Review overlapping non-key variable: ",
            Variable,
            ". If dplyr join was run without `by =`, this variable may have been used as an unintended join key."
          )
        ) %>%
        dplyr::select(Priority, Action)
    )
  }

  if (nrow(unresolved_duplicate_variables) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      unresolved_duplicate_variables %>%
        dplyr::mutate(
          Priority = "High",
          Action = paste0(
            "Resolve duplicated variable pair: ",
            Variable,
            "."
          )
        ) %>%
        dplyr::select(Priority, Action)
    )
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

    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      conflict_actions
    )
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

    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      suspicious_actions
    )
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
    KeyTypeMismatches = key_type_count,
    OverlappingVariables = nrow(overlapping_variables),
    UnresolvedDuplicateVariables = nrow(unresolved_duplicate_variables),
    VariableConflictCount = nrow(variable_conflicts),
    SuspiciousConflictCount = nrow(suspicious_conflicts),
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
    "\n\n",

    "Rows:\n",
    "Left: ",
    nrow(LeftData),
    "\n",
    "Right: ",
    nrow(RightData),
    "\n",
    "Merged: ",
    nrow(MergedData),
    "\n",
    "Row Inflation Factor: ",
    row_inflation_factor,
    "\n\n",

    "Unique Complete Key Combinations:\n",
    "Left: ",
    nrow(left_keys),
    "\n",
    "Right: ",
    nrow(right_keys),
    "\n",
    "Merged: ",
    nrow(merged_keys),
    "\n\n",

    "Rows With Missing Keys:\n",
    "Left: ",
    missing_key_rows_left,
    "\n",
    "Right: ",
    missing_key_rows_right,
    "\n",
    "Merged: ",
    missing_key_rows_merged,
    "\n\n",

    "Coverage:\n",
    "Matching: ",
    nrow(matching_keys),
    "\n",
    "Left Only: ",
    nrow(left_only_keys),
    "\n",
    "Right Only: ",
    nrow(right_only_keys),
    "\n",
    "Left Match Rate: ",
    left_match_rate,
    "%\n",
    "Right Match Rate: ",
    right_match_rate,
    "%\n\n",

    "Duplicate Complete Key Groups:\n",
    "Left: ",
    duplicate_key_groups_left,
    "\n",
    "Right: ",
    duplicate_key_groups_right,
    "\n",
    "Merged: ",
    duplicate_key_groups_merged,
    "\n\n",

    "Key Type Mismatches:\n",
    key_type_text,
    "\n\n",

    "Overlapping Variables Not In Keys:\n",
    overlap_text,
    "\n\n",

    "Unresolved Duplicate Variables:\n",
    unresolved_text,
    "\n\n",

    "Variable Conflicts:\n",
    nrow(variable_conflicts),
    "\n\n",

    "Suspicious Conflicts:\n",
    suspicious_text,
    "\n\n",

    "Detected Source Relationship:\n",
    relationship_type,
    "\n\n",

    "Ready For Analysis:\n",
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
