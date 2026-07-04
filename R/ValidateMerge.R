#' Validate a merged dataset against its source datasets
#'
#' Audit a merged dataset against the two source datasets used to create it.
#' This function is designed to catch common merge problems before analysis,
#' including incompatible key types, duplicate key combinations, unexpected
#' overlapping variables, unresolved `.x` / `.y` variables, row inflation, and
#' value conflicts between duplicated variables.
#'
#' This function does not perform a merge. It assumes the user has already
#' created a merged dataset and wants to check whether the merge can be trusted.
#'
#' Important behavior: duplicate-key and coverage checks are calculated using
#' complete key combinations only. Rows with missing merge keys are counted and
#' reported, but they are not treated as duplicate key groups. This avoids
#' false failures when a source file has many records with missing optional IDs.
#'
#' @param LeftData A data frame used as one source for the merge.
#' @param RightData A data frame used as the other source for the merge.
#' @param MergedData The merged data frame to audit.
#' @param keys Character vector of key variables intended to define the merge.
#'   Multiple keys are supported, such as `c("study_id", "TimePoint")`.
#' @param Keys Deprecated. Use `keys` instead.
#'
#' @return A list with merge validation results, including:
#' \describe{
#'   \item{SummaryText}{A plain-text summary of the merge audit.}
#'   \item{ReadyForAnalysis}{Logical value indicating whether major merge-integrity issues were absent.}
#'   \item{Summary}{One-row tibble with core merge metrics.}
#'   \item{SummaryTable}{Metric/value summary table with Status as the first row.}
#'   \item{SummaryDisplay}{Formatted kable table with color-coded status.}
#'   \item{Fingerprint}{Tibble comparing rows, columns, and unique key combinations across datasets.}
#'   \item{KeyTypes}{Tibble showing key variable classes before coercion.}
#'   \item{Checks}{Tibble summarizing validation checks and pass/warning/fail status.}
#'   \item{SuggestedActions}{Tibble with suggested review steps.}
#'   \item{Relationship}{Tibble describing key uniqueness in LeftData, RightData, and MergedData.}
#'   \item{IDCoverage}{List containing Matching, LeftOnly, and RightOnly complete key combinations.}
#'   \item{DuplicateKeys}{List containing duplicated complete-key rows from Left, Right, and Merged datasets.}
#'   \item{OverlappingVariables}{Tibble of variables present in both source datasets but not listed as keys.}
#'   \item{PotentialMergeRisk}{Tibble of overlap variables that could have affected an unspecified dplyr join.}
#'   \item{JoinAudit}{Tibble showing variables present in both source datasets and whether they were specified as merge keys.}
#'   \item{OverlapAudit}{Tibble auditing all source variables by source presence and key status.}
#'   \item{UnresolvedDuplicateVariables}{Tibble of `.x` / `.y` duplicate variable pairs still present in MergedData.}
#'   \item{DuplicateVariables}{Tibble with agreement, conflict counts, and classes for duplicated variables.}
#'   \item{SuspiciousConflicts}{Subset of duplicated variables with low agreement or mismatched classes.}
#'   \item{VariableConflicts}{Long-format tibble of record-level value conflicts.}
#'   \item{CoverageSummary}{Tibble summarizing complete-key coverage.}
#' }
#'
#' @export
ValidateMerge <- function(
  LeftData,
  RightData,
  MergedData,
  keys,
  Keys = lifecycle::deprecated()
) {
  if (lifecycle::is_present(Keys)) {
    lifecycle::deprecate_warn(
      "19.15.0",
      "ValidateMerge(Keys)",
      "ValidateMerge(keys)"
    )
    keys <- Keys
  }

  if (missing(keys) || length(keys) == 0) {
    stop("`keys` must be supplied as a character vector.", call. = FALSE)
  }

  if (!is.character(keys)) {
    stop("`keys` must be a character vector.", call. = FALSE)
  }

  Keys <- keys

  if (!is.data.frame(LeftData)) {
    stop("`LeftData` must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(RightData)) {
    stop("`RightData` must be a data.frame.", call. = FALSE)
  }

  if (!is.data.frame(MergedData)) {
    stop("`MergedData` must be a data.frame.", call. = FALSE)
  }

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
        x_name <- paste0(key, ".x")
        y_name <- paste0(key, ".y")

        if (x_name %in% names(MergedData) && y_name %in% names(MergedData)) {
          paste0(
            key,
            " appears to be present as ",
            x_name,
            " and ",
            y_name,
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
      "\n\nCheck whether the merge was performed with the intended `by =` variables.",
      call. = FALSE
    )
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

  key_type_count <- key_types |>
    dplyr::filter(
      LeftType != RightType |
        LeftType != MergedType |
        RightType != MergedType
    ) |>
    nrow()

  for (key in Keys) {
    LeftData[[key]] <- as.character(LeftData[[key]])
    RightData[[key]] <- as.character(RightData[[key]])
    MergedData[[key]] <- as.character(MergedData[[key]])
  }

  complete_key_rows <- function(data, keys) {
    stats::complete.cases(data[, keys, drop = FALSE])
  }

  LeftComplete <- LeftData[complete_key_rows(LeftData, Keys), , drop = FALSE]
  RightComplete <- RightData[complete_key_rows(RightData, Keys), , drop = FALSE]
  MergedComplete <- MergedData[complete_key_rows(MergedData, Keys), , drop = FALSE]

  missing_key_rows <- tibble::tibble(
    Dataset = c("Left", "Right", "Merged"),
    Rows = c(
      nrow(LeftData) - nrow(LeftComplete),
      nrow(RightData) - nrow(RightComplete),
      nrow(MergedData) - nrow(MergedComplete)
    )
  )

  left_keys <- LeftComplete |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys)))

  right_keys <- RightComplete |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys)))

  merged_keys <- MergedComplete |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys)))

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
      "Rows with complete keys",
      "Rows with missing keys",
      "Columns",
      "Unique complete key combinations"
    ),
    Left = c(
      nrow(LeftData),
      nrow(LeftComplete),
      nrow(LeftData) - nrow(LeftComplete),
      ncol(LeftData),
      nrow(left_keys)
    ),
    Right = c(
      nrow(RightData),
      nrow(RightComplete),
      nrow(RightData) - nrow(RightComplete),
      ncol(RightData),
      nrow(right_keys)
    ),
    Merged = c(
      nrow(MergedData),
      nrow(MergedComplete),
      nrow(MergedData) - nrow(MergedComplete),
      ncol(MergedData),
      nrow(merged_keys)
    )
  )

  duplicate_key_rows <- function(data, keys) {
    if (nrow(data) == 0) {
      return(data |> dplyr::mutate(.n = integer()))
    }

    data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::mutate(.n = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::filter(.n > 1)
  }

  duplicate_keys_left <- duplicate_key_rows(LeftComplete, Keys)
  duplicate_keys_right <- duplicate_key_rows(RightComplete, Keys)
  duplicate_keys_merged <- duplicate_key_rows(MergedComplete, Keys)

  duplicate_key_groups_left <- duplicate_keys_left |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys))) |>
    nrow()

  duplicate_key_groups_right <- duplicate_keys_right |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys))) |>
    nrow()

  duplicate_key_groups_merged <- duplicate_keys_merged |>
    dplyr::distinct(dplyr::across(dplyr::all_of(Keys))) |>
    nrow()

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
  ) |>
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
    Variable = all_source_vars,
    InLeft = Variable %in% names(LeftData),
    InRight = Variable %in% names(RightData),
    IsKey = Variable %in% Keys,
    Risk = dplyr::case_when(
      InLeft & InRight & !IsKey ~ "Potential merge risk",
      InLeft & InRight & IsKey ~ "Expected key overlap",
      TRUE ~ "None"
    )
  )

  x_vars <- grep("\\.x$", names(MergedData), value = TRUE)

  duplicate_pairs <- purrr::map_dfr(
    x_vars,
    function(x_var) {
      base_var <- sub("\\.x$", "", x_var)
      y_var <- paste0(base_var, ".y")

      if (y_var %in% names(MergedData)) {
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

  if (nrow(duplicate_pairs) > 0) {
    unresolved_duplicate_variables <- duplicate_pairs |>
      dplyr::select(Variable) |>
      dplyr::arrange(Variable)
  } else {
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

        agree <- ifelse(
          is.na(x) & is.na(y),
          TRUE,
          ifelse(is.na(x) | is.na(y), FALSE, x == y)
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
    ) |>
      dplyr::arrange(dplyr::desc(Conflicts), Variable)
  } else {
    duplicate_variables <- tibble::tibble(
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

  suspicious_conflicts <- duplicate_variables |>
    dplyr::filter(
      Agreement < 75 |
        LeftClass != RightClass
    )

  if (nrow(duplicate_pairs) > 0) {
    variable_conflicts <- purrr::map_dfr(
      seq_len(nrow(duplicate_pairs)),
      function(i) {
        x_var <- duplicate_pairs$XVar[i]
        y_var <- duplicate_pairs$YVar[i]

        tmp <- MergedData |>
          dplyr::select(
            dplyr::all_of(Keys),
            dplyr::all_of(c(x_var, y_var))
          )

        x <- as.character(tmp[[x_var]])
        y <- as.character(tmp[[y_var]])

        agree <- ifelse(
          is.na(x) & is.na(y),
          TRUE,
          ifelse(is.na(x) | is.na(y), FALSE, x == y)
        )

        idx <- !agree

        if (!any(idx)) {
          return(NULL)
        }

        tmp[idx, Keys, drop = FALSE] |>
          tibble::as_tibble() |>
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
    ) |>
      dplyr::arrange(Variable)
  } else {
    variable_conflicts <- tibble::tibble(
      !!!stats::setNames(
        rep(list(character()), length(Keys)),
        Keys
      ),
      Variable = character(),
      LeftValue = character(),
      RightValue = character(),
      ConflictType = character()
    )
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
    UniqueKeys = c(left_unique, right_unique, merged_unique),
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
        ifelse(merged_unique, "unique complete keys", "duplicate complete keys")
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

  row_inflation_factor <- ifelse(
    nrow(LeftData) > 0,
    round(nrow(MergedData) / nrow(LeftData), 3),
    NA_real_
  )

  row_count_change <- nrow(MergedData) - nrow(LeftData)

  missing_key_total <- sum(missing_key_rows$Rows)

  checks <- tibble::tibble(
    Check = c(
      "Key Types",
      "Missing Keys",
      "Duplicate Keys",
      "Coverage",
      "Row Count Change",
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
      row_count_change,
      row_inflation_factor,
      overlap_total,
      unresolved_duplicate_total,
      conflict_total,
      suspicious_conflict_total
    )
  ) |>
    dplyr::mutate(
      Status = dplyr::case_when(
        Check == "Key Types" & Count > 0 ~ "WARNING",
        Check == "Missing Keys" & Count > 0 ~ "WARNING",
        Check == "Duplicate Keys" & Count > 0 ~ "FAIL",
        Check == "Coverage" & Count > 0 ~ "WARNING",
        Check == "Row Count Change" & Count > 0 ~ "FAIL",
        Check == "Row Count Change" & Count < 0 ~ "WARNING",
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
          "Rows with missing key values were detected. These rows were excluded from duplicate-key and coverage checks.",
        Check == "Missing Keys" ~
          "No missing key values detected.",

        Check == "Duplicate Keys" & Count > 0 ~
          "Duplicate complete-key combinations were detected.",
        Check == "Duplicate Keys" ~
          "No duplicate complete-key combinations detected.",

        Check == "Coverage" & Count > 0 ~
          "Some complete-key combinations appear only in one source dataset.",
        Check == "Coverage" ~
          "All complete source key combinations overlap.",

        Check == "Row Count Change" & Count > 0 ~
          "MergedData has more rows than LeftData. This usually indicates row multiplication.",
        Check == "Row Count Change" & Count < 0 ~
          "MergedData has fewer rows than LeftData. Confirm that row loss was intentional.",
        Check == "Row Count Change" ~
          "MergedData has the same number of rows as LeftData.",

        Check == "Row Inflation" & Count > 2 ~
          "MergedData has more than twice as many rows as LeftData. This may indicate an accidental many-to-many merge.",
        Check == "Row Inflation" & Count > 1.05 ~
          "MergedData has more rows than expected. Review whether row multiplication was intentional.",
        Check == "Row Inflation" ~
          "No meaningful row inflation detected.",

        Check == "Overlapping Variables" & Count > 0 ~
          "Variables appear in both source datasets but were not specified as keys.",
        Check == "Overlapping Variables" ~
          "No non-key variables overlap across source datasets.",

        Check == "Unresolved Duplicate Variables" & Count > 0 ~
          "MergedData still contains unresolved .x / .y variable pairs.",
        Check == "Unresolved Duplicate Variables" ~
          "No unresolved .x / .y variable pairs detected.",

        Check == "Variable Conflicts" & Count > 0 ~
          "At least one .x / .y variable pair contains conflicting values.",
        Check == "Variable Conflicts" ~
          "No .x / .y value conflicts detected.",

        Check == "Suspicious Conflicts" & Count > 0 ~
          "At least one duplicated variable has low agreement or mismatched classes.",
        Check == "Suspicious Conflicts" ~
          "No low-agreement or class-mismatched duplicated variables detected.",

        TRUE ~ ""
      )
    )

  ready_for_analysis <- !any(checks$Status == "FAIL")

  checks <- dplyr::bind_rows(
    checks,
    tibble::tibble(
      Check = "Merge Readiness",
      Count = ifelse(ready_for_analysis, 0, 1),
      Status = ifelse(ready_for_analysis, "PASS", "FAIL"),
      Details = ifelse(
        ready_for_analysis,
        "No major merge-integrity blockers detected.",
        "Major merge-integrity blockers detected. Review failed checks before analysis."
      )
    )
  )

  suggested_actions <- tibble::tibble(
    Priority = character(),
    Action = character()
  )

  if (key_type_count > 0) {
    key_type_actions <- key_types |>
      dplyr::filter(
        LeftType != RightType |
          LeftType != MergedType |
          RightType != MergedType
      ) |>
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Review key type mismatch for ",
          Key,
          " (Left: ", LeftType,
          "; Right: ", RightType,
          "; Merged: ", MergedType,
          ")."
        )
      ) |>
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      key_type_actions
    )
  }

  if (missing_key_total > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      missing_key_rows |>
        dplyr::filter(Rows > 0) |>
        dplyr::mutate(
          Priority = "Medium",
          Action = paste0(
            "Review ",
            Rows,
            " row(s) with missing key values in ",
            Dataset,
            ". These were excluded from duplicate-key and coverage checks."
          )
        ) |>
        dplyr::select(Priority, Action)
    )
  }

  if (duplicate_key_groups_left > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "Investigate duplicate complete-key combinations in LeftData: ",
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
          "Investigate duplicate complete-key combinations in RightData: ",
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
          "Investigate duplicate complete-key combinations in MergedData: ",
          duplicate_key_groups_merged,
          " duplicated key group(s)."
        )
      )
    )
  }

  if (row_count_change > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      tibble::tibble(
        Priority = "High",
        Action = paste0(
          "MergedData gained ",
          row_count_change,
          " row(s) relative to LeftData. Check for duplicate keys or an unintended many-to-many join."
        )
      )
    )
  }

  if (nrow(potential_merge_risk) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      potential_merge_risk |>
        dplyr::mutate(
          Priority = "Medium",
          Action = paste0(
            "Review overlapping non-key variable: ",
            Variable,
            ". If dplyr join was run without `by =`, this variable may have been used as an unintended join key."
          )
        ) |>
        dplyr::select(Priority, Action)
    )
  }

  if (nrow(unresolved_duplicate_variables) > 0) {
    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      unresolved_duplicate_variables |>
        dplyr::mutate(
          Priority = "High",
          Action = paste0(
            "Resolve duplicated variable pair: ",
            Variable,
            ".x / ",
            Variable,
            ".y."
          )
        ) |>
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
    conflict_actions <- variable_conflicts |>
      dplyr::count(Variable, name = "Conflicts") |>
      dplyr::arrange(dplyr::desc(Conflicts), Variable) |>
      dplyr::mutate(
        Priority = "Medium",
        Action = paste0(
          "Review ",
          Conflicts,
          " conflict(s) for duplicated variable ",
          Variable,
          "."
        )
      ) |>
      dplyr::select(Priority, Action)

    suggested_actions <- dplyr::bind_rows(
      suggested_actions,
      conflict_actions
    )
  }

  if (nrow(suspicious_conflicts) > 0) {
    suspicious_actions <- suspicious_conflicts |>
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
      ) |>
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

  status <- dplyr::case_when(
    any(checks$Status == "FAIL") ~ "FAIL",
    any(checks$Status == "WARNING") ~ "WARNING",
    TRUE ~ "PASS"
  )

  summary_tbl <- tibble::tibble(
    LeftRows = nrow(LeftData),
    RightRows = nrow(RightData),
    MergedRows = nrow(MergedData),
    RowCountChange = row_count_change,
    RowInflationFactor = row_inflation_factor,
    LeftColumns = ncol(LeftData),
    RightColumns = ncol(RightData),
    MergedColumns = ncol(MergedData),
    LeftCompleteKeyRows = nrow(LeftComplete),
    RightCompleteKeyRows = nrow(RightComplete),
    MergedCompleteKeyRows = nrow(MergedComplete),
    LeftMissingKeyRows = nrow(LeftData) - nrow(LeftComplete),
    RightMissingKeyRows = nrow(RightData) - nrow(RightComplete),
    MergedMissingKeyRows = nrow(MergedData) - nrow(MergedComplete),
    LeftUniqueKeys = nrow(left_keys),
    RightUniqueKeys = nrow(right_keys),
    MergedUniqueKeys = nrow(merged_keys),
    MatchingKeys = nrow(matching_keys),
    LeftOnlyKeys = nrow(left_only_keys),
    RightOnlyKeys = nrow(right_only_keys),
    LeftMatchRate = left_match_rate,
    RightMatchRate = right_match_rate,
    DuplicateKeyGroups_Left = duplicate_key_groups_left,
    DuplicateKeyGroups_Right = duplicate_key_groups_right,
    DuplicateKeyGroups_Merged = duplicate_key_groups_merged,
    KeyTypeMismatches = key_type_count,
    OverlappingVariables = nrow(overlapping_variables),
    UnresolvedDuplicateVariables = nrow(unresolved_duplicate_variables),
    VariableConflictCount = nrow(variable_conflicts),
    SuspiciousConflictCount = nrow(suspicious_conflicts),
    Status = status,
    ReadyForAnalysis = ready_for_analysis
  )

  summary_table <- tibble::tibble(
    Metric = c(
      "Status",
      "Rows (left -> merged)",
      "Columns (left -> merged)",
      "Complete-key rows (left / right / merged)",
      "Missing-key rows (left / right / merged)",
      "Keys matched",
      "Left-only keys",
      "Right-only keys",
      "Duplicate key groups (left / right / merged)",
      "Row inflation factor",
      "Ready for analysis"
    ),
    Value = c(
      status,
      paste0(nrow(LeftData), " -> ", nrow(MergedData)),
      paste0(ncol(LeftData), " -> ", ncol(MergedData)),
      paste0(nrow(LeftComplete), " / ", nrow(RightComplete), " / ", nrow(MergedComplete)),
      paste0(
        nrow(LeftData) - nrow(LeftComplete),
        " / ",
        nrow(RightData) - nrow(RightComplete),
        " / ",
        nrow(MergedData) - nrow(MergedComplete)
      ),
      paste0(nrow(matching_keys), " / ", nrow(left_keys)),
      nrow(left_only_keys),
      nrow(right_only_keys),
      paste0(
        duplicate_key_groups_left,
        " / ",
        duplicate_key_groups_right,
        " / ",
        duplicate_key_groups_merged
      ),
      row_inflation_factor,
      ready_for_analysis
    )
  )

  make_summary_display <- function(summary_table, status) {
    status_background <- dplyr::case_when(
      status == "PASS" ~ "#2E7D32",
      status == "WARNING" ~ "#F9A825",
      status == "FAIL" ~ "#C62828",
      TRUE ~ "#616161"
    )

    status_text <- dplyr::case_when(
      status == "PASS" ~ "PASS",
      status == "WARNING" ~ "WARNING",
      status == "FAIL" ~ "FAIL",
      TRUE ~ status
    )

    display_table <- summary_table
    display_table$Value[display_table$Metric == "Status"] <- status_text

    if (requireNamespace("kableExtra", quietly = TRUE)) {
      display_table |>
        knitr::kable(
          escape = FALSE,
          align = c("l", "l")
        ) |>
        kableExtra::kable_styling(
          full_width = FALSE
        ) |>
        kableExtra::row_spec(
          1,
          bold = TRUE,
          color = "white",
          background = status_background
        )
    } else {
      knitr::kable(display_table)
    }
  }

  summary_display <- make_summary_display(summary_table, status)

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
    key_types |>
      dplyr::filter(
        LeftType != RightType |
          LeftType != MergedType |
          RightType != MergedType
      ) |>
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
      ) |>
      dplyr::pull(Text) |>
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
    "Status:\n",
    status, "\n\n",

    "Ready For Analysis:\n",
    ready_for_analysis, "\n\n",

    "Keys:\n",
    paste(Keys, collapse = ", "), "\n\n",

    "Rows:\n",
    "Left: ", nrow(LeftData), "\n",
    "Right: ", nrow(RightData), "\n",
    "Merged: ", nrow(MergedData), "\n",
    "Row Count Change: ", row_count_change, "\n",
    "Row Inflation Factor: ", row_inflation_factor, "\n\n",

    "Complete Key Rows:\n",
    "Left: ", nrow(LeftComplete), "\n",
    "Right: ", nrow(RightComplete), "\n",
    "Merged: ", nrow(MergedComplete), "\n\n",

    "Missing Key Rows:\n",
    "Left: ", nrow(LeftData) - nrow(LeftComplete), "\n",
    "Right: ", nrow(RightData) - nrow(RightComplete), "\n",
    "Merged: ", nrow(MergedData) - nrow(MergedComplete), "\n\n",

    "Unique Complete Key Combinations:\n",
    "Left: ", nrow(left_keys), "\n",
    "Right: ", nrow(right_keys), "\n",
    "Merged: ", nrow(merged_keys), "\n\n",

    "Coverage:\n",
    "Matching: ", nrow(matching_keys), "\n",
    "Left Only: ", nrow(left_only_keys), "\n",
    "Right Only: ", nrow(right_only_keys), "\n",
    "Left Match Rate: ", left_match_rate, "%\n",
    "Right Match Rate: ", right_match_rate, "%\n\n",

    "Duplicate Complete-Key Groups:\n",
    "Left: ", duplicate_key_groups_left, "\n",
    "Right: ", duplicate_key_groups_right, "\n",
    "Merged: ", duplicate_key_groups_merged, "\n\n",

    "Key Type Mismatches:\n",
    key_type_text, "\n\n",

    "Overlapping Variables Not In Keys:\n",
    overlap_text, "\n\n",

    "Unresolved Duplicate Variables:\n",
    unresolved_text, "\n\n",

    "Variable Conflicts:\n",
    nrow(variable_conflicts), "\n\n",

    "Suspicious Conflicts:\n",
    suspicious_text, "\n\n",

    "Detected Source Relationship:\n",
    relationship_type
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

  result <- list(
    SummaryText = summary_text,
    ReadyForAnalysis = ready_for_analysis,
    Status = status,
    Summary = summary_tbl,
    SummaryTable = summary_table,
    SummaryDisplay = summary_display,
    Fingerprint = fingerprint,
    KeyTypes = key_types,
    MissingKeyRows = missing_key_rows,
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
