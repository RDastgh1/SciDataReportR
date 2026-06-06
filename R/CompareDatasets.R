#' Compare two versions of a dataset
#'
#' Compare an old dataset and a new dataset using one or more key variables.
#' This function identifies record-level, variable-level, and cell-level changes
#' between dataset versions. It is useful when reviewing updated data extracts,
#' revised REDCap exports, cleaned spreadsheet versions, or vendor-delivered
#' dataset updates.
#'
#' The function always returns both detailed cell-level changes and summary-level
#' outputs. Cell-level changes are stored in long format, with one row per
#' changed value.
#'
#' Key variables are coerced to character internally before comparison because
#' IDs are often stored as numeric in one file and character in another. The
#' original key classes are preserved in the `KeyTypes` output.
#'
#' Variable names are also audited for common tibble/readxl name-repair suffixes
#' such as `...372`. These suffixes are ignored when identifying variables added
#' or removed, while raw variable names are still preserved in the output.
#'
#' If duplicate key combinations are detected, the function still runs, but
#' cell-level comparison is performed only for key combinations that are unique
#' in both datasets. Duplicate keys are returned separately and flagged in
#' `Checks`.
#'
#' @param OldData A data frame representing the earlier dataset version.
#' @param NewData A data frame representing the newer dataset version.
#' @param Keys Character vector of key variables used to align records across
#'   the two datasets. Multiple keys are supported, such as
#'   `c("study_id", "TimePoint")`.
#'
#' @return A list with dataset comparison results, including:
#' \describe{
#'   \item{SummaryText}{A plain-text summary of the dataset comparison.}
#'   \item{Summary}{One-row tibble with core comparison metrics.}
#'   \item{Fingerprint}{Tibble comparing rows, columns, and unique key combinations.}
#'   \item{KeyTypes}{Tibble showing key variable classes before coercion.}
#'   \item{Checks}{Tibble summarizing comparison checks and pass/warning/fail status.}
#'   \item{StructureChanges}{Tibble of variables added to or removed from NewData, using normalized variable names.}
#'   \item{AddedVariables}{Tibble of variables present in NewData but not OldData, using normalized variable names.}
#'   \item{RemovedVariables}{Tibble of variables present in OldData but not NewData, using normalized variable names.}
#'   \item{AddedRecords}{Tibble of key combinations present in NewData but not OldData.}
#'   \item{RemovedRecords}{Tibble of key combinations present in OldData but not NewData.}
#'   \item{DuplicateKeys}{List containing duplicated key rows from OldData and NewData.}
#'   \item{ComparisonKeys}{List describing matching keys, compared keys, and keys excluded from cell comparison due to duplicate key combinations.}
#'   \item{NameRepairAudit}{Tibble describing variables whose raw names differ after removing tibble-style `...number` suffixes.}
#'   \item{ComparisonVariableMap}{Tibble mapping normalized variable names to the raw OldData and NewData names used for cell-level comparison.}
#'   \item{ClassAudit}{Tibble comparing variable classes for common non-key variables.}
#'   \item{ModifiedValues}{Long-format tibble of cell-level value changes.}
#'   \item{VariableChangeSummary}{Tibble summarizing changes by variable.}
#'   \item{TopChangedVariables}{Top changed variables by number of modified values.}
#'   \item{SuspiciousChanges}{Tibble of high-change-rate or class-change variables.}
#' }
#'
#' @export
CompareDatasets <- function(
    OldData,
    NewData,
    Keys
) {

  # Validate inputs

  if (!is.data.frame(OldData)) {
    stop("OldData must be a data.frame.")
  }

  if (!is.data.frame(NewData)) {
    stop("NewData must be a data.frame.")
  }

  if (missing(Keys) || length(Keys) == 0) {
    stop("Keys must be supplied as a character vector.")
  }

  if (!is.character(Keys)) {
    stop("Keys must be a character vector.")
  }

  missing_old <- setdiff(Keys, names(OldData))
  missing_new <- setdiff(Keys, names(NewData))

  if (length(missing_old) > 0) {
    stop(
      "The following key variable(s) are missing from OldData: ",
      paste(missing_old, collapse = ", ")
    )
  }

  if (length(missing_new) > 0) {
    stop(
      "The following key variable(s) are missing from NewData: ",
      paste(missing_new, collapse = ", ")
    )
  }

  # Preserve original key classes before coercion

  key_types <- tibble::tibble(
    Key = Keys,
    OldType = purrr::map_chr(
      Keys,
      ~ paste(class(OldData[[.x]]), collapse = "/")
    ),
    NewType = purrr::map_chr(
      Keys,
      ~ paste(class(NewData[[.x]]), collapse = "/")
    )
  )

  key_type_count <- key_types %>%
    dplyr::filter(
      OldType != NewType
    ) %>%
    nrow()

  # Standardize key columns

  # Keys are converted to character because IDs are often numeric in one export
  # and character in another. The comparison is about record identity, not
  # preserving original key storage classes.
  for (key in Keys) {
    OldData[[key]] <- as.character(OldData[[key]])
    NewData[[key]] <- as.character(NewData[[key]])
  }

  # Prepare key tables

  old_keys <- OldData %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    )

  new_keys <- NewData %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    )

  matching_keys <- dplyr::inner_join(
    old_keys,
    new_keys,
    by = Keys
  )

  added_records <- dplyr::anti_join(
    new_keys,
    old_keys,
    by = Keys
  )

  removed_records <- dplyr::anti_join(
    old_keys,
    new_keys,
    by = Keys
  )

  # Detect duplicate keys

  duplicate_keys_old <- OldData %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(Keys))
    ) %>%
    dplyr::mutate(.n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.n > 1)

  duplicate_keys_new <- NewData %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(Keys))
    ) %>%
    dplyr::mutate(.n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.n > 1)

  duplicate_key_groups_old <- duplicate_keys_old %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    ) %>%
    nrow()

  duplicate_key_groups_new <- duplicate_keys_new %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    ) %>%
    nrow()

  duplicate_keys <- list(
    Old = duplicate_keys_old,
    New = duplicate_keys_new
  )

  duplicate_old_keys <- duplicate_keys_old %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    )

  duplicate_new_keys <- duplicate_keys_new %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    )

  duplicate_any_keys <- dplyr::bind_rows(
    duplicate_old_keys,
    duplicate_new_keys
  ) %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(Keys))
    )

  compared_keys <- dplyr::anti_join(
    matching_keys,
    duplicate_any_keys,
    by = Keys
  )

  excluded_duplicate_keys <- dplyr::inner_join(
    matching_keys,
    duplicate_any_keys,
    by = Keys
  )

  comparison_keys <- list(
    Matching = matching_keys,
    Compared = compared_keys,
    ExcludedDuplicateKeys = excluded_duplicate_keys
  )

  # Build fingerprint

  fingerprint <- tibble::tibble(
    Metric = c(
      "Rows",
      "Columns",
      "Unique Key Combinations"
    ),
    Old = c(
      nrow(OldData),
      ncol(OldData),
      nrow(old_keys)
    ),
    New = c(
      nrow(NewData),
      ncol(NewData),
      nrow(new_keys)
    )
  )

  # Audit variable names

  # Audit variable names

  old_raw_names <- names(OldData)
  old_base_names <- stringr::str_remove(
    old_raw_names,
    "\\.\\.\\.[0-9]+$"
  )

  old_base_names <- dplyr::if_else(
    old_base_names == "",
    old_raw_names,
    old_base_names
  )

  new_raw_names <- names(NewData)
  new_base_names <- stringr::str_remove(
    new_raw_names,
    "\\.\\.\\.[0-9]+$"
  )

  new_base_names <- dplyr::if_else(
    new_base_names == "",
    new_raw_names,
    new_base_names
  )

  old_variable_map <- tibble::tibble(
    OldVariable = old_raw_names,
    Variable = old_base_names
  )

  new_variable_map <- tibble::tibble(
    NewVariable = new_raw_names,
    Variable = new_base_names
  )

  old_variable_bases <- old_variable_map %>%
    dplyr::distinct(Variable)

  new_variable_bases <- new_variable_map %>%
    dplyr::distinct(Variable)

  added_variable_bases <- dplyr::anti_join(
    new_variable_bases,
    old_variable_bases,
    by = "Variable"
  )

  removed_variable_bases <- dplyr::anti_join(
    old_variable_bases,
    new_variable_bases,
    by = "Variable"
  )

  added_variables <- new_variable_map %>%
    dplyr::semi_join(
      added_variable_bases,
      by = "Variable"
    ) %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(
      NewVariables = paste(
        sort(unique(NewVariable)),
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Variable)

  removed_variables <- old_variable_map %>%
    dplyr::semi_join(
      removed_variable_bases,
      by = "Variable"
    ) %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(
      OldVariables = paste(
        sort(unique(OldVariable)),
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Variable)

  structure_changes <- dplyr::bind_rows(
    added_variables %>%
      dplyr::mutate(
        Change = "Added",
        SourceVariables = NewVariables
      ) %>%
      dplyr::select(
        Change,
        Variable,
        SourceVariables
      ),
    removed_variables %>%
      dplyr::mutate(
        Change = "Removed",
        SourceVariables = OldVariables
      ) %>%
      dplyr::select(
        Change,
        Variable,
        SourceVariables
      )
  )

  if (nrow(structure_changes) == 0) {
    structure_changes <- tibble::tibble(
      Change = character(),
      Variable = character(),
      SourceVariables = character()
    )
  }

  common_variable_bases <- dplyr::inner_join(
    old_variable_bases,
    new_variable_bases,
    by = "Variable"
  ) %>%
    dplyr::filter(
      !Variable %in% Keys
    )

  name_repair_audit <- purrr::map_dfr(
    common_variable_bases$Variable,
    function(variable) {

      old_vars <- old_variable_map %>%
        dplyr::filter(
          Variable == variable
        ) %>%
        dplyr::pull(OldVariable) %>%
        unique() %>%
        sort()

      new_vars <- new_variable_map %>%
        dplyr::filter(
          Variable == variable
        ) %>%
        dplyr::pull(NewVariable) %>%
        unique() %>%
        sort()

      tibble::tibble(
        Variable = variable,
        OldVariableCount = length(old_vars),
        NewVariableCount = length(new_vars),
        OldVariables = paste(old_vars, collapse = ", "),
        NewVariables = paste(new_vars, collapse = ", "),
        NameChanged = paste(old_vars, collapse = " | ") !=
          paste(new_vars, collapse = " | "),
        AmbiguousForComparison = length(old_vars) != 1 ||
          length(new_vars) != 1
      )
    }
  ) %>%
    dplyr::filter(
      NameChanged | AmbiguousForComparison
    ) %>%
    dplyr::arrange(Variable)

  if (nrow(name_repair_audit) == 0) {
    name_repair_audit <- tibble::tibble(
      Variable = character(),
      OldVariableCount = integer(),
      NewVariableCount = integer(),
      OldVariables = character(),
      NewVariables = character(),
      NameChanged = logical(),
      AmbiguousForComparison = logical()
    )
  }

  comparison_variable_map <- purrr::map_dfr(
    common_variable_bases$Variable,
    function(variable) {

      old_vars <- old_variable_map %>%
        dplyr::filter(
          Variable == variable
        ) %>%
        dplyr::pull(OldVariable) %>%
        unique() %>%
        sort()

      new_vars <- new_variable_map %>%
        dplyr::filter(
          Variable == variable
        ) %>%
        dplyr::pull(NewVariable) %>%
        unique() %>%
        sort()

      if (length(old_vars) == 1 && length(new_vars) == 1) {
        tibble::tibble(
          Variable = variable,
          OldVariable = old_vars,
          NewVariable = new_vars,
          OldClass = paste(class(OldData[[old_vars]]), collapse = "/"),
          NewClass = paste(class(NewData[[new_vars]]), collapse = "/"),
          ClassChanged = paste(class(OldData[[old_vars]]), collapse = "/") !=
            paste(class(NewData[[new_vars]]), collapse = "/"),
          NameChanged = old_vars != new_vars,
          Compared = TRUE,
          SkipReason = NA_character_
        )
      } else {
        tibble::tibble(
          Variable = variable,
          OldVariable = paste(old_vars, collapse = ", "),
          NewVariable = paste(new_vars, collapse = ", "),
          OldClass = NA_character_,
          NewClass = NA_character_,
          ClassChanged = NA,
          NameChanged = paste(old_vars, collapse = " | ") !=
            paste(new_vars, collapse = " | "),
          Compared = FALSE,
          SkipReason = "Variable name maps ambiguously after removing name-repair suffixes."
        )
      }
    }
  ) %>%
    dplyr::arrange(Variable)

  comparison_variable_map <- comparison_variable_map %>%
    dplyr::filter(
      !Variable %in% Keys
    )

  class_audit <- comparison_variable_map %>%
    dplyr::select(
      Variable,
      OldVariable,
      NewVariable,
      OldClass,
      NewClass,
      ClassChanged,
      NameChanged,
      Compared,
      SkipReason
    )

  class_change_count <- class_audit %>%
    dplyr::filter(
      Compared,
      ClassChanged
    ) %>%
    nrow()

  skipped_variable_count <- comparison_variable_map %>%
    dplyr::filter(
      !Compared
    ) %>%
    nrow()

  comparison_variables <- comparison_variable_map %>%
    dplyr::filter(
      Compared
    )

  # Compare cell-level values

  numeric_tolerance <- 1e-8

  if (nrow(comparison_variables) > 0 && nrow(compared_keys) > 0) {

    modified_values <- purrr::map_dfr(
      seq_len(nrow(comparison_variables)),
      function(i) {

        variable <- comparison_variables$Variable[i]
        old_variable <- comparison_variables$OldVariable[i]
        new_variable <- comparison_variables$NewVariable[i]

        old_tmp <- OldData %>%
          dplyr::semi_join(
            compared_keys,
            by = Keys
          ) %>%
          dplyr::select(
            dplyr::all_of(Keys),
            dplyr::all_of(old_variable)
          )

        new_tmp <- NewData %>%
          dplyr::semi_join(
            compared_keys,
            by = Keys
          ) %>%
          dplyr::select(
            dplyr::all_of(Keys),
            dplyr::all_of(new_variable)
          )

        names(old_tmp)[names(old_tmp) == old_variable] <- ".OldValue"
        names(new_tmp)[names(new_tmp) == new_variable] <- ".NewValue"

        comparison_tmp <- dplyr::inner_join(
          old_tmp,
          new_tmp,
          by = Keys
        )

        old_raw <- comparison_tmp$.OldValue
        new_raw <- comparison_tmp$.NewValue

        if (is.numeric(old_raw) && is.numeric(new_raw)) {
          same_value <- dplyr::case_when(
            is.na(old_raw) & is.na(new_raw) ~ TRUE,
            is.na(old_raw) | is.na(new_raw) ~ FALSE,
            abs(old_raw - new_raw) <= numeric_tolerance ~ TRUE,
            TRUE ~ FALSE
          )
        } else {
          old_chr <- as.character(old_raw)
          new_chr <- as.character(new_raw)

          same_value <- dplyr::case_when(
            is.na(old_chr) & is.na(new_chr) ~ TRUE,
            is.na(old_chr) | is.na(new_chr) ~ FALSE,
            old_chr == new_chr ~ TRUE,
            TRUE ~ FALSE
          )
        }

        idx <- !same_value

        if (!any(idx)) {
          return(NULL)
        }

        comparison_tmp[idx, Keys, drop = FALSE] %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            Variable = variable,
            OldVariable = old_variable,
            NewVariable = new_variable,
            OldValue = as.character(old_raw[idx]),
            NewValue = as.character(new_raw[idx]),
            OldClass = paste(class(old_raw), collapse = "/"),
            NewClass = paste(class(new_raw), collapse = "/"),
            ChangeType = dplyr::case_when(
              is.na(old_raw[idx]) & !is.na(new_raw[idx]) ~ "Old missing, New present",
              !is.na(old_raw[idx]) & is.na(new_raw[idx]) ~ "Old present, New missing",
              TRUE ~ "Different non-missing values"
            )
          )
      }
    )

  } else {

    modified_values <- tibble::tibble(
      !!!stats::setNames(
        rep(list(character()), length(Keys)),
        Keys
      ),
      Variable = character(),
      OldVariable = character(),
      NewVariable = character(),
      OldValue = character(),
      NewValue = character(),
      OldClass = character(),
      NewClass = character(),
      ChangeType = character()
    )
  }

  if (nrow(modified_values) == 0) {
    modified_values <- tibble::tibble(
      !!!stats::setNames(
        rep(list(character()), length(Keys)),
        Keys
      ),
      Variable = character(),
      OldVariable = character(),
      NewVariable = character(),
      OldValue = character(),
      NewValue = character(),
      OldClass = character(),
      NewClass = character(),
      ChangeType = character()
    )
  }

  # Summarize variable changes

  compared_record_count <- nrow(compared_keys)

  compared_cell_count <- compared_record_count * nrow(comparison_variables)

  modified_value_count <- nrow(modified_values)

  percent_values_modified <- ifelse(
    compared_cell_count > 0,
    round(100 * modified_value_count / compared_cell_count, 3),
    NA_real_
  )

  if (nrow(modified_values) > 0) {

    variable_change_summary <- modified_values %>%
      dplyr::count(
        Variable,
        name = "Changes"
      ) %>%
      dplyr::left_join(
        class_audit,
        by = "Variable"
      ) %>%
      dplyr::mutate(
        ComparedRecords = compared_record_count,
        PercentChanged = dplyr::if_else(
          ComparedRecords > 0,
          round(100 * Changes / ComparedRecords, 1),
          NA_real_
        )
      ) %>%
      dplyr::arrange(
        dplyr::desc(Changes),
        Variable
      )

  } else {

    variable_change_summary <- tibble::tibble(
      Variable = character(),
      Changes = integer(),
      OldVariable = character(),
      NewVariable = character(),
      OldClass = character(),
      NewClass = character(),
      ClassChanged = logical(),
      NameChanged = logical(),
      Compared = logical(),
      SkipReason = character(),
      ComparedRecords = integer(),
      PercentChanged = numeric()
    )
  }

  top_changed_variables <- variable_change_summary %>%
    dplyr::slice_head(n = 10)

  # Suspicious changes

  suspicious_change_threshold <- 10

  suspicious_from_change_rate <- variable_change_summary %>%
    dplyr::filter(
      PercentChanged >= suspicious_change_threshold
    ) %>%
    dplyr::mutate(
      Reason = paste0(
        "Changed in at least ",
        suspicious_change_threshold,
        "% of compared records."
      )
    )

  suspicious_from_class <- class_audit %>%
    dplyr::filter(
      Compared,
      ClassChanged
    ) %>%
    dplyr::mutate(
      Changes = NA_integer_,
      ComparedRecords = compared_record_count,
      PercentChanged = NA_real_,
      Reason = "Variable class changed between dataset versions."
    ) %>%
    dplyr::select(
      Variable,
      Changes,
      OldVariable,
      NewVariable,
      OldClass,
      NewClass,
      ClassChanged,
      NameChanged,
      Compared,
      SkipReason,
      ComparedRecords,
      PercentChanged,
      Reason
    )

  suspicious_changes <- dplyr::bind_rows(
    suspicious_from_change_rate,
    suspicious_from_class
  ) %>%
    dplyr::distinct(
      Variable,
      Reason,
      .keep_all = TRUE
    ) %>%
    dplyr::arrange(
      Variable,
      Reason
    )

  if (nrow(suspicious_changes) == 0) {
    suspicious_changes <- tibble::tibble(
      Variable = character(),
      Changes = integer(),
      OldVariable = character(),
      NewVariable = character(),
      OldClass = character(),
      NewClass = character(),
      ClassChanged = logical(),
      NameChanged = logical(),
      Compared = logical(),
      SkipReason = character(),
      ComparedRecords = integer(),
      PercentChanged = numeric(),
      Reason = character()
    )
  }

  # Build checks

  variables_changed_count <- if (modified_value_count > 0) {
    dplyr::n_distinct(modified_values$Variable)
  } else {
    0
  }

  added_variable_count <- nrow(added_variables)
  removed_variable_count <- nrow(removed_variables)
  added_record_count <- nrow(added_records)
  removed_record_count <- nrow(removed_records)
  duplicate_key_total <- duplicate_key_groups_old + duplicate_key_groups_new
  excluded_duplicate_key_count <- nrow(excluded_duplicate_keys)
  suspicious_change_count <- nrow(suspicious_changes)
  name_repair_change_count <- nrow(name_repair_audit)

  checks <- tibble::tibble(
    Check = c(
      "Key Types",
      "Duplicate Keys",
      "Records Added",
      "Records Removed",
      "Variables Added",
      "Variables Removed",
      "Name Repair Differences",
      "Variables Skipped From Cell Comparison",
      "Class Changes",
      "Values Modified",
      "Suspicious Changes",
      "Duplicate Keys Excluded From Cell Comparison"
    ),
    Count = c(
      key_type_count,
      duplicate_key_total,
      added_record_count,
      removed_record_count,
      added_variable_count,
      removed_variable_count,
      name_repair_change_count,
      skipped_variable_count,
      class_change_count,
      modified_value_count,
      suspicious_change_count,
      excluded_duplicate_key_count
    ),
    Status = dplyr::case_when(
      Check == "Key Types" & Count > 0 ~ "WARNING",
      Check == "Duplicate Keys" & Count > 0 ~ "FAIL",
      Check == "Records Added" & Count > 0 ~ "WARNING",
      Check == "Records Removed" & Count > 0 ~ "WARNING",
      Check == "Variables Added" & Count > 0 ~ "WARNING",
      Check == "Variables Removed" & Count > 0 ~ "WARNING",
      Check == "Name Repair Differences" & Count > 0 ~ "WARNING",
      Check == "Variables Skipped From Cell Comparison" & Count > 0 ~ "WARNING",
      Check == "Class Changes" & Count > 0 ~ "WARNING",
      Check == "Values Modified" & Count > 0 ~ "WARNING",
      Check == "Suspicious Changes" & Count > 0 ~ "WARNING",
      Check == "Duplicate Keys Excluded From Cell Comparison" & Count > 0 ~ "WARNING",
      TRUE ~ "PASS"
    ),
    Details = dplyr::case_when(
      Check == "Key Types" & Count > 0 ~
        "At least one key has different storage classes across dataset versions. Keys were coerced to character for comparison.",
      Check == "Key Types" ~
        "Key storage classes match across dataset versions.",

      Check == "Duplicate Keys" & Count > 0 ~
        "Duplicate key combinations were detected. Cell-level comparison only used keys unique in both datasets.",
      Check == "Duplicate Keys" ~
        "No duplicate key combinations detected.",

      Check == "Records Added" & Count > 0 ~
        "Some key combinations are present in NewData but not OldData.",
      Check == "Records Added" ~
        "No added records detected.",

      Check == "Records Removed" & Count > 0 ~
        "Some key combinations are present in OldData but not NewData.",
      Check == "Records Removed" ~
        "No removed records detected.",

      Check == "Variables Added" & Count > 0 ~
        "NewData contains variables not present in OldData after ignoring tibble-style name-repair suffixes.",
      Check == "Variables Added" ~
        "No added variables detected after ignoring tibble-style name-repair suffixes.",

      Check == "Variables Removed" & Count > 0 ~
        "OldData contains variables not present in NewData after ignoring tibble-style name-repair suffixes.",
      Check == "Variables Removed" ~
        "No removed variables detected after ignoring tibble-style name-repair suffixes.",

      Check == "Name Repair Differences" & Count > 0 ~
        "Some common variables have different raw names after tibble-style name repair.",
      Check == "Name Repair Differences" ~
        "No tibble-style name-repair differences detected among common variables.",

      Check == "Variables Skipped From Cell Comparison" & Count > 0 ~
        "Some variables were skipped because normalized names mapped ambiguously to multiple raw columns.",
      Check == "Variables Skipped From Cell Comparison" ~
        "No variables were skipped because of ambiguous normalized names.",

      Check == "Class Changes" & Count > 0 ~
        "At least one common variable changed class between versions.",
      Check == "Class Changes" ~
        "No class changes detected among compared variables.",

      Check == "Values Modified" & Count > 0 ~
        "At least one cell-level value changed among compared records.",
      Check == "Values Modified" ~
        "No cell-level value changes detected among compared records.",

      Check == "Suspicious Changes" & Count > 0 ~
        "High change-rate or class-change variables were detected.",
      Check == "Suspicious Changes" ~
        "No high change-rate or class-change variables detected.",

      Check == "Duplicate Keys Excluded From Cell Comparison" & Count > 0 ~
        "Some matching key combinations were excluded from cell comparison because they were duplicated in OldData or NewData.",
      Check == "Duplicate Keys Excluded From Cell Comparison" ~
        "No duplicate matching keys were excluded from cell comparison.",

      TRUE ~ ""
    )
  )

  # Build summary

  old_match_rate <- ifelse(
    nrow(old_keys) > 0,
    round(100 * nrow(matching_keys) / nrow(old_keys), 1),
    NA_real_
  )

  new_match_rate <- ifelse(
    nrow(new_keys) > 0,
    round(100 * nrow(matching_keys) / nrow(new_keys), 1),
    NA_real_
  )

  schema_change_count <- added_variable_count + removed_variable_count

  summary_tbl <- tibble::tibble(
    OldRows = nrow(OldData),
    NewRows = nrow(NewData),
    RowDifference = nrow(NewData) - nrow(OldData),
    OldColumns = ncol(OldData),
    NewColumns = ncol(NewData),
    ColumnDifference = ncol(NewData) - ncol(OldData),
    OldUniqueKeys = nrow(old_keys),
    NewUniqueKeys = nrow(new_keys),
    MatchingKeys = nrow(matching_keys),
    OldMatchRate = old_match_rate,
    NewMatchRate = new_match_rate,
    AddedRecords = added_record_count,
    RemovedRecords = removed_record_count,
    AddedVariables = added_variable_count,
    RemovedVariables = removed_variable_count,
    SchemaChangeCount = schema_change_count,
    NameRepairDifferences = name_repair_change_count,
    VariablesSkippedFromCellComparison = skipped_variable_count,
    CommonVariablesMapped = nrow(comparison_variable_map),
    CommonVariablesCompared = nrow(comparison_variables),
    ComparedKeys = nrow(compared_keys),
    ComparedCells = compared_cell_count,
    DuplicateKeyGroups_Old = duplicate_key_groups_old,
    DuplicateKeyGroups_New = duplicate_key_groups_new,
    ExcludedDuplicateKeys = excluded_duplicate_key_count,
    KeyTypeMismatches = key_type_count,
    ClassChanges = class_change_count,
    ModifiedValues = modified_value_count,
    PercentValuesModified = percent_values_modified,
    VariablesChanged = variables_changed_count,
    SuspiciousChanges = suspicious_change_count
  )

  # Build summary text

  added_variable_text <- if (added_variable_count > 0) {
    added_preview <- added_variables$Variable[
      seq_len(min(10, added_variable_count))
    ]

    if (added_variable_count > 10) {
      paste0(
        paste(added_preview, collapse = ", "),
        ", ... [",
        added_variable_count - 10,
        " more]"
      )
    } else {
      paste(added_preview, collapse = ", ")
    }
  } else {
    "None"
  }

  removed_variable_text <- if (removed_variable_count > 0) {
    removed_preview <- removed_variables$Variable[
      seq_len(min(10, removed_variable_count))
    ]

    if (removed_variable_count > 10) {
      paste0(
        paste(removed_preview, collapse = ", "),
        ", ... [",
        removed_variable_count - 10,
        " more]"
      )
    } else {
      paste(removed_preview, collapse = ", ")
    }
  } else {
    "None"
  }

  key_type_text <- if (key_type_count > 0) {
    key_types %>%
      dplyr::filter(
        OldType != NewType
      ) %>%
      dplyr::mutate(
        Text = paste0(
          Key,
          " (Old: ",
          OldType,
          "; New: ",
          NewType,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "; ")
  } else {
    "None"
  }

  most_changed_text <- if (nrow(top_changed_variables) > 0) {
    top_changed_variables %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          Changes,
          "; ",
          PercentChanged,
          "%)"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = ", ")
  } else {
    "None"
  }

  suspicious_text <- if (nrow(suspicious_changes) > 0) {
    suspicious_changes %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " [",
          Reason,
          "]"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "; ")
  } else {
    "None"
  }

  summary_text <- paste0(
    "DATASET COMPARISON SUMMARY\n\n",
    "Keys:\n",
    paste(Keys, collapse = ", "), "\n\n",

    "Rows:\n",
    "Old: ", nrow(OldData), "\n",
    "New: ", nrow(NewData), "\n",
    "Difference: ", nrow(NewData) - nrow(OldData), "\n\n",

    "Columns:\n",
    "Old: ", ncol(OldData), "\n",
    "New: ", ncol(NewData), "\n",
    "Difference: ", ncol(NewData) - ncol(OldData), "\n\n",

    "Unique Key Combinations:\n",
    "Old: ", nrow(old_keys), "\n",
    "New: ", nrow(new_keys), "\n",
    "Matching: ", nrow(matching_keys), "\n",
    "Old Match Rate: ", old_match_rate, "%\n",
    "New Match Rate: ", new_match_rate, "%\n\n",

    "Records:\n",
    "Added: ", added_record_count, "\n",
    "Removed: ", removed_record_count, "\n\n",

    "Variables:\n",
    "Added: ", added_variable_count, "\n",
    "Removed: ", removed_variable_count, "\n",
    "Schema Change Count: ", schema_change_count, "\n",
    "Name Repair Differences: ", name_repair_change_count, "\n",
    "Variables Skipped From Cell Comparison: ", skipped_variable_count, "\n",
    "Added Variables: ", added_variable_text, "\n",
    "Removed Variables: ", removed_variable_text, "\n\n",

    "Key Type Mismatches:\n",
    key_type_text, "\n\n",

    "Duplicate Key Groups:\n",
    "Old: ", duplicate_key_groups_old, "\n",
    "New: ", duplicate_key_groups_new, "\n",
    "Excluded Matching Keys: ", excluded_duplicate_key_count, "\n\n",

    "Cell-Level Changes:\n",
    "Compared Cells: ", compared_cell_count, "\n",
    "Modified Values: ", modified_value_count, "\n",
    "Percent Values Modified: ", percent_values_modified, "%\n",
    "Variables Changed: ", variables_changed_count, "\n",
    "Most Changed Variables: ", most_changed_text, "\n\n",

    "Suspicious Changes:\n",
    suspicious_text
  )

  # Return result

  result <- list(
    SummaryText = summary_text,
    Summary = summary_tbl,
    Fingerprint = fingerprint,
    KeyTypes = key_types,
    Checks = checks,
    StructureChanges = structure_changes,
    AddedVariables = added_variables,
    RemovedVariables = removed_variables,
    AddedRecords = added_records,
    RemovedRecords = removed_records,
    DuplicateKeys = duplicate_keys,
    ComparisonKeys = comparison_keys,
    NameRepairAudit = name_repair_audit,
    ComparisonVariableMap = comparison_variable_map,
    ClassAudit = class_audit,
    ModifiedValues = modified_values,
    VariableChangeSummary = variable_change_summary,
    TopChangedVariables = top_changed_variables,
    SuspiciousChanges = suspicious_changes
  )

  return(result)
}
