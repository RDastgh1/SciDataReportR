#' Combine Two Codebooks with Conflict Detection
#'
#' @description
#' `CombineCodebooks` compares two codebook data frames (an old version and a new version), identifies added or removed variables and columns, detects cell-by-cell differences, and produces a combined codebook.  For records without any differences, it returns a single "Combined" row; for records with differences, it returns both the "Old" and "New" rows, flagged as conflicts.
#'
#' @param OldCodebook A data.frame or tibble representing the old codebook.  Each row must correspond to a single variable entry.
#' @param NewCodebook A data.frame or tibble representing the new codebook.  Must have the same structure (column names) as `OldCodebook`, though extra or missing columns will be handled.
#' @param key A string giving the name of the key column to identify variables (e.g. "Variable").  Defaults to "Variable".
#'
#' @details
#' This function:
#' * Builds a `combined_df` that includes all columns from both inputs:
#'   - For non-conflicting records, a single row with `Version = "Combined"`.
#'   - For conflicting records, two rows: one with `Version = "Old"` and one with `Version = "New"`.
#' * Flags each row with `CHECKFORCONFLICTS` (0 for combined, 1 for old/new conflict rows).
#'
#' @return A list with elements:
#' * `added_variables`: character vector of keys present only in `NewCodebook`.
#' * `removed_variables`: character vector of keys present only in `OldCodebook`.
#' * `columns_added`: character vector of column names present only in `NewCodebook`.
#' * `columns_removed`: character vector of column names present only in `OldCodebook`.
#' * `value_differences`: tibble of cell-level differences (`RowID`, `Field`, `OldValue`, `NewValue`).
#' * `combined_df`: tibble containing the merged codebook with versions and conflict flag.
#' @import dplyr
#' @import tidyr
#' @export
CombineCodebooks <- function(OldCodebook, NewCodebook, key = "Variable") {
  # Internal implementation as provided...
  # 1) Tag rows so we can merge 1:1 even if 'key' repeats
  df_old <- dplyr::mutate(OldCodebook, RowID = seq_len(nrow(OldCodebook)))
  df_new <- dplyr::mutate(NewCodebook, RowID = seq_len(nrow(NewCodebook)))

  # 2) Added / removed variables & columns
  added_variables   <- setdiff(df_new[[key]], df_old[[key]])
  removed_variables <- setdiff(df_old[[key]], df_new[[key]])
  cols_old          <- names(df_old)
  cols_new          <- names(df_new)
  columns_added     <- setdiff(cols_new, cols_old)
  columns_removed   <- setdiff(cols_old, cols_new)

  # 3) Common fields for diffing
  common_fields <- intersect(cols_old, cols_new) %>% setdiff(c("RowID", key))

  # 4) Cell-by-cell diffs on common fields
  df_old_long <- df_old %>%
    mutate(dplyr::across(dplyr::all_of(common_fields), as.character)) %>%
    tidyr::pivot_longer(dplyr::all_of(common_fields), names_to = "Field", values_to = "OldValue") %>%
    dplyr::select(RowID, Field, OldValue)

  df_new_long <- df_new %>%
    mutate(dplyr::across(dplyr::all_of(common_fields), as.character)) %>%
    tidyr::pivot_longer(dplyr::all_of(common_fields), names_to = "Field", values_to = "NewValue") %>%
    dplyr::select(RowID, Field, NewValue)

  merged_vals <- dplyr::full_join(df_old_long, df_new_long, by = c("RowID", "Field"))
  value_differences <- merged_vals %>%
    dplyr::filter(
      OldValue != NewValue |
        (is.na(OldValue) & !is.na(NewValue)) |
        (!is.na(OldValue) & is.na(NewValue))
    ) %>%
    dplyr::arrange(RowID, Field)

  # 5) Determine conflicting rows
  conflict_rids <- unique(value_differences$RowID)

  # 6) Prepare superset of all fields
  all_fields <- setdiff(union(cols_old, cols_new), "RowID")

  # 7) Pad dataframes to same structure
  missing_old <- setdiff(all_fields, cols_old)
  if (length(missing_old)) df_old[missing_old] <- NA
  missing_new <- setdiff(all_fields, cols_new)
  if (length(missing_new)) df_new[missing_new] <- NA

  # 8) Merge with suffixes
  merged_join <- dplyr::full_join(df_old, df_new, by = "RowID", suffix = c(".old", ".new"))

  # 9) Coalesce all fields
  fields_to_coalesce <- union(key, all_fields)
  for (fld in fields_to_coalesce) {
    merged_join[[fld]] <- dplyr::coalesce(
      merged_join[[paste0(fld, ".new")]],
      merged_join[[paste0(fld, ".old")]]
    )
  }

  # 10a) Non-conflicts: single Combined row
  no_conflicts <- merged_join %>%
    dplyr::filter(!RowID %in% conflict_rids) %>%
    dplyr::select(dplyr::all_of(fields_to_coalesce)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Version = "Combined")

  # 10b) Conflicts: Old and New rows
  old_conflicts <- df_old %>%
    dplyr::filter(RowID %in% conflict_rids) %>%
    dplyr::select(dplyr::all_of(fields_to_coalesce)) %>%
    dplyr::mutate(Version = "Old")

  new_conflicts <- df_new %>%
    dplyr::filter(RowID %in% conflict_rids) %>%
    dplyr::select(dplyr::all_of(fields_to_coalesce)) %>%
    dplyr::mutate(Version = "New")

  # 11) Combine and flag
  combined_df <- dplyr::bind_rows(no_conflicts, old_conflicts, new_conflicts) %>%
    dplyr::arrange(.data[[key]], Version) %>%
    dplyr::mutate(CHECKFORCONFLICTS = if_else(Version == "Combined", 0L, 1L))

  # 12) Return results
  list(
    added_variables   = added_variables,
    removed_variables = removed_variables,
    columns_added     = columns_added,
    columns_removed   = columns_removed,
    value_differences = value_differences,
    combined_df       = combined_df
  )
}
