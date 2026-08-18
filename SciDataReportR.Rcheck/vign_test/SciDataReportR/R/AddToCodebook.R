#' Add a new variable to a codebook
#'
#' Adds one variable entry to a codebook while preserving its existing schema.
#' In addition to the standard codebook fields, named values supplied through
#' `...` populate user-defined columns. A new named `...` column is added to
#' the codebook (with `NA` for existing rows) and produces a warning.
#'
#' The function uses the supplied codebook as the source of truth for coded
#' metadata. For every supplied field except `Variable`, `Label`, and `Notes`, it warns when
#' a non-missing value has not previously appeared in that column or has a
#' different storage type. These warnings are advisory: the row is still added
#' and R may promote the column type while binding the new row.
#'
#' @param codebook A data frame representing the codebook. It must contain a
#'   `Variable` column.
#' @param VariableName A single, non-missing, non-empty character variable name.
#'   It must not already appear in `codebook$Variable`.
#' @param VariableLabel A single label for the variable. Defaults to
#'   `VariableName` when `NA`.
#' @param VariableType A single value for the `Type` column.
#' @param VariableCategory A single value for the `Category` column.
#' @param VariableRecode A single value for the `Recode` column.
#' @param VariableCode A single value for the `Code` column.
#' @param VariableExclude A single value for the `Exclude` column.
#' @param VariableNotes A single value for the `Notes` column.
#' @param ... Named, single atomic values for user-defined codebook columns.
#'   Names matching existing columns populate them. New names create a column
#'   and warn. Standard fields (`Variable`, `Label`, `Type`, `Category`,
#'   `Recode`, `Code`, `Exclude`, and `Notes`) must be supplied through their
#'   corresponding formal arguments.
#'
#' @return A data frame representing the updated codebook with the new variable
#'   added.
#'
#' @examples
#' \donttest{
#' # An existing codebook, including the user-defined `Domain` column
#' codebook <- data.frame(
#'   Variable = c("sex", "visit", "mmse"),
#'   Label = c("Sex assigned at birth", "Study visit", "MMSE total score"),
#'   Type = c("Categorical", "Categorical", "Double"),
#'   Category = c("Demographics", "Design", "Cognition"),
#'   Recode = c(1, 0, 0),
#'   Code = c("0 = Female; 1 = Male", "1 = Baseline; 2 = Follow-up", NA),
#'   Exclude = c(FALSE, FALSE, FALSE),
#'   Notes = NA_character_,
#'   Domain = c("Clinical", "Study", "Clinical")
#' )
#'
#' # Each call returns the updated codebook, so entries chain
#' codebook <- AddToCodebook(
#'   codebook, "age", "Age at enrollment", "Double", "Demographics",
#'   Domain = "Clinical"
#' )
#'
#' # `Source` is a new column: created with a warning, back-filled with NA
#' codebook <- AddToCodebook(
#'   codebook, "site", "Enrolling site", "Categorical", "Design",
#'   Source = "REDCap"
#' )
#'
#' # Kept, but warned about: they conflict with the established schema
#' codebook <- AddToCodebook(
#'   codebook, "participant_id", "Participant identifier",
#'   VariableRecode = 4, VariableExclude = 3
#' )
#'
#' # The result as a formatted table
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(codebook, height = "320px", full_width = TRUE)
#' )))
#' }
#'
#' @param CB \strong{Deprecated} (since 19.15.0). Use \code{codebook} instead.
#' @export
AddToCodebook <- function(codebook,
    VariableName,
    VariableLabel = NA,
    VariableType = NA,
    VariableCategory = NA,
    VariableRecode = NA,
    VariableCode = NA,
    VariableExclude = NA,
    VariableNotes = NA,
    CB = lifecycle::deprecated(),
    ...) {
  VariableTypeSupplied <- !missing(VariableType)
  VariableCategorySupplied <- !missing(VariableCategory)
  VariableRecodeSupplied <- !missing(VariableRecode)
  VariableCodeSupplied <- !missing(VariableCode)
  VariableExcludeSupplied <- !missing(VariableExclude)
  VariableNotesSupplied <- !missing(VariableNotes)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(CB)) {
    lifecycle::deprecate_warn("19.15.0", "AddToCodebook(CB)", "AddToCodebook(codebook)")
    codebook <- CB
  }

  if (missing(codebook)) {
    stop("`codebook` must be supplied.", call. = FALSE)
  }
  if (!is.data.frame(codebook)) {
    stop("`codebook` must be a data frame.", call. = FALSE)
  }
  if (!"Variable" %in% names(codebook)) {
    stop("`codebook` must contain a `Variable` column.", call. = FALSE)
  }

  ValidateCodebookValue <- function(value, argument) {
    if (is.null(value) || !is.atomic(value) || length(value) != 1) {
      stop(sprintf("`%s` must be a single atomic value or NA.", argument), call. = FALSE)
    }
  }

  ValidateCodebookValue(VariableName, "VariableName")
  if (!is.character(VariableName) || is.na(VariableName) || !nzchar(trimws(VariableName))) {
    stop("`VariableName` must be a single, non-missing, non-empty character string.", call. = FALSE)
  }
  if (VariableName %in% as.character(codebook$Variable)) {
    stop(sprintf("`VariableName` already exists in `codebook$Variable`: %s", VariableName), call. = FALSE)
  }

  legacy_values <- list(
    Variable = VariableName,
    Label = VariableLabel,
    Type = VariableType,
    Category = VariableCategory,
    Recode = VariableRecode,
    Code = VariableCode,
    Exclude = VariableExclude,
    Notes = VariableNotes
  )
  legacy_supplied <- c(
    Variable = TRUE,
    Label = TRUE,
    Type = VariableTypeSupplied,
    Category = VariableCategorySupplied,
    Recode = VariableRecodeSupplied,
    Code = VariableCodeSupplied,
    Exclude = VariableExcludeSupplied,
    Notes = VariableNotesSupplied
  )

  for (field in names(legacy_values)) {
    ValidateCodebookValue(legacy_values[[field]], paste0("Variable", field))
  }
  if (is.na(VariableLabel)) {
    legacy_values$Label <- VariableName
  }

  custom_values <- list(...)
  custom_names <- names(custom_values)
  if (length(custom_values) > 0) {
    if (is.null(custom_names) || any(is.na(custom_names) | !nzchar(custom_names))) {
      stop("All `...` values must be named.", call. = FALSE)
    }
    if (anyDuplicated(custom_names)) {
      stop("`...` column names must be unique.", call. = FALSE)
    }
    standard_names <- names(legacy_values)
    duplicated_standard <- intersect(custom_names, standard_names)
    if (length(duplicated_standard) > 0) {
      stop(
        sprintf(
          "Supply standard field(s) through their formal argument(s), not `...`: %s",
          paste(duplicated_standard, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    for (field in custom_names) {
      ValidateCodebookValue(custom_values[[field]], field)
    }
  }

  original_columns <- names(codebook)
  legacy_new_columns <- names(legacy_values)[legacy_supplied & !names(legacy_values) %in% original_columns]
  custom_new_columns <- setdiff(custom_names, original_columns)

  for (field in custom_new_columns) {
    warning(
      sprintf("`%s` is not an existing codebook column; adding it with NA for existing rows.", field),
      call. = FALSE
    )
  }

  # Add requested legacy fields and user-defined fields after all old columns.
  for (field in c(legacy_new_columns, custom_new_columns)) {
    codebook[[field]] <- NA
  }

  WarnIfSchemaMismatch <- function(field, value) {
    if (field %in% c("Variable", "Label", "Notes") || is.na(value) || !field %in% original_columns) {
      return(invisible(NULL))
    }

    existing <- codebook[[field]][seq_len(nrow(codebook))]
    non_missing <- !is.na(existing)
    if (!any(non_missing)) {
      return(invisible(NULL))
    }
    existing <- existing[non_missing]

    value_seen <- any(vapply(seq_along(existing), function(i) {
      identical(value, existing[i])
    }, logical(1)))
    type_matches <- identical(typeof(value), typeof(existing)) &&
      identical(class(value), class(existing))

    diagnostics <- character(0)
    if (!value_seen) {
      diagnostics <- c(diagnostics, "has not previously appeared in this column")
    }
    if (!type_matches) {
      diagnostics <- c(diagnostics, "has a different storage type or class from existing values")
    }
    if (length(diagnostics) > 0) {
      warning(
        sprintf("`%s` %s.", field, paste(diagnostics, collapse = " and ")),
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  supplied_values <- legacy_values[legacy_supplied]
  supplied_values <- c(supplied_values, custom_values)
  for (field in names(supplied_values)) {
    WarnIfSchemaMismatch(field, supplied_values[[field]])
  }

  all_columns <- names(codebook)
  row_values <- stats::setNames(as.list(rep(NA, length(all_columns))), all_columns)
  row_values[names(supplied_values)] <- supplied_values
  NewRow <- as.data.frame(row_values, stringsAsFactors = FALSE)

  plyr::rbind.fill(codebook, NewRow)
}
