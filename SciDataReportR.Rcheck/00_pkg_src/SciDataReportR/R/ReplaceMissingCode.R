#' Replace Missing Codes with NA
#'
#' This function replaces specified missing codes in a data frame with `NA` values based on a given variable codebook.
#'
#' @details
#' Study databases rarely leave missing values empty. They record *why* the
#' value is missing, using sentinel numbers outside the plausible range - `-7`
#' for "not applicable", `-8` for "refused", `-9` for "not asked", `999` for
#' "unknown". Left in place, those codes are read as real measurements: a mean
#' age climbs into the hundreds, and a variable coded `-7`/`-8` gets a negative
#' mean that no one notices until much later.
#'
#' This converts every code listed in the codebook to `NA`, so the values are
#' excluded from statistics and counted as missing.
#'
#' @section Several markers for one variable:
#' A variable usually has more than one missing code, since the whole point of
#' sentinels is to distinguish reasons. List them all in a single
#' `MissingCode` cell, separated by commas or semicolons:
#'
#' ```
#' Variable   MissingCode
#' age        -7, -8, -9
#' income     999; 9999
#' ```
#'
#' Whitespace around the separators is ignored, so `-7,-8` and `-7, -8` behave
#' identically. When a codebook is assembled in R rather than read from a
#' spreadsheet, `MissingCode` may also be a list-column of vectors.
#'
#' Codes are matched on the column's own scale: numerically for numeric
#' columns, and as text for character and factor columns, so markers like
#' `"Refused"` or `"Unknown"` work. When such a marker is a factor level, the
#' level is removed along with the values, so it does not linger as an empty
#' category in later tables.
#'
#' Variables whose `MissingCode` is blank are left untouched, and a codebook
#' naming variables that are not in `data` warns rather than failing - a
#' codebook usually describes more columns than any one analysis selects.
#'
#' @param data A data frame containing the data.
#' @param codebook A data frame containing the variable codebook. It must have
#'   columns `Variable` and `MissingCode`. `Variable` names the column in
#'   `data`; `MissingCode` holds the code, or several codes separated by commas
#'   or semicolons, to replace with `NA`. Rows with a missing or blank
#'   `MissingCode` are skipped.
#'
#' @return A data frame with specified missing codes replaced by `NA`.
#'
#' @seealso [CreateVariableTypesTemplate()], which generates a codebook with a
#'   `MissingCode` column ready to fill in, and [RevalueData()], which applies
#'   missing codes as part of the full relabelling workflow.
#'
#' @examples
#' # `age` uses three different codes to record three different reasons
#' df <- data.frame(
#'   id     = 1:6,
#'   age    = c(34, 999, 52, -7, 41, -8),
#'   score  = c(10, -9, -9, 15, 20, 12),
#'   status = c("Active", "Unknown", "Active", "Withdrawn", "Unknown", "Active")
#' )
#'
#' # Several markers for one variable go in one cell
#' codebook <- data.frame(
#'   Variable    = c("age", "score", "status"),
#'   MissingCode = c("999, -7, -8", "-9", "Unknown")
#' )
#'
#' # Before: the sentinels are averaged in as if they were ages
#' df
#' mean(df$age)
#'
#' # After: every listed code becomes NA
#' cleaned <- ReplaceMissingCode(df, codebook)
#' cleaned
#' mean(cleaned$age, na.rm = TRUE)
#' colSums(is.na(cleaned))
#'
#' # Blank codes are skipped
#' ReplaceMissingCode(
#'   df,
#'   data.frame(Variable = c("id", "age"), MissingCode = c(NA, "999, -7, -8"))
#' )
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param VariableCodebook \strong{Deprecated} (since 19.15.0). Use \code{codebook} instead.
#' @export
ReplaceMissingCode <- function(data,
    codebook,
    DataFrame = lifecycle::deprecated(),
    VariableCodebook = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "ReplaceMissingCode(DataFrame)", "ReplaceMissingCode(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(VariableCodebook)) {
    lifecycle::deprecate_warn("19.15.0", "ReplaceMissingCode(VariableCodebook)", "ReplaceMissingCode(codebook)")
    codebook <- VariableCodebook
  }
  if (!missing(codebook)) VariableCodebook <- codebook

  if (!is.data.frame(VariableCodebook)) {
    stop("`codebook` must be a data frame.", call. = FALSE)
  }
  for (required in c("Variable", "MissingCode")) {
    if (!required %in% names(VariableCodebook)) {
      stop(sprintf("`codebook` must contain a `%s` column.", required), call. = FALSE)
    }
  }

  # Filter out rows where MissingCode is NA
  CodeSubset <- VariableCodebook %>% dplyr::filter(!is.na(MissingCode))

  # One cell can hold several markers, separated by commas or semicolons, and
  # a list-column of codes is accepted for codebooks assembled in R.
  ParseMissingCodes <- function(value) {
    if (is.list(value)) value <- unlist(value, use.names = FALSE)
    value <- as.character(value)
    value <- unlist(strsplit(value, "[,;]"), use.names = FALSE)
    value <- trimws(value)
    value[!is.na(value) & nzchar(value)]
  }

  UnmatchedVariables <- character(0)

  # seq_len, not 1:nrow: a codebook where nothing carries a missing code is a
  # no-op, not an error.
  for (i in seq_len(nrow(CodeSubset))) {
    Var <- as.character(CodeSubset$Variable[i])
    MissingCodes <- ParseMissingCodes(CodeSubset$MissingCode[[i]])

    if (length(MissingCodes) == 0) next

    # A codebook routinely outlives the column selection it was written for.
    if (!Var %in% names(DataFrame)) {
      UnmatchedVariables <- c(UnmatchedVariables, Var)
      next
    }

    Column <- DataFrame[[Var]]

    for (MC in MissingCodes) {
      if (is.numeric(Column)) {
        # Codes are compared on the column's own scale; a non-numeric code
        # cannot match a numeric column and is left alone.
        Numeric <- suppressWarnings(as.numeric(MC))
        Hits <- if (is.na(Numeric)) FALSE else !is.na(Column) & Column == Numeric
      } else {
        Hits <- !is.na(Column) & as.character(Column) == MC
      }

      if (any(Hits)) Column[Hits] <- NA

      # A factor level that means "missing" should not survive as an empty
      # level in downstream tables.
      if (is.factor(Column) && MC %in% levels(Column)) {
        Column <- factor(Column, levels = setdiff(levels(Column), MC))
      }
    }

    DataFrame[[Var]] <- Column
  }

  if (length(UnmatchedVariables) > 0) {
    warning(
      "Missing codes were listed for variable(s) not present in `data`: ",
      paste(unique(UnmatchedVariables), collapse = ", "),
      call. = FALSE
    )
  }

  return(DataFrame)
}
