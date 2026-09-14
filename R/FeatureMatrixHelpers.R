#' Validate a wide feature matrix
#'
#' @param data A data frame with samples in rows.
#' @param variables Character vector of numeric feature columns.
#'
#' @return Invisibly returns the validated feature names.
#' @keywords internal
#' @noRd
ValidateFeatureMatrix <- function(data, variables = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (is.null(variables)) {
    variables <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  if (!is.character(variables) || length(variables) == 0) {
    stop("`variables` must identify at least one numeric feature column.", call. = FALSE)
  }
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("These `variables` are not columns in `data`: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  non_numeric <- variables[!vapply(data[variables], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("These `variables` must be numeric: ", paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  invisible(variables)
}
