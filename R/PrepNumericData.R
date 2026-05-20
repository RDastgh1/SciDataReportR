#' Prepare numeric data safely for analysis
#'
#' Safely coerces selected variables in a data frame to numeric format while
#' preserving column names and replacing non-finite values (`Inf`, `-Inf`,
#' `NaN`) with `NA`.
#'
#' This function is useful for preparing real-world biomedical, omics, and
#' survey datasets for downstream statistical analyses such as correlations,
#' PCA, clustering, regression, and factor analysis.
#'
#' Character and factor variables are converted using
#' `as.numeric(as.character(x))` to avoid incorrect factor level coercion.
#'
#' @param df A data frame.
#' @param variables Character vector of variable names to process.
#'   Defaults to all columns.
#'
#' @return A data frame with selected variables converted to numeric and
#'   non-finite values replaced with `NA`.
#'
#' @examples
#' df <- data.frame(
#'   x = c("1", "2", "3"),
#'   y = c(1, 2, Inf),
#'   z = factor(c("4", "5", "6"))
#' )
#'
#' PrepNumericData(df)
#'
#' @export
PrepNumericData <- function(df, variables = names(df)) {

  # ---------------------------
  # Validate inputs
  # ---------------------------
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }

  variables <- intersect(as.character(variables), names(df))

  if (length(variables) == 0) {
    stop("No valid variables found in df.")
  }


# Prepare data

  out <- df[, variables, drop = FALSE]
  names_out <- names(out)

  out <- as.data.frame(
    lapply(out, function(z) {

      # Convert non-numeric variables safely
      if (!is.numeric(z)) {
        z <- suppressWarnings(
          as.numeric(as.character(z))
        )
      }

      # Replace Inf, -Inf, NaN with NA
      z[!is.finite(z)] <- NA_real_

      z
    }),
    check.names = FALSE
  )

  names(out) <- names_out

  out
}
