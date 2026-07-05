#' Replace Missing Codes with NA
#'
#' This function replaces specified missing codes in a data frame with `NA` values based on a given variable codebook.
#'
#' @param data A data frame containing the data.
#' @param codebook A data frame containing the variable codebook. It should have columns `Variable` and `MissingCode`, where `Variable` specifies the variable name in the data frame and `MissingCode` specifies the code to be replaced with `NA`.
#'
#' @return A data frame with specified missing codes replaced by `NA`.
#'
#' @examples
#' # A data frame that uses sentinel values to encode missingness
#' df <- data.frame(
#'   id    = 1:6,
#'   age   = c(34, 999, 52, 999, 41, 29),
#'   score = c(10, -9, -9, 15, 20, 12)
#' )
#'
#' # A codebook mapping each variable to its missing code(s)
#' codebook <- data.frame(
#'   Variable    = c("age", "score"),
#'   MissingCode = c("999", "-9")
#' )
#'
#' # Before: sentinel codes still present
#' df
#'
#' # After: sentinel codes replaced with NA
#' ReplaceMissingCode(df, codebook)
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

  # Filter out rows where MissingCode is NA
  CodeSubset <- VariableCodebook %>% dplyr::filter(!is.na(MissingCode))

  # Loop through each row of the codebook
  for (i in 1:nrow(CodeSubset)) {
    Var <- CodeSubset$Variable[i]  # The variable name
    MissingCodes <- strsplit(CodeSubset$MissingCode[i], ",\\s*")[[1]]  # Split by commas and handle spaces

    # Loop through each code in MissingCodes and replace with NA
    for (MC in MissingCodes) {
      MC <- as.numeric(MC)  # Convert to numeric in case the missing code is represented as a number
      DataFrame[[Var]][DataFrame[[Var]] == MC] <- NA
    }
  }

  return(DataFrame)
}
