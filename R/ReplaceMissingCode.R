#' Replace Missing Codes with NA
#'
#' This function replaces specified missing codes in a data frame with `NA` values based on a given variable codebook.
#'
#' @param DataFrame A data frame containing the data.
#' @param VariableCodebook A data frame containing the variable codebook. It should have columns `Variable` and `MissingCode`, where `Variable` specifies the variable name in the data frame and `MissingCode` specifies the code to be replaced with `NA`.
#'
#' @return A data frame with specified missing codes replaced by `NA`.

#' @export
ReplaceMissingCode <- function(DataFrame, VariableCodebook) {
  CodeSubset <- VariableCodebook %>% dplyr::filter(!is.na(MissingCode))
  for (i in 1:nrow(CodeSubset)) {
    Var <- CodeSubset$Variable[i]
    MC <- CodeSubset$MissingCode[i]
    DataFrame[[Var]][DataFrame[[Var]] == MC] <- NA
  }
  return(DataFrame)
}
