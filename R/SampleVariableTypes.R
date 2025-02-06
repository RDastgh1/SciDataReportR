#' Example Dataset: SampleVariableTypes
#'
#' An example of a modified VariableTypes file to be used to Revalue SampleData
#'
#' @format A data frame with 11 columns and 138 rows
#' \describe{
#'   \item{Variable}{Variable name, as listed in column name of data file}
#'   \item{Label}{Variable Label, as desired in tables and figures}
#'   \item{Type}{How the variable should be treated, e.g. continuous(double) or categorical}
#'   \item{Category}{Optional for future filtering: what category type this variable belongs to}
#'   \item{Recode}{0 or 1, whether or not this variable should be recoded}
#'   \item{Code}{If the variable should be recoded, how it should be recoded}
#'   \item{Notes}{Optional}
#'   \item{Exclude}{Optional for filtering}
#'   \item{Subcategory}{Optional additional categories}
#'   \item{Include}{optional for filtering}
#'   \item{Missingcode}{optional: if a value should be considered NA}
#' }
#' @source Exported from CreateVariableTypes(SampleData, "SampleData.csv") and modified in excel
#' @examples
#' data(SampleVariableTypes)
#' head(SampleVariableTypes)
"SampleVariableTypes"
