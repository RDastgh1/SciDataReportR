#' Create a Template for Variable Types
#'
#' This function generates a data frame that summarizes the types and labels of the variables in a given data frame.
#' It optionally saves this summary to a CSV file.
#'
#' @param DataFrame A data frame containing the variables to be summarized.
#' @param CSVFileName A string specifying the path and name of the CSV file to save the summary.
#'                    If NULL (the default), the CSV file will not be created.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Variable}{The names of the variables in the input data frame.}
#'   \item{Label}{The labels of the variables, if available; otherwise, the variable names.}
#'   \item{Type}{The data types of the variables, converted to more user-friendly descriptions.}
#'   \item{Category}{A placeholder column for categorizing variables (default is NA).}
#'   \item{Recode}{A placeholder column for recoding information (default is NA).}
#'   \item{Code}{A placeholder column for code information (default is NA).}
#'   \item{Notes}{A placeholder column for any additional notes (default is an empty string).}
#'   \item{Exclude}{A placeholder column for exclusion flags (default is NA).}
#' }
#'
#' @examples
#' df <- data.frame(
#'   num = c(1.1, 2.2),
#'   int = c(1L, 2L),
#'   fact = factor(c("A", "B")),
#'   char = c("a", "b"),
#'   date = as.Date(c("2021-01-01", "2021-01-02"))
#' )
#' CreateVariableTypesTemplate(df)
#' CreateVariableTypesTemplate(df, "variable_types.csv")
#'
#' @importFrom sjlabelled get_label
#' @importFrom utils write.csv
#' @export
CreateVariableTypesTemplate <- function(DataFrame, CSVFileName = NULL) {

  # Get the classes of the variables in the DataFrame
  Types <- sapply(DataFrame, class)

  # Convert to a more user-friendly factor
  Types <- factor(Types, levels = c("numeric", "integer", "factor", "character", "Date", "logical"),
                  labels = c("Double", "Double", "Categorical", "String", "Date", "Categorical"))

  # Get labels of the variables, use column names as default
  DataLabels <- sjlabelled::get_label(DataFrame, def.value = colnames(DataFrame))

  # Create the VariableTypes data frame
  VariableTypes <- data.frame(Variable = colnames(DataFrame),
                              Label = DataLabels,
                              Type = Types,
                              Category = NA,
                              Recode = NA,
                              Code = NA,
                              Notes = "",
                              Exclude = NA,
                              MissingCode = "")

  # If a CSV file name is provided, write the data frame to a CSV file
  if (!is.null(CSVFileName)) {
    write.csv(VariableTypes, CSVFileName, row.names = FALSE)
  }

  return(VariableTypes)
}
