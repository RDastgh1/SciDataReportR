#' Create a Template for Variable Types
#'
#' This function generates a data frame that summarizes the types and labels of the variables in a given data frame.
#' It optionally saves this summary to a CSV file.
#'
#' @param DataFrame A data frame containing the variables to be summarized.
#' @param CSVFileName A string specifying the path and name of the CSV file to save the summary.
#'                    If NULL (the default), the CSV file will not be created.
#' @param GuessCategorical A logical variable specifying if the function should guess what variables are categorical based on having <= 5 unique values
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
CreateVariableTypesTemplate <- function(DataFrame, CSVFileName = NULL, GuessCategorical = TRUE) {

  # Get the classes of the variables in the DataFrame
  Types <- sapply(DataFrame, class)
  Types <- factor(Types, levels = c("numeric", "integer", "factor",
                                    "character", "Date", "logical"), labels = c("Double",
                                                                                "Double", "Categorical", "String", "Date", "Categorical"))

  # Decide if categorical Values need to be calculated based on whether there were any categorical variables or not.
  if(GuessCategorical == T || (is.null(GuessCategorical) & sum(Types == "Categorical", na.rm = T) == 0)){

    # Get the number of unique values for each. If it's <=5, decide that it is categorical
    unique_counts <- sapply(DataFrame, function(x) length(unique(x)))

    Types[unique_counts <= 5] <- "Categorical"
  }

  # Get variable labels (if available)
  DataLabels <- sjlabelled::get_label(DataFrame, def.value = colnames(DataFrame))

  # Create the template dataframe
  VariableTypes <- data.frame(Variable = colnames(DataFrame),
                              Label = DataLabels, Type = Types, Category = NA, Recode = NA,
                              Code = NA, Notes = "", Exclude = NA, MissingCode = "")

  # Handle labelled factors (e.g., variables loaded from SPSS using `haven`)
  labelled_factors <- sapply(DataFrame, sjlabelled::is_labelled)

  for (i in which(labelled_factors)) {
    var_name <- colnames(DataFrame)[i]

    # Set Recode to 1 for labelled factors
    VariableTypes$Recode[i] <- 1

    # Get the labels for the factor and create a recoding string
    labels <- levels(as.factor(DataFrame[[var_name]]))
    codes <- sjlabelled::get_labels(DataFrame[[var_name]])

    recode_str <- paste0(paste(labels, codes, sep = "=",
                               collapse = "; "))
    # Set the Code column to show recoding
    VariableTypes$Code[i] <- recode_str
  }

  # Optionally save the output to a CSV file
  if (!is.null(CSVFileName)) {
    write.csv(VariableTypes, CSVFileName, row.names = FALSE)
  }

  return(VariableTypes)
}
