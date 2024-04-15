#' Add a new variable to a codebook
#'
#' This function adds a new variable entry to an existing codebook.
#'
#' @param CB A data frame representing the codebook.
#' @param VariableName The name of the variable.
#' @param VariableLabel The label of the variable (default: same as VariableName).
#' @param VariableType The type of the variable (e.g., numeric, categorical).
#' @param VariableCategory The category to which the variable belongs.
#' @param VariableRecode The recode information for the variable.
#' @param VariableCode The code associated with the variable.
#' @param VariableExclude A flag indicating whether the variable should be excluded.
#' @param VariableNotes Any additional notes or comments about the variable.
#'
#' @return A data frame representing the updated codebook with the new variable added.
#'
#' @examples
#' # Create an empty codebook
#' codebook <- data.frame(Variable = character(0), Label = character(0),
#'                        Type = character(0), Category = character(0),
#'                        Recode = character(0), Code = character(0),
#'                        Exclude = logical(0), Notes = character(0))
#'
#' # Add a new variable to the codebook
#' codebook <- AddToCodebook(codebook, "Age", "Age of participants", "numeric", "Demographics")
#'
#' @export
AddToCodebook <- function(CB, VariableName, VariableLabel = NA, VariableType = NA, VariableCategory = NA, VariableRecode = NA, VariableCode = NA, VariableExclude = NA, VariableNotes = NA) {

  # Set default variable label if NA
  if (is.na(VariableLabel)) {
    VariableLabel <- VariableName
  }

  # Create a new row with variable information
  NewRow <- data.frame(Variable = VariableName, Label = VariableLabel, Type = VariableType, Category = VariableCategory,
                       Recode = VariableRecode, Code = VariableCode, Exclude = VariableExclude, Notes = VariableNotes)

  # Append the new row to the codebook
  CB <- plyr::rbind.fill(CB, NewRow)

  return(CB)
}
