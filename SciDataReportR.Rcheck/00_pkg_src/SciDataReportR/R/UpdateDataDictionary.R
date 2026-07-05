#' Update an existing data dictionary with new variables and types
#'
#' This function updates an existing data dictionary with new variables and their corresponding types.
#'
#' @param OldDataDictionary The existing data dictionary to be updated.
#' @param NewDataFrame The new data frame containing additional variables.
#'
#' @return A list containing the new variables and the updated data dictionary.
#'
#' @export
UpdateDataDictionary <- function(OldDataDictionary, NewDataFrame) {

  # Identify new variables not in the existing data dictionary
  NewDataVars <- colnames(NewDataFrame)
  VariableTypesVars <- OldDataDictionary$Variable
  NewVars <- NewDataVars[!(NewDataVars %in% VariableTypesVars)]

  # Create a data frame for new variables
  NewVariableTypes <- data.frame(Variable = NewVars, Label = NewVars, Type = NA, Category = NA, Recode = NA)

  # Fill in types for each new variable
  for (i in seq_along(NewVars)) {
    Var <- NewVariableTypes$Variable[i]
    VarType <- class(NewDataFrame[[Var]])
    if (VarType == "numeric") {
      NewVariableTypes$Type[i] <- "Double"
    } else if (VarType == "character") {
      NewVariableTypes$Type[i] <- "Categorical"
    }
  }

  # Combine existing and new variable types
  CombinedVariableTypes <- dplyr::full_join(OldDataDictionary, NewVariableTypes, by = "Variable")

  # Return the new additions and the updated data dictionary
  return(list(NewAdditions = NewVariableTypes, NewDataDictionary = CombinedVariableTypes))
}
