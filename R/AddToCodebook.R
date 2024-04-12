AddToCodebook <- function(CB, VariableName, VariableLabel = NA, VariableType = NA, VariableCategory = NA, VariableRecode = NA, VariableCode = NA, VariableExclude = NA, VariableNotes = NA){
  if(is.na(VariableLabel)){
    VariableLabel = VariableName
  }

  NewRow = data.frame(Variable = VariableName, Label = VariableLabel, Type = VariableType, Category = VariableCategory,
                      Recode = VariableRecode, Code = VariableCode, Exclude = VariableExclude, Notes = "VariableNotes")

  CB <- plyr::rbind.fill(CB, NewRow)

  return(CB)
}
