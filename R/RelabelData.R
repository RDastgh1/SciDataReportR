ReLabelData<- function(Data, Vars, VariableTypes){
  #' @importFrom labelled var_label

   for (v in Vars){
    if(v %in% VariableTypes$Variable){
      var_label(Data[[v]])<- VariableTypes$Label[VariableTypes$Variable %in% v]
    }
  }
  return(Data)
}
