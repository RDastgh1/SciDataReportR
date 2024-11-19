#' Revalue Data
#'
#'
#' This function revalues variables in a dataset according to the specifications provided in a codebook.
#'
#' @param DatatoRevalue The dataset to be revalued.
#' @param VarTypes A data frame containing information about the variables and how they should be revalued.
#'                 It should have columns: Variable (variable names), Recode (yes/no for recoding),
#'                 Code (the revalue codes separated by ","), Type (Categorical, Double, or Ordinal), and Label (optional label for the variable).
#' @param missingVal The value to use for missing values (default is -999).
#' @param splitchar The character used to split codes in the Code column (default is ";").
#' @return A list containing the revalued dataset, a warning list, and a list of recoded variables.
#' @export
RevalueData <- function(DatatoRevalue, VarTypes, missingVal = -999, splitchar = ";") {
  library(sjlabelled)

  # Remove NA variables from VarTypes
  VarTypes <- VarTypes %>% filter(!is.na(Variable))

  # Initialize the revalued data with the original data
  RevaluedData <- DatatoRevalue

  # Get variables in data and codebook
  VariablesinData <- colnames(DatatoRevalue)
  VariablesinVarTypes <- VarTypes$Variable
  vars <- VariablesinData[VariablesinData %in% VariablesinVarTypes]
  variablesinData_notinCodebook <- VariablesinData[VariablesinData %notin% VariablesinVarTypes]
  variablesinCodebook_notinData <- VariablesinVarTypes[VariablesinVarTypes %notin% VariablesinData]

  # Initialize lists for warnings and recoded variables
  warninglist <- list()
  recodedvars <- list()

  # Check if there is a Recode column
  if (!("Recode" %in% colnames(VarTypes))) {
    VarTypes$Recode <- "NA"
  }

  # Check to see if there is a "Missing" Column, and if there is, add what the missing value should be
  if (!("Missing" %in% colnames(VarTypes))) {
    VarTypes$Missing <- missingVal
  }else{
    VarTypes$Missing[is.na(VarTypes$Missing)] <- missingVal
  }

  # Iterate over variables
  for (var in vars) {

    # Missing?

    mV <- VarTypes$Missing[VarTypes$Variable == var]
    x<- RevaluedData[[var]]
    x[x == mV] <- NA
    RevaluedData[[var]] <-x

    # Recode?
    rc <- VarTypes$Recode[VarTypes$Variable == var]
    if (!is.na(rc) && (rc == "yes" || rc == 1)) {
      recodedvars <- append(recodedvars, var)
      newcode <- VarTypes$Code[VarTypes$Variable == var]
      lookup <- data.frame(t(sapply(strsplit(as.character(newcode), splitchar)[[1]],
                                    function(x) strsplit(x, "=")[[1]])), row.names = NULL)
      lookup$X1 <- trimws(lookup$X1)
      lookup$X2 <- trimws(lookup$X2)
      freqs1 <- as.factor(RevaluedData[[var]]) %>% summary()
      nameList <- lookup$X2
      names(nameList) <- lookup$X1
      d <- sjlabelled::set_labels(RevaluedData[[var]], labels = nameList, force.labels = TRUE, force.values = TRUE)
      RevaluedData[[var]] <- sjlabelled::as_label(d)
      freqs2 <- as.factor(RevaluedData[[var]]) %>% summary()
      IsSame <- if_else(length(base::setdiff(freqs1, freqs2)) > 0, 0, 1)
      if (!IsSame) {
        warninglist <- append(warninglist, paste(var, ": final frequencies not consistent"))
      }
    } else {
      # Don't recode, but set to appropriate type

      type <- VarTypes$Type[VarTypes$Variable == var]
      if(!is.na(type)){
      if (type %in% c("Categorical", "categorical", "factor", "Factor")) {
        RevaluedData[[var]] <- to_factor(RevaluedData[[var]])
      }
      if (type %in% c("Double", "Numeric", "Numerical", "numeric", "numerical", "double")) {
        RevaluedData[[var]] <- as.numeric(RevaluedData[[var]])
      }
      if (type %in% c("Ordinal", "ordinal", "ordered factor") {
        RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
      }
      }else{
        RevaluedData[[var]] <- RevaluedData[[var]]
      }
    }

    # Fix the label
    newVarLabel <- VarTypes$Label[VarTypes$Variable == var]
    if (!is.null(newVarLabel)) {
      set_label(RevaluedData[var]) <- newVarLabel
    }
  }

  return(list(RevaluedData = RevaluedData, warninglist = warninglist, recodedvars = recodedvars))
}
