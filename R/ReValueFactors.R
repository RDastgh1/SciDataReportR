#' Revalue Factors
#'
#' This function revalues factor variables in a dataset according to the specifications provided in a codebook.
#'
#' @param DatatoRevalue The dataset to be revalued.
#' @param VarTypes A data frame containing information about the variables and how they should be revalued.
#'                 It should have columns: Variable (variable names), Recode (yes/no for recoding),
#'                 and Code (the revalue codes separated by "=" and ","). DO NOT USE COMMAS ANYWHERE ELSE IN THIS COLUMN.
#' @return The revalued dataset.
#' @export
ReValueFactors <- function(DatatoRevalue, VarTypes) {
  # Initialize the revalued data with the original data
  RevaluedData <- DatatoRevalue

  # Get variables in data that are also in VarTypes
  vars <- VarTypes$Variable[VarTypes$Variable %in% colnames(RevaluedData)]

  # Iterate over variables
  for (var in vars) {
    rc <- VarTypes$Recode[VarTypes$Variable == var]
    if (!is.na(rc) && rc == "yes") {
      # Split the code and create lookup table
      lookup <- data.frame(t(sapply(strsplit(as.character(VarTypes$Code[VarTypes$Variable == var]), ",")[[1]],
                                    function(x) strsplit(x, "=")[[1]])), row.names = NULL)
      lookup$X1 <- trimws(lookup$X1)
      lookup$X2 <- trimws(lookup$X2)

      # Match original values to revalue codes
      newvarind <- match(RevaluedData[[var]], lookup$X1, nomatch = -999)
      newvarmatch <- as.character(RevaluedData[[var]])
      newvarmatch[newvarind != -999] <- lookup$X2[newvarind[newvarind != -999]]

      # Revalue the factor variable
      RevaluedData[[var]] <- as.factor(newvarmatch)
    }
  }

  return(RevaluedData)
}
