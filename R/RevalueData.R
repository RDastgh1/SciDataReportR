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
  VarTypes <- VarTypes %>% filter(!is.na(Variable))
  RevaluedData <- DatatoRevalue
  VariablesinData <- colnames(DatatoRevalue)
  VariablesinVarTypes <- VarTypes$Variable
  vars <- VariablesinData[VariablesinData %in% VariablesinVarTypes]
  warninglist <- list()
  recodedvars <- list()

  if (!("Recode" %in% colnames(VarTypes))) {
    VarTypes$Recode <- "NA"
  }
  if (!("Missing" %in% colnames(VarTypes))) {
    VarTypes$Missing <- missingVal
  } else {
    VarTypes$Missing[is.na(VarTypes$Missing)] <- missingVal
  }

  for (var in vars) {
    mV <- VarTypes$Missing[VarTypes$Variable == var]
    x <- RevaluedData[[var]]
    x[x == mV] <- NA
    RevaluedData[[var]] <- x
    rc <- VarTypes$Recode[VarTypes$Variable == var]

    if (!is.na(rc) && (rc == "yes" || rc == 1)) {
      recodedvars <- append(recodedvars, var)
      newcode <- VarTypes$Code[VarTypes$Variable == var]

      # Try to create lookup and catch errors
      lookup <- tryCatch({
        df <- data.frame(t(sapply(strsplit(as.character(newcode), splitchar)[[1]], function(x) {
          strsplit(x, "=")[[1]]
        })), row.names = NULL)

        df$X1 <- tryCatch(trimws(df$X1), error = function(e) {
          warninglist <<- append(warninglist, paste(var, ": recoding failed at trimws(X1) -", e$message))
          return(NULL)
        })

        df$X2 <- tryCatch(trimws(df$X2), error = function(e) {
          warninglist <<- append(warninglist, paste(var, ": recoding failed at trimws(X2) -", e$message))
          return(NULL)
        })

        df
      }, error = function(e) {
        warninglist <<- append(warninglist, paste(var, ": recoding failed during lookup creation -", e$message))
        return(NULL)
      })

      if (!is.null(lookup)) {
        freqs1 <- as.factor(RevaluedData[[var]]) %>% summary()
        nameList <- lookup$X2
        names(nameList) <- lookup$X1
        d <- sjlabelled::set_labels(RevaluedData[[var]],
                                    labels = nameList, force.labels = TRUE, force.values = TRUE)
        RevaluedData[[var]] <- sjlabelled::as_label(d)
        freqs2 <- as.factor(RevaluedData[[var]]) %>% summary()

        if (length(base::setdiff(freqs1, freqs2)) > 0) {
          warninglist <- append(warninglist, paste(var, ": final frequencies not consistent"))
        }
      }
    } else {
      type <- VarTypes$Type[VarTypes$Variable == var]
      if (!is.na(type)) {
        if (type %in% c("Categorical", "categorical", "factor", "Factor")) {
          RevaluedData[[var]] <- to_factor(RevaluedData[[var]])
        }
        if (type %in% c("Double", "Numeric", "Numerical", "numeric", "numerical", "double")) {
          RevaluedData[[var]] <- as.numeric(RevaluedData[[var]])
        }
        if (type %in% c("Ordinal", "ordinal", "ordered factor")) {
          RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        }
      }
    }

    newVarLabel <- VarTypes$Label[VarTypes$Variable == var]
    if (!is.null(newVarLabel)) {
      set_label(RevaluedData[var]) <- newVarLabel
    }
  }

  return(list(RevaluedData = RevaluedData, warninglist = warninglist, recodedvars = recodedvars))
}
