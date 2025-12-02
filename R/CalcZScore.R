#' Calculate Z-scores and return data + parameters
#'
#' @param df Data frame with variables to standardize
#' @param variables Character vector of variable names. If NULL, uses
#'   SciDataReportR::getNumVars(df).
#' @param names_prefix Prefix to prepend to variable names (default "Z_").
#' @param RetainLabels Logical; if TRUE and Hmisc is available, copy labels.
#' @param RenameLabels Logical; if TRUE, apply the same prefix to labels.
#'
#' @return An object of class "ZScoreObj", a list with:
#'   - ZScores: data frame of z-scored variables only
#'   - DataWithZ: original df + z-score variables
#'   - Parameters: data frame with Variable, Mean, SD
#'
CalcZScore <- function(df,
                       variables = NULL,
                       names_prefix = "Z_",
                       RetainLabels = TRUE,
                       RenameLabels = TRUE) {
  # Determine variables -----------------------------------------------
  if (is.null(variables)) {
    if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
      stop("SciDataReportR is required to auto-detect numeric variables via getNumVars().")
    }
    variables <- SciDataReportR::getNumVars(df, Ordinal = FALSE)
  }

  # Ensure variables exist in df
  missing_vars <- setdiff(variables, names(df))
  if (length(missing_vars) > 0) {
    stop("The following variables are not in df: ",
         paste(missing_vars, collapse = ", "))
  }

  # Drop any non-numeric variables (in case user passed weird stuff)
  is_num <- vapply(df[variables], is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric variables from z-score calculation: ",
            paste(variables[!is_num], collapse = ", "))
    variables <- variables[is_num]
  }

  if (length(variables) == 0) {
    stop("No numeric variables available to z-score.")
  }

  #Compute parameters -------------------------------------------------
  means <- vapply(df[variables], function(x) mean(x, na.rm = TRUE), numeric(1))
  sds   <- vapply(df[variables], function(x) stats::sd(x, na.rm = TRUE), numeric(1))
  ns    <- vapply(df[variables], function(x) sum(!is.na(x)), integer(1))

  params <- data.frame(
    Variable = variables,
    N        = ns,
    Mean     = means,
    SD       = sds,
    stringsAsFactors = FALSE
  )
  # Compute Z-scores ---------------------------------------------------
  # Handle SD == 0 or NA by returning all-NA for that variable
  z_list <- mapply(
    FUN = function(x, m, s) {
      if (is.na(s) || s == 0) {
        return(rep(NA_real_, length(x)))
      }
      (x - m) / s
    },
    df[variables],
    means,
    sds,
    SIMPLIFY = FALSE
  )

  z_df <- as.data.frame(z_list, stringsAsFactors = FALSE)
  names(z_df) <- paste0(names_prefix, variables)

#Handle labels ------------------------------------------------------
  if (RetainLabels && requireNamespace("Hmisc", quietly = TRUE)) {
    for (v in variables) {
      original_label <- Hmisc::label(df[[v]])
      if (!is.null(original_label) && nzchar(original_label)) {
        new_name <- paste0(names_prefix, v)
        if (RenameLabels) {
          Hmisc::label(z_df[[new_name]]) <- paste0(names_prefix, original_label)
        } else {
          Hmisc::label(z_df[[new_name]]) <- original_label
        }
      }
    }
  }

#Combine and return -------------------------------------------------
  combined <- cbind(df, z_df)

  out <- list(
    ZScores    = z_df,
    CombinedScores  = combined,
    Parameters = params
  )
  class(out) <- c("ZScoreObj", class(out))
  out
}
