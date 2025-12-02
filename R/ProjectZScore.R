
#' Project Z-scores using external parameters or reference data
#'
#' @param df Data frame on which to project Z-scores.
#' @param variables Character vector; if NULL, use all variables for which
#'   parameters exist (intersect with df).
#' @param parameters Source of parameters, interpreted according to
#'   ParameterInputType:
#'   - "df_parameter": a data frame with cols Variable, Mean, SD
#'   - "ZScoreObj": output object from CalcZScore()
#'   - "ExternalDataframe": a raw reference data frame; parameters will be
#'       estimated via CalcZScore() on that frame.
#' @param ParameterInputType One of "df_parameter", "ZScoreObj",
#'   "ExternalDataframe".
#' @param names_prefix Prefix for projected Z-score variable names.
#' @param RetainLabels Logical; if TRUE and Hmisc available, copy labels
#'   from df to new variables.
#' @param RenameLabels Logical; if TRUE, prefix labels the same way as names.
#'
#' @return List with same structure as CalcZScore():
#'   - ZScores: projected z-score variables only
#'   - DataWithZ: df with projected z-scores appended
#'   - Parameters: parameter data frame actually used for projection
#'
ProjectZScore <- function(df,
                          variables = NULL,
                          parameters,
                          ParameterInputType = c("df_parameter",
                                                 "ZScoreObj",
                                                 "ExternalDataframe"),
                          names_prefix = "Z_",
                          RetainLabels = TRUE,
                          RenameLabels = TRUE) {
  ParameterInputType <- match.arg(ParameterInputType)

  # 1. Resolve parameter data frame --------------------------------------
  param_df <- switch(
    ParameterInputType,
    "df_parameter" = {
      if (!is.data.frame(parameters)) {
        stop("For ParameterInputType = 'df_parameter', `parameters` must be a data frame.")
      }
      parameters
    },
    "ZScoreObj" = {
      if (!is.list(parameters) || is.null(parameters$Parameters)) {
        stop("For ParameterInputType = 'ZScoreObj', `parameters` must be a CalcZScore() output.")
      }
      parameters$Parameters
    },
    "ExternalDataframe" = {
      # Here `parameters` is actually an external reference data frame
      ref_df <- parameters
      cs <- CalcZScore(
        df            = ref_df,
        variables     = variables,
        names_prefix  = names_prefix,
        RetainLabels  = FALSE,   # labels of reference df not needed here
        RenameLabels  = FALSE
      )
      cs$Parameters
    }
  )

  # Basic sanity checks for param_df
  required_cols <- c("Variable", "Mean", "SD")
  missing_cols <- setdiff(required_cols, names(param_df))
  if (length(missing_cols) > 0) {
    stop("Parameter data frame must contain columns: ",
         paste(required_cols, collapse = ", "),
         ". Missing: ",
         paste(missing_cols, collapse = ", "))
  }

  # 2. Determine variables to project ------------------------------------
  if (is.null(variables)) {
    # "do all the variables that it has data on"
    variables <- intersect(param_df$Variable, names(df))
    if (length(variables) == 0) {
      stop("No overlap between parameter variables and df columns.")
    }
  } else {
    # Ensure requested vars exist in both df and parameter table
    missing_in_df <- setdiff(variables, names(df))
    if (length(missing_in_df) > 0) {
      stop("The following variables are not in df: ",
           paste(missing_in_df, collapse = ", "))
    }
    missing_in_params <- setdiff(variables, param_df$Variable)
    if (length(missing_in_params) > 0) {
      stop("The following variables have no parameters: ",
           paste(missing_in_params, collapse = ", "))
    }
  }

  # Drop non-numeric variables from df if user passed them
  is_num <- vapply(df[variables], is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric variables from z-score projection: ",
            paste(variables[!is_num], collapse = ", "))
    variables <- variables[is_num]
  }
  if (length(variables) == 0) {
    stop("No numeric variables available to project z-scores.")
  }

  # Align parameters to variable order
  param_df_sub <- param_df[match(variables, param_df$Variable), ]
  means <- param_df_sub$Mean
  sds   <- param_df_sub$SD

  # 3. Compute projected Z-scores ----------------------------------------
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

  # 4. Handle labels (from *df*, not from reference) ---------------------
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

  combined <- cbind(df, z_df)

  out <- list(
    ZScores    = z_df,
    DataWithZ  = combined,
    Parameters = param_df_sub
  )
  class(out) <- c("ZScoreObj", class(out))
  out
}
