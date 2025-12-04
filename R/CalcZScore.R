#' Calculate Z-scores (or standardized scores) and return data + parameters
#'
#' @param df Data frame with variables to standardize.
#' @param variables Character vector of variable names. If NULL, uses
#'   SciDataReportR::getNumVars(df).
#' @param names_prefix Prefix to prepend to variable names (default "Z_").
#' @param RetainLabels Logical; if TRUE and Hmisc is available, copy labels.
#' @param RenameLabels Logical; if TRUE, apply the same prefix to labels.
#' @param center Logical; if TRUE, subtract the mean.
#' @param scale Logical; if TRUE, divide by the SD.
#'
#' @return An object of class "ZScoreObj", a list with:
#'   - ZScores: data frame of standardized variables only
#'   - DataWithZ: original df + standardized variables
#'   - Parameters: data frame with Variable, N, Mean, SD
#'   - Center: logical flag used
#'   - Scale: logical flag used
#' @export
CalcZScore <- function(df,
                       variables = NULL,
                       names_prefix = "Z_",
                       RetainLabels = TRUE,
                       RenameLabels = TRUE,
                       center = TRUE,
                       scale = TRUE) {
  # 1. Determine variables -------------------------------------------------
  if (is.null(variables)) {
    if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
      stop("SciDataReportR is required to auto-detect numeric variables via getNumVars().")
    }
    variables <- SciDataReportR::getNumVars(df, Ordinal = FALSE)
  }

  missing_vars <- setdiff(variables, names(df))
  if (length(missing_vars) > 0) {
    stop("The following variables are not in df: ",
         paste(missing_vars, collapse = ", "))
  }

  is_num <- vapply(df[variables], is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric variables from standardization: ",
            paste(variables[!is_num], collapse = ", "))
    variables <- variables[is_num]
  }

  if (length(variables) == 0) {
    stop("No numeric variables available to standardize.")
  }

  # 2. Compute parameters --------------------------------------------------
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

  # 3. Compute standardized scores ----------------------------------------
  transform_fun <- function(x, m, s) {
    # protect against 0/NA SD
    if (is.na(s) || s == 0) s <- NA_real_

    # start from raw x
    out <- x

    if (center) {
      out <- out - m
    }
    if (scale) {
      # if s is NA after protection, result is NA
      if (is.na(s)) {
        out[] <- NA_real_
      } else {
        out <- out / s
      }
    }
    out
  }

  z_list <- mapply(
    FUN = transform_fun,
    df[variables],
    means,
    sds,
    SIMPLIFY = FALSE
  )

  z_df <- as.data.frame(z_list, stringsAsFactors = FALSE)
  names(z_df) <- paste0(names_prefix, variables)

  # 4. Handle labels -------------------------------------------------------
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

  # 5. Combine and return --------------------------------------------------
  combined <- cbind(df, z_df)

  out <- list(
    ZScores    = z_df,
    DataWithZ  = combined,
    Parameters = params,
    Center     = center,
    Scale      = scale
  )
  class(out) <- c("ZScoreObj", class(out))
  out
}
