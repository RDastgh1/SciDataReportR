#' Calculate Z-scores (or standardized scores) and return data + parameters
#'
#' Standardizes each variable to a common scale and, critically, returns the
#' constants used to do it so the identical transformation can be replayed on
#' other data later.
#'
#' @section The equation:
#' For each variable, every value is expressed as its distance from that
#' variable's mean, measured in standard deviations:
#'
#' \deqn{z_i = \frac{x_i - \bar{x}}{s}}{z = (x - mean) / sd}
#'
#' where \eqn{\bar{x}} and \eqn{s} are the mean and standard deviation of the
#' variable, both computed with `na.rm = TRUE`. The result has mean 0 and
#' standard deviation 1, so `z = -2` means the same thing - two standard
#' deviations below average - no matter which variable it came from.
#'
#' `center = FALSE` drops the \eqn{-\bar{x}} term and `scale = FALSE` drops the
#' division by \eqn{s}, which is how the "Center Only" and "Scale Only"
#' options elsewhere in the package are expressed.
#'
#' @section Why the parameters are returned:
#' The \eqn{\bar{x}} and \eqn{s} used for each variable are stored in
#' `Parameters`, and this is the whole point of returning an object rather
#' than just a standardized data frame.
#'
#' A z-score is only interpretable relative to the sample it was computed
#' from. If a validation cohort, a follow-up visit, or a new site is
#' standardized against *its own* mean and SD, then every group is centered at
#' 0 by construction and any real difference between them is scaled away - the
#' clinical cohort looks exactly like the reference cohort. Worse, the
#' resulting columns are not comparable to the original ones even though they
#' carry the same names.
#'
#' Passing this object to [ProjectZScore()] applies the *frozen* training
#' constants to the new data instead, so new observations land on the original
#' scale and group differences survive. The same principle is why the
#' clustering pipelines store a `ZScoreObject` alongside the fitted model: a
#' projected participant must pass through exactly the scaling the model was
#' trained on.
#'
#' @section What the example demonstrates:
#' The first pair of density plots shows why standardizing helps at all: before
#' it, each biomarker sits at its own location, so no single axis is meaningful
#' for all of them and they cannot be compared or pooled; after it, every
#' variable is centered at 0 with a standard deviation of 1, so they share one
#' axis and a value of -2 means the same thing everywhere.
#'
#' The second part splits the cohort in half and scores the same participants
#' two ways. Freezing the training mean and SD and applying them through
#' [ProjectZScore()] keeps whatever offset the new sample really has.
#' Re-standardizing the new half against itself forces every variable back to
#' mean 0, so the new cohort can never differ from the training cohort no
#' matter what the measured values were. That difference is the information
#' loss the stored parameters exist to prevent.
#'
#' @seealso [ProjectZScore()] to apply stored parameters to new data, and
#'   [CreateNormativeTScoreModel()] when the reference values should also be
#'   adjusted for covariates such as age or education.
#'
#' @param data Data frame with variables to standardize.
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
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' vars_Biomarkers <- c("AXL", "Adiponectin", "Ferritin", "MMP7", "tau", "p_tau")
#'
#' z_obj <- CreateZScoreObject(df_Labelled, variables = vars_Biomarkers)
#'
#' # The stored centering and scaling constants
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(z_obj$Parameters, full_width = TRUE)
#' )))
#'
#' # Before standardizing: each biomarker has its own measurement scale.
#' PlotContinuousDistributions(
#'   df_Labelled, variables = vars_Biomarkers, ncol = 3
#' )
#'
#' # After standardizing: every variable shares the same SD scale.
#' PlotContinuousDistributions(
#'   z_obj$ZScores, variables = names(z_obj$ZScores), ncol = 3
#' )
#'
#' # Split the cohort and score the second half two different ways
#' set.seed(11)
#' rows_Train <- sample(nrow(df_Labelled), floor(nrow(df_Labelled) / 2))
#' df_Train <- df_Labelled[rows_Train, ]
#' df_New <- df_Labelled[-rows_Train, ]
#'
#' z_Train <- CreateZScoreObject(df_Train, variables = vars_Biomarkers)
#'
#' # Right: apply the frozen training parameters
#' z_Projected <- ProjectZScore(
#'   df_New,
#'   variables = vars_Biomarkers,
#'   parameters = z_Train,
#'   ParameterInputType = "ZScoreObj"
#' )
#'
#' # Wrong: re-standardize the new data against itself
#' z_Restandardized <- CreateZScoreObject(df_New, variables = vars_Biomarkers)
#'
#' df_Compare <- data.frame(
#'   Variable = vars_Biomarkers,
#'   TrainingMean = round(z_Train$Parameters$Mean, 3),
#'   TrainingSD = round(z_Train$Parameters$SD, 3),
#'   Projected = round(vapply(
#'     z_Projected$ZScores, function(x) mean(as.numeric(x), na.rm = TRUE),
#'     numeric(1)
#'   ), 3),
#'   Restandardized = round(vapply(
#'     z_Restandardized$ZScores, function(x) mean(as.numeric(x), na.rm = TRUE),
#'     numeric(1)
#'   ), 3),
#'   row.names = NULL
#' )
#'
#' # Every `Restandardized` mean is 0 by construction
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_Compare, full_width = TRUE)
#' )))
#' }
#'
#' @param df \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
CreateZScoreObject <- function(data,
    variables = NULL,
    names_prefix = "Z_",
    RetainLabels = TRUE,
    RenameLabels = TRUE,
    center = TRUE,
    scale = TRUE,
    df = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(df)) {
    lifecycle::deprecate_warn("19.15.0", "CreateZScoreObject(df)", "CreateZScoreObject(data)")
    data <- df
  }
  if (!missing(data)) df <- data

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

#' @description `CalcZScore()` has been superseded by `CreateZScoreObject()`.
#'   It remains available as a backwards-compatible alias and returns the same
#'   reusable Z-score object.
#' @rdname CreateZScoreObject
#' @param ... Arguments passed to [CreateZScoreObject()].
#' @export
CalcZScore <- function(...) {
  CreateZScoreObject(...)
}
