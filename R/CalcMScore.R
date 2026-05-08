#' Calculate M-scores using median and MAD
#'
#' Calculate robust standardized scores using the median as the center and
#' median absolute deviation as the scale. M-scores are useful for skewed,
#' heavy-tailed, or outlier-prone variables where standard z-scores may be
#' overly influenced by the mean and standard deviation.
#'
#' An M-score is calculated as:
#'
#' `(x - median(x)) / MAD(x)`
#'
#' By default, `stats::mad()` uses a scaling constant of `1.4826`, which makes
#' MAD approximately comparable to the standard deviation when data are normally
#' distributed.
#'
#' @param df Data frame with variables to robust-standardize.
#' @param variables Character vector of variable names. If `NULL`, uses
#'   `SciDataReportR::getNumVars(df, Ordinal = FALSE)`.
#' @param names_prefix Prefix to prepend to variable names. Default is `"M_"`.
#' @param RetainLabels Logical; if `TRUE` and `Hmisc` is available, copy labels.
#' @param RenameLabels Logical; if `TRUE`, apply the same prefix to labels.
#' @param center Logical; if `TRUE`, subtract the median.
#' @param scale Logical; if `TRUE`, divide by the MAD.
#' @param constant Numeric scaling constant passed to `stats::mad()`. Default is
#'   `1.4826`.
#'
#' @return An object of class `"MScoreObj"`, a list with:
#'   - `MScores`: data frame of M-score variables only
#'   - `DataWithM`: original `df` plus M-score variables
#'   - `Parameters`: data frame with `Variable`, `N`, `Median`, and `MAD`
#'   - `Center`: logical flag used
#'   - `Scale`: logical flag used
#'   - `Constant`: MAD scaling constant used
#' @export
CalcMScore <- function(df,
                       variables = NULL,
                       names_prefix = "M_",
                       RetainLabels = TRUE,
                       RenameLabels = TRUE,
                       center = TRUE,
                       scale = TRUE,
                       constant = 1.4826) {

  # Validate inputs

  if (!is.data.frame(df)) {
    stop("df must be a data frame.")
  }

  if (is.null(variables)) {
    if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
      stop("SciDataReportR is required to auto-detect numeric variables via getNumVars().")
    }

    variables <- SciDataReportR::getNumVars(df, Ordinal = FALSE)
  }

  if (!is.character(variables)) {
    stop("variables must be a character vector of variable names.")
  }

  missing_vars <- setdiff(variables, names(df))

  if (length(missing_vars) > 0) {
    stop(
      "The following variables are not in df: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  is_num <- vapply(df[variables], is.numeric, logical(1))

  if (!all(is_num)) {
    warning(
      "Dropping non-numeric variables from M-score calculation: ",
      paste(variables[!is_num], collapse = ", ")
    )

    variables <- variables[is_num]
  }

  if (length(variables) == 0) {
    stop("No numeric variables available to calculate M-scores.")
  }

  if (!is.character(names_prefix) || length(names_prefix) != 1 || is.na(names_prefix)) {
    stop("names_prefix must be a single character string.")
  }

  if (!is.logical(RetainLabels) || length(RetainLabels) != 1 || is.na(RetainLabels)) {
    stop("RetainLabels must be TRUE or FALSE.")
  }

  if (!is.logical(RenameLabels) || length(RenameLabels) != 1 || is.na(RenameLabels)) {
    stop("RenameLabels must be TRUE or FALSE.")
  }

  if (!is.logical(center) || length(center) != 1 || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }

  if (!is.logical(scale) || length(scale) != 1 || is.na(scale)) {
    stop("scale must be TRUE or FALSE.")
  }

  if (!is.numeric(constant) || length(constant) != 1 || is.na(constant) || constant <= 0) {
    stop("constant must be a single positive numeric value.")
  }

  # Compute parameters

  medians <- vapply(
    df[variables],
    function(x) stats::median(x, na.rm = TRUE),
    numeric(1)
  )

  mads <- vapply(
    df[variables],
    function(x) {
      med <- stats::median(x, na.rm = TRUE)

      stats::mad(
        x,
        center = med,
        constant = constant,
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  ns <- vapply(
    df[variables],
    function(x) sum(!is.na(x)),
    integer(1)
  )

  params <- data.frame(
    Variable = variables,
    N = ns,
    Median = medians,
    MAD = mads,
    stringsAsFactors = FALSE
  )

  # Compute M-scores

  transform_fun <- function(x, med, mad_value) {
    if (is.na(mad_value) || mad_value == 0) {
      mad_value <- NA_real_
    }

    out <- x

    if (center) {
      out <- out - med
    }

    if (scale) {
      if (is.na(mad_value)) {
        out[] <- NA_real_
      } else {
        out <- out / mad_value
      }
    }

    out
  }

  mscore_list <- mapply(
    FUN = transform_fun,
    df[variables],
    medians,
    mads,
    SIMPLIFY = FALSE
  )

  mscore_df <- as.data.frame(mscore_list, stringsAsFactors = FALSE)
  names(mscore_df) <- paste0(names_prefix, variables)

  # Handle labels

  if (RetainLabels && requireNamespace("Hmisc", quietly = TRUE)) {
    for (v in variables) {
      original_label <- Hmisc::label(df[[v]])

      if (!is.null(original_label) && nzchar(original_label)) {
        new_name <- paste0(names_prefix, v)

        if (RenameLabels) {
          Hmisc::label(mscore_df[[new_name]]) <- paste0(names_prefix, original_label)
        } else {
          Hmisc::label(mscore_df[[new_name]]) <- original_label
        }
      }
    }
  }

  # Return result

  combined <- cbind(df, mscore_df)

  out <- list(
    MScores = mscore_df,
    DataWithM = combined,
    Parameters = params,
    Center = center,
    Scale = scale,
    Constant = constant
  )

  class(out) <- c("MScoreObj", class(out))

  out
}
