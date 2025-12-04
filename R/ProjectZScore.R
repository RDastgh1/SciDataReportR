#' Project standardized scores onto new data using external parameters
#'
#' @param df Data frame on which to project scores.
#' @param variables Character vector; if NULL, project onto all variables
#'   for which parameters exist and that are present in df.
#' @param parameters Source of parameters, interpreted by ParameterInputType:
#'   - "df_parameter": a data frame with cols Variable, N, Mean, SD
#'   - "ZScoreObj": output object from CalcZScore()
#'   - "ExternalDataframe": a raw reference data frame; parameters are
#'       estimated via CalcZScore() on that frame.
#' @param ParameterInputType One of "df_parameter", "ZScoreObj",
#'   "ExternalDataframe".
#' @param names_prefix Prefix for projected variable names.
#' @param RetainLabels Logical; if TRUE and Hmisc available, copy labels
#'   from df to new variables.
#' @param RenameLabels Logical; if TRUE, prefix labels the same way as names.
#' @param center Logical; used when ParameterInputType is "df_parameter"
#'   or "ExternalDataframe". Ignored for "ZScoreObj" (it uses stored flags).
#' @param scale Logical; same logic as `center`.
#'
#' @return List with same structure as CalcZScore():
#'   - ZScores: projected standardized variables only
#'   - DataWithZ: df with projected scores appended
#'   - Parameters: parameter data frame actually used for projection
#'   - Center, Scale: flags used
#' @export
Project_ZScore <- function(df,
                          variables = NULL,
                          parameters,
                          ParameterInputType = c("df_parameter",
                                                 "ZScoreObj",
                                                 "ExternalDataframe"),
                          names_prefix = "Z_",
                          RetainLabels = TRUE,
                          RenameLabels = TRUE,
                          center = TRUE,
                          scale  = TRUE) {

  ParameterInputType <- match.arg(ParameterInputType)

  # 1. Resolve parameter data frame and center/scale flags -----------------

  if (ParameterInputType == "df_parameter") {

    if (!is.data.frame(parameters)) {
      stop("For ParameterInputType = 'df_parameter', `parameters` must be a data frame.")
    }
    param_df <- parameters
    param_center <- center
    param_scale  <- scale

  } else if (ParameterInputType == "ZScoreObj") {

    if (!is.list(parameters) || is.null(parameters$Parameters)) {
      stop("For ParameterInputType = 'ZScoreObj', `parameters` must be a CalcZScore() output.")
    }
    param_df <- parameters$Parameters
    param_center <- isTRUE(parameters$Center)
    param_scale  <- isTRUE(parameters$Scale)

    # if user passes conflicting center/scale, warn
    if (center != param_center || scale != param_scale) {
      warning("center/scale arguments differ from those used in the ZScoreObj; ",
              "using the flags stored in the ZScoreObj (Center = ",
              param_center, ", Scale = ", param_scale, ").")
    }

  } else if (ParameterInputType == "ExternalDataframe") {

    ref_df <- parameters
    cs <- CalcZScore(
      df            = ref_df,
      variables     = variables,
      names_prefix  = names_prefix,
      RetainLabels  = FALSE,
      RenameLabels  = FALSE,
      center        = center,
      scale         = scale
    )
    param_df     <- cs$Parameters
    param_center <- cs$Center
    param_scale  <- cs$Scale

  } else {
    stop("Unknown ParameterInputType.")
  }

  required_cols <- c("Variable", "N", "Mean", "SD")
  missing_cols <- setdiff(required_cols, names(param_df))
  if (length(missing_cols) > 0) {
    stop("Parameter data frame must contain columns: ",
         paste(required_cols, collapse = ", "),
         ". Missing: ",
         paste(missing_cols, collapse = ", "))
  }

  # 2. Determine variables to project -------------------------------------

  if (is.null(variables)) {
    variables <- intersect(param_df$Variable, names(df))
    if (length(variables) == 0) {
      stop("No overlap between parameter variables and df columns.")
    }
  } else {
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

  is_num <- vapply(df[variables], is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric variables from projection: ",
            paste(variables[!is_num], collapse = ", "))
    variables <- variables[is_num]
  }
  if (length(variables) == 0) {
    stop("No numeric variables available to project scores.")
  }

  # Align params to variable order
  param_df_sub <- param_df[match(variables, param_df$Variable), ]
  means <- param_df_sub$Mean
  sds   <- param_df_sub$SD

  # 3. Compute projected scores -------------------------------------------

  transform_fun <- function(x, m, s, center_flag, scale_flag) {
    out <- x
    if (center_flag) {
      out <- out - m
    }
    if (scale_flag) {
      if (is.na(s) || s == 0) {
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
    MoreArgs = list(center_flag = param_center, scale_flag = param_scale),
    SIMPLIFY = FALSE
  )

  z_df <- as.data.frame(z_list, stringsAsFactors = FALSE)
  names(z_df) <- paste0(names_prefix, variables)

  # 4. Handle labels from df ----------------------------------------------

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
    Parameters = param_df_sub,
    Center     = param_center,
    Scale      = param_scale
  )
  class(out) <- c("ZScoreObj", class(out))
  out
}
