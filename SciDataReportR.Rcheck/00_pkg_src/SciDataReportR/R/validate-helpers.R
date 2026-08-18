# Internal input-validation helpers. These exist so that every function in the
# package reports the same user error the same way, naming the argument the user
# actually typed rather than an internal (often deprecated) variable name.

# Canonical check that a supplied variable set names columns in `data`.
#
# Returns the variables as a unique character vector, so callers can replace the
# old `unique(intersect(as.character(vars), names(Data)))` idiom outright. That
# idiom silently dropped misspelled names, which produced a smaller results
# matrix -- or, for covariates, statistics that were reported as adjusted while
# actually being unadjusted.
#
# `arg_name` is the user-facing argument name (e.g. "predictor_vars"), which is
# what appears in the error. `data_arg` is the user-facing data argument name.
ScidrValidateVariables <- function(data,
    variables,
    arg_name = "variables",
    data_arg = "data",
    allow_null = TRUE,
    unique_only = TRUE) {
  if (!is.data.frame(data)) {
    stop("`", data_arg, "` must be a data frame.", call. = FALSE)
  }
  if (is.null(variables) || length(variables) == 0) {
    if (!allow_null) {
      stop("`", arg_name, "` must name at least one column in `", data_arg, "`.",
           call. = FALSE)
    }
    return(character(0))
  }
  if (!is.character(variables)) {
    variables <- tryCatch(
      as.character(variables),
      error = function(e) {
        stop("`", arg_name, "` must be a character vector of column names.",
             call. = FALSE)
      }
    )
  }
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars)) {
    stop("Variables not found in `", data_arg, "`: ",
         paste(unique(missing_vars), collapse = ", "),
         " (supplied to `", arg_name, "`).", call. = FALSE)
  }
  if (unique_only) unique(variables) else variables
}

# Single-variable convenience wrapper, for arguments like `interVar` or
# `group_var` that name exactly one column.
ScidrValidateVariable <- function(data,
    variable,
    arg_name,
    data_arg = "data",
    allow_null = FALSE) {
  if (is.null(variable) || length(variable) == 0) {
    if (allow_null) return(NULL)
    stop("`", arg_name, "` must name a column in `", data_arg, "`.", call. = FALSE)
  }
  if (length(variable) != 1L) {
    stop("`", arg_name, "` must name exactly one column in `", data_arg,
         "`; received ", length(variable), ".", call. = FALSE)
  }
  ScidrValidateVariables(data, variable, arg_name, data_arg, allow_null = FALSE)
}
