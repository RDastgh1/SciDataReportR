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

# Prepare only inputs whose public role is categorical. Warning state lives in
# the calling stack, never in the caller's data or package-global options.
ScidrPrepareCategorical <- function(data, variables, arg_name) {
  vars_Categorical <- ScidrValidateVariables(data, variables, arg_name)
  if (!length(vars_Categorical)) return(data)

  caller <- parent.frame()
  warning_state <- NULL
  for (frame in rev(sys.frames()[-length(sys.frames())])) {
    if (exists(".ScidrCategoricalWarnings", envir = frame, inherits = FALSE)) {
      warning_state <- get(".ScidrCategoricalWarnings", envir = frame, inherits = FALSE)
      break
    }
  }
  if (is.null(warning_state)) warning_state <- new.env(parent = emptyenv())
  assign(".ScidrCategoricalWarnings", warning_state, envir = caller)

  for (variable in vars_Categorical) {
    x <- data[[variable]]
    function_name <- paste(deparse(sys.call(sys.parent())[[1L]]), collapse = "")
    context <- paste0(function_name, "(): `", arg_name, " = \"", variable, "\"`")
    if (!is.null(dim(x)) || !(is.factor(x) || is.character(x) ||
        is.logical(x) || (is.numeric(x) && !is.complex(x) &&
        !inherits(x, c("Date", "POSIXt", "difftime"))))) {
      stop(context, " requires a categorical vector; unsupported column type.", call. = FALSE)
    }
    observed <- unique(x[!is.na(x)])
    n_groups <- length(observed)
    if (!n_groups) {
      stop(context, " has no observed non-missing values for categorical grouping.", call. = FALSE)
    }
    examples <- paste(utils::head(as.character(observed), 5L), collapse = ", ")
    warning_text <- NULL
    if (is.numeric(x)) {
      valid_codes <- all(is.finite(observed)) && all(observed == trunc(observed))
      if (!valid_codes || n_groups > 10L) {
        stop(context, " is numeric with ", n_groups, " observed non-missing values (examples: ",
          examples, "). Categorical inputs are converted automatically only for finite whole-number codes with at most 10 observed values. ",
          "If these values intentionally represent categories, explicitly convert `", variable,
          "` to a factor before calling this function.", call. = FALSE)
      }
      warning_text <- paste0(context, " is numeric with ", n_groups,
        " observed non-missing values (examples: ", examples,
        "); converting locally to a factor. Few unique values do not establish categorical meaning. ",
        "Explicitly convert this column to a factor if categorical treatment is intended.")
      label <- attr(x, "label", exact = TRUE)
      data[[variable]] <- factor(x, levels = sort(observed))
      if (!is.null(label)) attr(data[[variable]], "label") <- label
    } else if (n_groups > 10L) {
      warning_text <- paste0(context, " has ", n_groups,
        " observed categories; more than 10 groups may make plots or comparisons difficult to interpret. Continuing with the declared categories.")
    }
    if (!is.null(warning_text)) {
      # The same column can occur in multiple categorical roles or wrappers.
      if (!isTRUE(warning_state[[variable]])) {
        warning_state[[variable]] <- TRUE
        warning(warning_text, call. = FALSE)
      }
    }
  }
  data
}
