# Internal display-label and ordinal helpers. These deliberately operate on
# copies so presentation choices never alter the caller's data frame.

ScidrDisplayLabels <- function(data, variables, Relabel = TRUE) {
  if (!is.logical(Relabel) || length(Relabel) != 1 || is.na(Relabel)) {
    stop("Relabel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(variables) || any(!variables %in% names(data))) {
    stop("variables must name columns in data.", call. = FALSE)
  }
  if (!Relabel) {
    labels <- variables
    categorical_prefix <- ".scidr_ordinal_categorical_"
    continuous_prefix <- ".scidr_ordinal_continuous_"
    is_categorical <- startsWith(variables, categorical_prefix)
    is_continuous <- startsWith(variables, continuous_prefix)
    labels[is_categorical] <- paste0(
      sub(categorical_prefix, "", variables[is_categorical], fixed = TRUE),
      " (categorical)"
    )
    labels[is_continuous] <- paste0(
      sub(continuous_prefix, "", variables[is_continuous], fixed = TRUE),
      " (continuous)"
    )
    return(stats::setNames(labels, variables))
  }

  labels <- sjlabelled::get_label(data[variables], def.value = variables)
  labels <- as.character(labels)
  labels[is.na(labels) | !nzchar(labels)] <- variables[is.na(labels) | !nzchar(labels)]
  stats::setNames(labels, variables)
}

ScidrApplyDisplayLabels <- function(data, variables, Relabel = TRUE) {
  labels <- ScidrDisplayLabels(data, variables, Relabel)
  for (var in variables) {
    data[[var]] <- sjlabelled::set_label(data[[var]], labels[[var]])
  }
  data
}

ScidrMatchOrdinalTreatment <- function(TreatOrdinalAs = "Categorical") {
  match.arg(
    TreatOrdinalAs,
    c("Categorical", "Continuous", "Both", "Exclude")
  )
}

ScidrIsOrdinal <- function(x) {
  is.ordered(x) || identical(attr(x, "scidr_type"), "ordinal")
}

ScidrOrdinalScores <- function(x) {
  if (!ScidrIsOrdinal(x)) return(NULL)

  levels_x <- levels(x)
  scores <- attr(x, "scidr_ordinal_scores")
  source <- attr(x, "scidr_ordinal_score_source")

  if (!is.null(scores) && !is.null(names(scores)) && all(levels_x %in% names(scores))) {
    return(list(
      scores = unname(as.numeric(scores[levels_x])),
      source = if (is.null(source)) "codebook" else source
    ))
  }

  list(scores = seq_along(levels_x), source = "rank")
}

ScidrOrdinalAsNumeric <- function(x) {
  score_info <- ScidrOrdinalScores(x)
  if (is.null(score_info)) return(x)

  out <- unname(score_info$scores[as.integer(x)])
  out[is.na(x)] <- NA_real_
  out <- sjlabelled::set_label(out, sjlabelled::get_label(x, def.value = NULL))
  attr(out, "scidr_ordinal_score_source") <- score_info$source
  out
}

ScidrAsPlainFactor <- function(x) {
  if (is.factor(x)) {
    return(factor(as.character(x), levels = levels(x)))
  }
  factor(x)
}

ScidrPrepareOrdinal <- function(data, variables, TreatOrdinalAs = "Categorical") {
  TreatOrdinalAs <- ScidrMatchOrdinalTreatment(TreatOrdinalAs)
  ordinal_vars <- variables[vapply(data[variables], ScidrIsOrdinal, logical(1))]

  if (TreatOrdinalAs == "Both") {
    stop(
      "TreatOrdinalAs = 'Both' is only supported by summary and comparison tables.",
      call. = FALSE
    )
  }
  if (TreatOrdinalAs == "Exclude") {
    variables <- setdiff(variables, ordinal_vars)
  } else if (TreatOrdinalAs == "Continuous") {
    for (var in ordinal_vars) data[[var]] <- ScidrOrdinalAsNumeric(data[[var]])
  }

  list(data = data, variables = variables, ordinal_variables = ordinal_vars,
       treatment = TreatOrdinalAs)
}

ScidrExpandOrdinalForTable <- function(data, variables,
                                       TreatOrdinalAs = "Categorical",
                                       Relabel = TRUE) {
  TreatOrdinalAs <- ScidrMatchOrdinalTreatment(TreatOrdinalAs)
  prep <- ScidrPrepareOrdinal(
    data, variables,
    if (TreatOrdinalAs == "Both") "Categorical" else TreatOrdinalAs
  )
  data <- prep$data
  variables <- prep$variables

  if (TreatOrdinalAs != "Both") {
    variable_map <- stats::setNames(as.list(variables), variables)
    return(list(data = data, variables = variables, variable_map = variable_map,
                treatment = TreatOrdinalAs))
  }

  ordinal_vars <- prep$ordinal_variables
  if (!length(ordinal_vars)) {
    variable_map <- stats::setNames(as.list(variables), variables)
    return(list(data = data, variables = variables, variable_map = variable_map,
                treatment = TreatOrdinalAs))
  }

  display_labels <- ScidrDisplayLabels(data, ordinal_vars, Relabel)
  expanded <- character(0)
  variable_map <- list()
  for (var in variables) {
    if (!var %in% ordinal_vars) {
      expanded <- c(expanded, var)
      variable_map[[var]] <- var
      next
    }
    cat_var <- paste0(".scidr_ordinal_categorical_", var)
    cont_var <- paste0(".scidr_ordinal_continuous_", var)
    data[[cat_var]] <- data[[var]]
    data[[cont_var]] <- ScidrOrdinalAsNumeric(data[[var]])
    data[[cat_var]] <- sjlabelled::set_label(
      data[[cat_var]], paste0(display_labels[[var]], " (categorical)")
    )
    data[[cont_var]] <- sjlabelled::set_label(
      data[[cont_var]], paste0(display_labels[[var]], " (continuous)")
    )
    expanded <- c(expanded, cat_var, cont_var)
    variable_map[[var]] <- c(cat_var, cont_var)
  }

  list(data = data, variables = expanded, variable_map = variable_map,
       treatment = TreatOrdinalAs)
}
