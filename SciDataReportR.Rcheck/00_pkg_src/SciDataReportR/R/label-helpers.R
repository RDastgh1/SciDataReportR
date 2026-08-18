# Internal display-label helpers. These operate on copies so presentation
# choices never alter the caller's data frame.

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
