#' Prepare ordinal variables for analysis
#'
#' Applies a consistent ordinal-treatment policy to selected variables. Ordinal
#' score mappings recorded by [RevalueData()] are used when available;
#' otherwise ordered-factor ranks are used.
#'
#' @section The four policies:
#' `"Continuous"` replaces each level with its rank. `"Categorical"` and
#' `"Exclude"` both leave the values alone; they differ in whether the variable
#' is offered to the analysis at all. `"Both"` produces one column of each.
#'
#' With `ReturnMetadata = TRUE` the derived column names come back alongside
#' the data, so downstream functions know which column carries which treatment
#' without having to guess from the naming convention.
#'
#' @param data The data frame containing the variables.
#' @param variables Character vector of variables to consider. If `NULL`, all
#'   columns are considered.
#' @param TreatOrdinalAs How ordinal variables are handled: `"Continuous"`,
#'   `"Categorical"`, `"Both"`, or `"Exclude"`.
#' @param Relabel Logical; when `TreatOrdinalAs = "Both"`, apply descriptive
#'   labels to the derived categorical and continuous variables.
#' @param ReturnMetadata Logical; if `FALSE` (default), return only the
#'   transformed data frame. If `TRUE`, return a list containing the data,
#'   selected variables, ordinal variables, variable map, and treatment.
#'
#' @return A transformed data frame, or a metadata list when
#'   `ReturnMetadata = TRUE`.
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id = 1:6,
#'   Severity = factor(
#'     c("None", "Mild", "Severe", "Mild", "Moderate", "None"),
#'     levels = c("None", "Mild", "Moderate", "Severe"), ordered = TRUE
#'   ),
#'   Education = factor(
#'     c("HighSchool", "College", "Graduate", "College", "Graduate", "HighSchool"),
#'     levels = c("HighSchool", "College", "Graduate"), ordered = TRUE
#'   )
#' )
#'
#' # The four policies, applied to the same ordinal variable
#' df_Both <- ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Both")
#'
#' df_Policies <- data.frame(
#'   id = df$id,
#'   Original = as.character(df$Severity),
#'   Continuous = ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Continuous")$Severity,
#'   Categorical = as.character(
#'     ConvertOrdinalToNumeric(df, TreatOrdinalAs = "Categorical")$Severity
#'   ),
#'   Both_Categorical = as.character(df_Both$.scidr_ordinal_categorical_Severity),
#'   Both_Continuous = df_Both$.scidr_ordinal_continuous_Severity
#' )
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_Policies, full_width = TRUE)
#' )))
#'
#' # The derived column names, returned alongside the data
#' metadata <- ConvertOrdinalToNumeric(
#'   df, TreatOrdinalAs = "Both", ReturnMetadata = TRUE
#' )
#'
#' df_Map <- data.frame(
#'   Variable = rep(names(metadata$variable_map),
#'                  lengths(metadata$variable_map)),
#'   DerivedColumn = unlist(metadata$variable_map, use.names = FALSE),
#'   IsOrdinal = rep(names(metadata$variable_map),
#'                   lengths(metadata$variable_map)) %in%
#'     metadata$ordinal_variables
#' )
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_Map, full_width = TRUE)
#' )))
#' }
#'
#' @importFrom sjlabelled get_label set_label
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @export
ConvertOrdinalToNumeric <- function(data,
    variables = NULL,
    TreatOrdinalAs = c("Continuous", "Categorical", "Both", "Exclude"),
    Relabel = TRUE,
    ReturnMetadata = FALSE,
    Data = lifecycle::deprecated(),
    Variables = lifecycle::deprecated()) {
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "ConvertOrdinalToNumeric(Data)", "ConvertOrdinalToNumeric(data)")
    data <- Data
  }
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "ConvertOrdinalToNumeric(Variables)", "ConvertOrdinalToNumeric(variables)")
    variables <- Variables
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (is.null(variables)) variables <- names(data)
  if (!is.character(variables)) {
    stop("`variables` must be a character vector or NULL.", call. = FALSE)
  }
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars)) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if (!is.logical(Relabel) || length(Relabel) != 1 || is.na(Relabel)) {
    stop("`Relabel` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(ReturnMetadata) || length(ReturnMetadata) != 1 || is.na(ReturnMetadata)) {
    stop("`ReturnMetadata` must be TRUE or FALSE.", call. = FALSE)
  }
  TreatOrdinalAs <- match.arg(TreatOrdinalAs)

  IsOrdinal <- function(x) {
    is.ordered(x) || identical(attr(x, "scidr_type"), "ordinal")
  }
  OrdinalAsNumeric <- function(x) {
    if (!IsOrdinal(x)) return(x)

    levels_x <- levels(x)
    scores <- attr(x, "scidr_ordinal_scores")
    score_source <- attr(x, "scidr_ordinal_score_source")
    if (!is.null(scores) && !is.null(names(scores)) && all(levels_x %in% names(scores))) {
      scores <- unname(as.numeric(scores[levels_x]))
      score_source <- if (is.null(score_source)) "codebook" else score_source
    } else {
      scores <- seq_along(levels_x)
      score_source <- "rank"
    }

    out <- unname(scores[as.integer(x)])
    out[is.na(x)] <- NA_real_
    out <- sjlabelled::set_label(out, sjlabelled::get_label(x, def.value = NULL))
    attr(out, "scidr_ordinal_score_source") <- score_source
    out
  }

  ordinal_variables <- variables[vapply(data[variables], IsOrdinal, logical(1))]
  output_variables <- variables
  variable_map <- stats::setNames(as.list(variables), variables)

  if (TreatOrdinalAs == "Continuous") {
    for (var in ordinal_variables) data[[var]] <- OrdinalAsNumeric(data[[var]])
  } else if (TreatOrdinalAs == "Exclude") {
    output_variables <- setdiff(variables, ordinal_variables)
    variable_map[ordinal_variables] <- rep(list(character(0)), length(ordinal_variables))
  } else if (TreatOrdinalAs == "Both" && length(ordinal_variables)) {
    display_labels <- ScidrDisplayLabels(data, ordinal_variables, Relabel)
    output_variables <- character(0)
    variable_map <- list()
    for (var in variables) {
      if (!var %in% ordinal_variables) {
        output_variables <- c(output_variables, var)
        variable_map[[var]] <- var
        next
      }
      categorical_var <- paste0(".scidr_ordinal_categorical_", var)
      continuous_var <- paste0(".scidr_ordinal_continuous_", var)
      data[[categorical_var]] <- sjlabelled::set_label(
        data[[var]], paste0(display_labels[[var]], " (categorical)")
      )
      data[[continuous_var]] <- sjlabelled::set_label(
        OrdinalAsNumeric(data[[var]]), paste0(display_labels[[var]], " (continuous)")
      )
      output_variables <- c(output_variables, categorical_var, continuous_var)
      variable_map[[var]] <- c(categorical_var, continuous_var)
    }
  }

  if (!ReturnMetadata) return(data)
  list(
    data = data,
    variables = output_variables,
    ordinal_variables = ordinal_variables,
    variable_map = variable_map,
    treatment = TreatOrdinalAs
  )
}
