# Codebook metadata is resolved independently of immutable detection evidence.
ScidrAnomalyCodebook <- function(codebook) {
  if (is.null(codebook)) return(NULL)
  if (!is.data.frame(codebook) || anyDuplicated(names(codebook)) || !"Variable" %in% names(codebook))
    stop("codebook must be a dataframe with a Variable column and unique column names.")
  for (field in intersect(c("Variable", "Label", "Category"), names(codebook))) {
    if (!is.character(codebook[[field]]) && !is.factor(codebook[[field]]) && !all(is.na(codebook[[field]]))) stop("codebook$", field, " must contain text.")
  }
  codebook <- as.data.frame(codebook)
  codebook$Variable <- as.character(codebook$Variable)
  if (anyNA(codebook$Variable) || any(!nzchar(trimws(codebook$Variable)))) stop("codebook$Variable must contain nonmissing, nonblank names.")
  # Only these fields are read, so duplicates are judged on them alone. A data
  # dictionary usually carries far more -- value ranges, missing counts, notes --
  # and two rows for one variable that disagree only there say nothing
  # contradictory about how the variable is labelled, categorised or typed.
  fields <- unique(c("Label", "Category", intersect(c("Type", "TimeUnit", "TimeFormat"), names(codebook))))
  for (field in fields) {
    if (!field %in% names(codebook)) codebook[[field]] <- NA_character_
    if (!is.list(codebook[[field]])) {
      codebook[[field]] <- trimws(as.character(codebook[[field]]))
      codebook[[field]][is.na(codebook[[field]]) | !nzchar(codebook[[field]])] <- NA_character_
    }
  }
  codebook <- unique(codebook[c("Variable", fields)])
  duplicates <- unique(codebook$Variable[duplicated(codebook$Variable)])
  if (length(duplicates)) {
    Stated <- function(variable, field) {
      values <- unique(codebook[[field]][codebook$Variable == variable])
      values[!is.na(values)]
    }
    disagreeing <- lapply(duplicates, function(variable)
      fields[vapply(fields, function(field) length(Stated(variable, field)) > 1L, logical(1))])
    conflicted <- duplicates[lengths(disagreeing) > 0L]
    if (length(conflicted)) {
      Describe <- function(i) paste0(conflicted[i], " (", paste(vapply(
        disagreeing[lengths(disagreeing) > 0L][[i]], function(field)
          paste0(field, ": ", paste(Stated(conflicted[i], field), collapse = " vs ")),
        character(1)), collapse = "; "), ")")
      shown <- seq_len(min(5L, length(conflicted)))
      stop("Conflicting duplicate codebook entries: ",
        paste(vapply(shown, Describe, character(1)), collapse = "; "),
        if (length(conflicted) > length(shown)) paste0("; and ", length(conflicted) - length(shown), " more"),
        ". Keep one row per variable, or drop the disagreeing field.", call. = FALSE)
    }
    # No row disagrees with another about anything that is actually stated, so
    # the duplicates differ only in what one of them leaves empty. A blank says
    # nothing rather than contradicting something, so the rows are merged into
    # the one that carries the information.
    merged <- codebook[!duplicated(codebook$Variable), , drop = FALSE]
    for (variable in duplicates) for (field in fields) {
      values <- Stated(variable, field)
      if (length(values)) merged[[field]][merged$Variable == variable] <- values
    }
    codebook <- merged
  }
  codebook[unique(c("Variable", "Label", "Category", intersect(c("Type", "TimeUnit", "TimeFormat"), names(codebook))))]
}

ScidrAnomalyVariableMetadata <- function(source, codebook = NULL, saved = source$Metadata$VariableMetadata) {
  variables <- setdiff(names(source$ReviewData), ".Row")
  labels <- unname(source$Metadata$Labels[variables])
  if (!length(labels)) labels <- rep(NA_character_, length(variables))
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- variables[missing]
  attached <- vapply(source$ReviewData[variables], function(x) {
    label <- attr(x, "label", exact = TRUE)
    length(label) > 0 && !is.na(label[1]) && nzchar(label[1])
  }, logical(1))
  label_source <- ifelse(labels == variables & !attached, "name fallback", "attached label")
  categories <- rep("Uncategorized", length(variables))
  category_source <- rep("missing", length(variables))
  if (!is.null(saved)) {
    index <- match(variables, saved$Variable)
    take <- !is.na(index) & saved$LabelSource[index] == "codebook"
    labels[take] <- saved$Label[index[take]]
    label_source[take] <- "codebook"
    take <- !is.na(index) & !is.na(saved$Category[index]) & nzchar(saved$Category[index])
    categories[take] <- saved$Category[index[take]]
    category_source[take] <- if ("CategorySource" %in% names(saved)) saved$CategorySource[index[take]] else ifelse(categories[take] == "Uncategorized", "missing", "saved")
  }
  if (!is.null(codebook)) {
    index <- match(variables, codebook$Variable)
    take <- !is.na(index) & !is.na(codebook$Label[index])
    labels[take] <- codebook$Label[index[take]]
    label_source[take] <- "codebook"
    take <- !is.na(index) & !is.na(codebook$Category[index])
    categories[take] <- codebook$Category[index[take]]
    category_source[take] <- "codebook"
  }
  relabel <- source$Metadata$Relabel
  if (is.null(relabel)) relabel <- !identical(source$Metadata$Call$Relabel, FALSE)
  if (!relabel) labels <- variables
  df_Metadata <- tibble::tibble(Variable = variables, Label = labels, Category = categories, LabelSource = label_source, CategorySource = category_source)
  diagnostics <- list()
  Add <- function(issue, values) {
    if (length(values)) diagnostics[[length(diagnostics) + 1L]] <<- tibble::tibble(Issue = issue, Detail = paste(values, collapse = ", "))
  }
  if (!is.null(codebook)) {
    Add("codebook_unmatched_variables", setdiff(variables, codebook$Variable))
    Add("codebook_unused_entries", setdiff(codebook$Variable, variables))
    Add("codebook_missing_labels", variables[label_source == "name fallback"])
    Add("codebook_missing_categories", variables[category_source == "missing"])
  }
  list(Variables = df_Metadata, Diagnostics = dplyr::bind_rows(c(list(tibble::tibble(Issue = character(), Detail = character())), diagnostics)))
}

ScidrAnomalyAttachCodebook <- function(source, codebook) {
  resolved <- ScidrAnomalyVariableMetadata(source, codebook)
  source$Metadata$VariableMetadata <- resolved$Variables
  source$Metadata$CodebookDiagnostics <- resolved$Diagnostics
  if (!is.null(codebook)) {
    supplied <- resolved$Variables$Variable[resolved$Variables$LabelSource == "codebook"]
    source$Diagnostics <- source$Diagnostics[!(source$Diagnostics$Issue == "missing_label" & source$Diagnostics$Detail %in% supplied), ]
    df_Diagnostics <- resolved$Diagnostics
    df_Diagnostics$Method <- rep(NA_character_, nrow(df_Diagnostics))
    df_Diagnostics$Stratum <- rep(NA_integer_, nrow(df_Diagnostics))
    df_Diagnostics$Rows <- rep(list(integer()), nrow(df_Diagnostics))
    source$Diagnostics <- dplyr::bind_rows(source$Diagnostics, df_Diagnostics)
  }
  source
}

ScidrReviewCodebook <- function(object, codebook = NULL) {
  catalogs <- diagnostics <- list()
  for (analysis in names(object$Sources)) {
    source <- object$Sources[[analysis]]
    saved <- object$Metadata$VariableMetadata[[analysis]]
    if (is.null(saved)) saved <- source$Metadata$VariableMetadata
    resolved <- ScidrAnomalyVariableMetadata(source, codebook, saved)
    if (!is.null(codebook) && "Type" %in% names(codebook)) {
      requested <- tolower(trimws(as.character(codebook$Type)))
      selected <- codebook$Variable %in% source$Metadata$Variables & requested %in% c("clock", "duration")
      current <- source$Metadata$TimeConversions
      actual <- if (is.null(current)) rep(NA_character_, nrow(codebook)) else current$TimeType[match(codebook$Variable, current$Variable)]
      different <- which(selected & (is.na(actual) | actual != requested))
      if (length(different)) resolved$Diagnostics <- dplyr::bind_rows(resolved$Diagnostics, tibble::tibble(
        Issue = "temporal_declarations_not_applied", Detail = paste0(paste(codebook$Variable[different], collapse = ", "),
          ": codebook temporal declarations do not change saved models in a review view. Rerun DetectDataAnomalies() to apply them.")))
    }
    catalogs[[analysis]] <- resolved$Variables
    if (nrow(resolved$Diagnostics)) {
      resolved$Diagnostics$Analysis <- analysis
      diagnostics[[analysis]] <- resolved$Diagnostics
    }
  }
  object$Metadata$VariableMetadata <- catalogs
  if (length(diagnostics)) object$Diagnostics <- dplyr::distinct(dplyr::bind_rows(object$Diagnostics, dplyr::bind_rows(diagnostics)))
  object
}

ScidrReviewCatalog <- function(object, analysis) {
  catalog <- object$Metadata$VariableMetadata[[analysis]]
  source <- object$Sources[[analysis]]
  if (is.null(catalog) || !all(setdiff(names(source$ReviewData), ".Row") %in% catalog$Variable)) catalog <- ScidrAnomalyVariableMetadata(source, saved = catalog)$Variables
  labels <- unname(source$Metadata$Labels[catalog$Variable])
  take <- catalog$LabelSource != "codebook" & !is.na(labels) & nzchar(labels)
  catalog$Label[take] <- labels[take]
  if (identical(source$Metadata$Relabel, FALSE) || (is.null(source$Metadata$Relabel) && identical(source$Metadata$Call$Relabel, FALSE))) catalog$Label <- catalog$Variable
  catalog
}

ScidrReviewDisplaySource <- function(object, analysis) {
  source <- object$Sources[[analysis]]
  catalog <- ScidrReviewCatalog(object, analysis)
  source$Metadata$Labels[catalog$Variable] <- catalog$Label
  source$Metadata$VariableMetadata <- catalog
  source
}

ScidrAnomalyDisplayFields <- function(df_Data, catalog) {
  if (is.null(catalog)) return(df_Data)
  for (pair in list(c("Variable", "Label"), c("OtherVariable", "OtherLabel"))) {
    if (!pair[1] %in% names(df_Data)) next
    index <- match(df_Data[[pair[1]]], catalog$Variable)
    if (!pair[2] %in% names(df_Data)) df_Data[[pair[2]]] <- rep(NA_character_, nrow(df_Data))
    take <- !is.na(index)
    df_Data[[pair[2]]][take] <- catalog$Label[index[take]]
  }
  df_Data
}

ScidrReviewInspectionColumns <- function(source, row) {
  df_Values <- source$FeatureSensitivity
  effect <- "ScoreDecrease"
  if (is.null(df_Values) || !nrow(df_Values)) {
    df_Values <- source$FeatureContributions
    effect <- "Contribution"
  }
  if (is.null(df_Values) || !nrow(df_Values)) return(character())
  df_Values <- df_Values[df_Values$.Row == row & is.finite(df_Values[[effect]]) & df_Values[[effect]] > 0, ]
  head(df_Values$Variable[order(-df_Values[[effect]], df_Values$Variable)], 5L)
}

ScidrReviewCategoryFields <- function(object, df_Data) {
  categories <- basis <- variables <- conditions <- keys <- character(nrow(df_Data))
  for (analysis in unique(df_Data$Analysis)) {
    catalog <- ScidrReviewCatalog(object, analysis)
    positions <- which(df_Data$Analysis == analysis)
    df_Data[positions, ] <- ScidrAnomalyDisplayFields(df_Data[positions, , drop = FALSE], catalog)
    for (i in positions) {
      fields <- intersect(unique(c(df_Data$Variable[i], df_Data$OtherVariable[i])), catalog$Variable)
      suggested <- !length(fields)
      if (suggested && (if ("Method" %in% names(df_Data)) df_Data$Method[i] == "isolation" else
          any(object$Evidence$IssueID == df_Data$IssueID[i] & object$Evidence$Method == "isolation"))) {
        fields <- ScidrReviewInspectionColumns(object$Sources[[analysis]], df_Data$SourceRow[i])
        fields <- intersect(fields, catalog$Variable)
      }
      members <- if (length(fields)) unique(catalog$Category[match(fields, catalog$Variable)]) else "Whole record\u2014category unspecified"
      categories[i] <- paste(members, collapse = "; ")
      keys[i] <- as.character(jsonlite::toJSON(unname(members), auto_unbox = FALSE))
      basis[i] <- if (!length(fields)) "Whole record category unspecified" else if (suggested) "Suggested inspection categories" else "Affected variable categories"
      variables[i] <- paste(vapply(fields, function(variable) {
        index <- match(variable, catalog$Variable)
        paste0(catalog$Label[index], " (", variable, "): ", catalog$Category[index])
      }, character(1)), collapse = "; ")
      if ("Conditions" %in% names(df_Data)) {
        condition_vars <- unique(unlist(lapply(df_Data$Conditions[[i]], function(x) x$column), use.names = FALSE))
        condition_vars <- intersect(condition_vars, catalog$Variable)
        conditions[i] <- paste(vapply(condition_vars, function(variable) {
          index <- match(variable, catalog$Variable)
          paste0(catalog$Label[index], " (", variable, "): ", catalog$Category[index])
        }, character(1)), collapse = "; ")
      }
    }
  }
  df_Data$CodebookCategory <- categories
  df_Data$CodebookCategoryKeys <- keys
  df_Data$CategoryBasis <- basis
  df_Data$CategoryVariables <- variables
  df_Data$ConditionCategories <- conditions
  df_Data
}

ScidrReviewPresentationFields <- function() c("Label", "OtherLabel", "WhyFlagged", "ComparisonDetails",
  "CodebookCategory", "CodebookCategoryKeys", "CategoryBasis", "CategoryVariables", "ConditionCategories")

ScidrAnomalyMethods <- function() data.frame(
  Key = c("missingness", "conditional", "isolation", "dates", "local", "distance", "cellwise", "change", "trajectory", "transition"),
  Method = c("Unexpected missingness", "Conditional pattern", "Isolation forest", "Date screening", "Local density", "Robust multivariate distance", "Cellwise consistency", "Change rate", "Trajectory residual", "Categorical transition"),
  Meaning = c("Blank fields where comparable records almost always carry a value. Comparable means a rule learned from another column's values, or a group of variables that go missing together, in which this record carries most members. A blank explained by neither is not reported.",
    "A value or response unusual within a subgroup learned from other variables, or within the whole column.",
    "An unusual overall combination of values. Column screening provides inspection leads.",
    "Unusual dates or signed intervals between date columns.",
    "Records less densely surrounded by similar records than their neighbors.",
    "Numeric combinations far from the robust center, accounting for correlations.",
    "Individual numeric values inconsistent with estimates based on related columns.",
    "Unusual changes per unit time between observed visits.",
    "Observations far from an in-sample fitted longitudinal trajectory.",
    "Uncommon consecutive responses relative to other participants' transitions."), stringsAsFactors = FALSE)

# Review severity. Scores are not comparable across methods -- a missingness
# observed-fraction and a conditional boundary distance measure different things,
# and several methods emit native flags with no record score at all -- so the
# order is driven by what kind of finding it is, and the within-method percentile
# only separates findings of the same type.
ScidrAnomalySeverity <- function() data.frame(
  Key = c("inconsistent", "date", "date_interval", "conditional", "distribution",
    "change", "transition", "missingness", "coded", "combination"),
  Severity = c(1L, 2L, 2L, 3L, 4L, 5L, 5L, 6L, 7L, 8L),
  Rationale = c(
    "Value contradicts related fields in the same record",
    "Date outside its learned boundaries",
    "Interval between two dates outside its learned range",
    "Value outside the boundary learned for its subgroup",
    "Value extreme within its column",
    "Change between visits outside its learned range",
    "Uncommon consecutive responses",
    "Blank where comparable records carry a value",
    "Rare value in a small set of numeric codes",
    "Exploratory whole-record combination with no attributed field"),
  stringsAsFactors = FALSE)

# One canonical taxonomy. Both the per-evidence and per-issue type come from
# ScidrAnomalyType(), so the explorer, the workbook and the queue cannot disagree.
ScidrAnomalyTypes <- function() data.frame(
  Key = c("distribution", "conditional", "combination", "missingness", "inconsistent",
    "date", "date_interval", "change", "transition", "coded"),
  Type = c("Distribution outlier", "Conditional outlier", "Uncommon combination",
    "Unexpected missingness", "Inconsistent with related fields", "Implausible date",
    "Implausible date interval", "Unexpected change over time", "Uncommon transition",
    "Rare coded value"),
  Definition = c("A value that is unusual within its own column, across this analysis group.",
    "A value that is unusual for records sharing the comparison group learned from other columns.",
    "No single field looks wrong, but the overall combination of values is unusual.",
    "A field is blank where comparable records carry a value: records matching a rule learned from another column, or records carrying the same group of variables that go missing together.",
    "A numeric value that disagrees with estimates built from related numeric columns.",
    "A date outside the boundaries learned from the other values in its column.",
    "A signed interval between two date columns outside the learned boundaries.",
    "A change between visits, or a departure from a fitted trajectory, larger than usual.",
    "A consecutive pair of responses that few other participants record.",
    "A rare value in a column holding only a small set of numeric codes."),
  stringsAsFactors = FALSE)

# Rare numeric codes outrank the method-derived type. The longitudinal methods
# are typed from the method before the no-field test: a change finding is
# attributed to a record rather than a column, but it is still a change, not the
# exploratory whole-record combination that isolation screening produces.
ScidrAnomalyTypeKey <- function(method, kind, field, conditioned, coded = FALSE) {
  if (isTRUE(coded)) "coded" else if (method == "missingness") "missingness" else
    if (method %in% c("change", "trajectory")) "change" else
    if (method == "transition") "transition" else
    if (!isTRUE(field)) "combination" else if (method == "cellwise") "inconsistent" else
    if (identical(kind, "date_interval")) "date_interval" else if (method == "dates") "date" else
    if (method %in% c("change", "trajectory") || identical(kind, "interval")) "change" else
    if (method == "transition") "transition" else
    if (isTRUE(conditioned)) "conditional" else "distribution"
}

ScidrAnomalyTypeFromKey <- function(key) {
  df_Types <- ScidrAnomalyTypes()
  df_Types$Type[match(key, df_Types$Key)]
}

ScidrAnomalyType <- function(method, kind, field, conditioned, coded = FALSE)
  ScidrAnomalyTypeFromKey(ScidrAnomalyTypeKey(method, kind, field, conditioned, coded))

# Severity as the page shows it: one row per tier, naming the types that sit
# there and why that tier outranks the next.
ScidrAnomalySeverityTable <- function() {
  df_Severity <- ScidrAnomalySeverity()
  df_Severity$Type <- ScidrAnomalyTypeFromKey(df_Severity$Key)
  parts <- lapply(split(df_Severity, df_Severity$Severity), function(df_Tier) data.frame(
    Tier = df_Tier$Severity[1], Types = paste(df_Tier$Type, collapse = "; "),
    Meaning = paste(unique(df_Tier$Rationale), collapse = "; "), stringsAsFactors = FALSE))
  df_Table <- do.call(rbind, parts)
  df_Table[order(df_Table$Tier), ]
}

# A value can be so far outside its column that no judgement is needed to see it
# is wrong: a BMI of 276,249 beside a maximum of 48.5. Severity by kind of
# finding cannot express that, because the kind is the same whether the value is
# slightly or absurdly out of range. The test is scale free and robust: the
# distance from the column median, as a multiple of the largest such distance
# among the other records. A value at least this many times further out than
# every other record is promoted ahead of the tiers.
ScidrAnomalyExtremeFactor <- function() 10

# NA when the column cannot support the test: too few comparable records, a
# non-numeric or temporal column, or no spread among the other values.
ScidrAnomalyExtremeRatio <- function(values, index, minimum = 10L) {
  if (!is.numeric(values) || length(index) != 1L || is.na(index)) return(NA_real_)
  value <- values[index]
  if (!is.finite(value)) return(NA_real_)
  others <- values[-index]
  others <- others[is.finite(others)]
  if (length(others) < minimum) return(NA_real_)
  centre <- stats::median(others)
  span <- max(abs(others - centre))
  if (!is.finite(span) || span <= 0) return(NA_real_)
  abs(value - centre) / span
}
