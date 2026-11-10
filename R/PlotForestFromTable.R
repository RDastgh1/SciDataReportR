#' Create a Forest Plot from Univariate Regression Tables
#'
#' This function generates a forest plot from the results of
#' [MakeUnivariateRegressionTable()].
#'
#' @param UnivariateRegressionTables Either the full list returned by
#'   [MakeUnivariateRegressionTable()] (its `Results` dataframe is used
#'   directly), or a dataframe with the `Results` columns. Passing a dataframe
#'   lets you filter, reorder, or relabel `Results` before plotting; required
#'   columns are `OutcomeLabel`, `TermLabel`, `Estimate`, `ConfLow`,
#'   `ConfHigh`, and `PValue` (`Significant` and `ReferenceValue` are
#'   recomputed if absent). Lists created by older package versions (without a
#'   `Results` element) are still supported.
#' @param PVar Significance column: "PValue" (default) or "PValueFDR".
#' @param Pthresh Significance cutoff, default 0.05, applied per estimate.
#' @param show_interaction_p Annotate interaction p-values when Interaction
#'   rows are present; simple slopes are displayed side by side by level.
#' @param facet_scales Facet scales: "fixed", "free", "free_x", or "free_y".
#' @param xlim Optional two-element x-axis display range; estimates are retained.
#' @section Grouped inputs:
#' Group levels are dodged within rows and coloured with SciDataPalette;
#' nonsignificant estimates have reduced opacity. A named list of Results
#' tables or regression outputs is bound using list names as Group labels.
#' Facet labels wrap at 25 characters.
#' @param pSize Numeric. Size of the points in the plot. Default is 2.
#' @param Flip Logical. If `FALSE`, outcomes are facets and predictors/terms
#'   are rows. If `TRUE`, predictors/terms are facets and outcomes are rows.
#' @section Reading the plot:
#' The forest plot earns its place when several predictors are screened against
#' several outcomes at once: one panel per outcome, the same predictors down
#' every panel, so a predictor that matters for one outcome and not the others
#' is visible in a single glance.
#'
#' `Flip = TRUE` swaps the roles - one panel per predictor, outcomes down the
#' rows. Use it to ask "what does age predict?" rather than "what predicts
#' tau?": the same estimates, organized around the other question.
#'
#' Passing the `Results` dataframe instead of the whole object means the plot
#' can be filtered first, for instance to the associations that survive FDR
#' correction across the whole screen. Standardized estimates put every
#' predictor on the same axis, which is what makes the widths of the intervals
#' comparable across rows. A binary outcome gives odds ratios, and the
#' reference line moves from 0 to 1.
#'
#' @return A ggplot object representing the forest plot.
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Several predictors screened against several outcomes
#' urt <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = c("AXL", "tau", "p_tau", "Ferritin"),
#'   predictor_vars = c("age", "sex", "Adiponectin", "Cortisol", "Insulin")
#' )
#'
#' PlotForestFromTable(urt)
#'
#' # One panel per predictor instead
#' PlotForestFromTable(urt, Flip = TRUE)
#'
#' # Filtered to associations surviving FDR correction
#' urt$Results$PValueFDR <- ApplyFDRCorrection(urt$Results$PValue)
#' PlotForestFromTable(urt$Results[urt$Results$PValueFDR < 0.05, ])
#'
#' # Standardized estimates
#' urt_Std <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = c("AXL", "tau", "p_tau", "Ferritin"),
#'   predictor_vars = c("age", "sex", "Adiponectin", "Cortisol", "Insulin"),
#'   Standardize = TRUE
#' )
#' PlotForestFromTable(urt_Std)
#'
#' # A binary outcome: odds ratios
#' urt_Logistic <- MakeUnivariateRegressionTable(
#'   data = Labelled,
#'   outcome_vars = "Diagnosis",
#'   predictor_vars = c("age", "sex", "AXL", "tau", "p_tau")
#' )
#' PlotForestFromTable(urt_Logistic)
#' }
#' @export
PlotForestFromTable <- function(UnivariateRegressionTables, pSize = 2, Flip = FALSE,
    PVar = "PValue", Pthresh = 0.05, show_interaction_p = FALSE,
    facet_scales = c("fixed", "free", "free_x", "free_y"), xlim = NULL) {

  if (!is.logical(Flip) || length(Flip) != 1 || is.na(Flip)) {
    stop("Flip must be TRUE or FALSE.")
  }

  facet_scales <- match.arg(facet_scales)
  if (!PVar %in% c("PValue", "PValueFDR")) stop("PVar must be PValue or PValueFDR.")
  if (length(Pthresh) != 1L || !is.finite(Pthresh) || Pthresh <= 0 || Pthresh >= 1) stop("Invalid Pthresh.")
  if (!is.null(xlim) && (length(xlim) != 2L || any(!is.finite(xlim)) || xlim[1] >= xlim[2])) stop("xlim must be two increasing finite values.")
  df_Combined <- ScidrForestPlotResolveResults(UnivariateRegressionTables)
  interaction_rows <- if ("Term" %in% names(df_Combined)) df_Combined$Term == "Interaction" else rep(FALSE, nrow(df_Combined))
  df_Interactions <- df_Combined[which(interaction_rows), , drop = FALSE]
  if (show_interaction_p && nrow(df_Interactions)) {
    df_Combined <- df_Combined[which(!interaction_rows), , drop = FALSE]
    if (all(c("Predictor", "Level") %in% names(df_Combined))) {
      df_Combined$TermLabel <- if ("PredictorLabel" %in% names(df_Combined)) df_Combined$PredictorLabel else df_Combined$Predictor
      df_Combined$Group <- df_Combined$Level
    }
  }
  df_Combined <- df_Combined %>% dplyr::filter(!is.na(Estimate))
  if (!nrow(df_Combined)) stop("No estimates available to plot.")
  if (!PVar %in% names(df_Combined)) stop("Requested p-value column is unavailable: ", PVar)
  df_Combined$Significant <- !is.na(df_Combined[[PVar]]) & df_Combined[[PVar]] <= Pthresh
  grouped <- "Group" %in% names(df_Combined) && length(unique(stats::na.omit(df_Combined$Group))) > 1L
  if (grouped) {
    df_Combined$Group <- factor(df_Combined$Group, levels = unique(df_Combined$Group))
  }
  if (!"ReferenceValue" %in% names(df_Combined)) {
    df_Combined$ReferenceValue <- if ("EffectType" %in% names(df_Combined)) {
      ifelse(df_Combined$EffectType == "Odds ratio", 1, 0)
    } else {
      0
    }
  }

  outcome_levels <- unique(df_Combined$OutcomeLabel)
  term_levels <- unique(df_Combined$TermLabel)

  if (Flip) {
    df_Combined$PlotRow <- factor(
      df_Combined$OutcomeLabel,
      levels = rev(outcome_levels)
    )
    df_Combined$PlotFacet <- factor(
      df_Combined$TermLabel,
      levels = term_levels
    )
    y_label <- "Outcome"
  } else {
    df_Combined$PlotRow <- factor(
      df_Combined$TermLabel,
      levels = rev(term_levels)
    )
    df_Combined$PlotFacet <- factor(
      df_Combined$OutcomeLabel,
      levels = outcome_levels
    )
    y_label <- "Variable"
  }
  df_ReferenceLines <- df_Combined %>%
    dplyr::select(PlotFacet, ReferenceValue) %>%
    distinct()

  x_label <- if (all(df_Combined$ReferenceValue == 1)) {
    "Odds ratio"
  } else if (all(df_Combined$ReferenceValue == 0)) {
    "Estimate"
  } else {
    "Estimate / odds ratio"
  }

  # Create the forest plot
  p <- df_Combined %>%
    ggplot(aes(x = Estimate, y = PlotRow, color = Significant)) +
    geom_point(size = pSize) +  # Plot the point for the estimate
    geom_errorbar(aes(xmin = ConfLow, xmax = ConfHigh), orientation = "y", width = 0.2) +  # Add the error bars
    facet_wrap(~PlotFacet, nrow = 1, scales = facet_scales, labeller = ggplot2::label_wrap_gen(width = 25)) +  # Facet by outcome or term
    geom_vline(
      data = df_ReferenceLines,
      aes(xintercept = ReferenceValue),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank()
    ) +
    scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "black")) +
    theme(legend.position = "none")

  if (grouped) {
    dodge <- ggplot2::position_dodge(width = 0.65, orientation = "y")
    p$mapping <- ggplot2::aes(x = Estimate, y = PlotRow, colour = Group, group = Group, alpha = Significant)
    p$layers[[1]]$position <- dodge
    p$layers[[2]]$position <- dodge
    p <- suppressMessages(p + ggplot2::scale_colour_manual(values = stats::setNames(
      SciDataPalette(nlevels(df_Combined$Group)), levels(df_Combined$Group)))) +
      ggplot2::scale_alpha_manual(values = c("FALSE" = 0.35, "TRUE" = 1), name = "Significant") +
      ggplot2::theme(legend.position = "right")
  }
  if (show_interaction_p && nrow(df_Interactions)) {
    if (!PVar %in% names(df_Interactions)) stop("Interaction p-value column is unavailable.")
    key <- function(df_Input) paste(df_Input$Outcome, df_Input$Predictor, sep = "\r")
    df_Notes <- df_Combined[!duplicated(paste(df_Combined$PlotRow, df_Combined$PlotFacet)), , drop = FALSE]
    selected <- match(key(df_Notes), key(df_Interactions))
    df_Notes$Note <- paste0("Interaction p = ", format.pval(df_Interactions[[PVar]][selected], digits = 2))
    p <- p + ggplot2::geom_text(data = df_Notes,
      ggplot2::aes(x = Inf, y = PlotRow, label = Note), inherit.aes = FALSE,
      hjust = 1.02, vjust = -1, size = 3)
  }
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  return(p)
}

#' @description `plotForestFromTable()` was renamed to `PlotForestFromTable()`
#'   in SciDataReportR 20.5.0 to match the package's `Plot*` naming
#'   convention. It remains available as a backwards-compatible synonym.
#' @rdname PlotForestFromTable
#' @export
plotForestFromTable <- function(UnivariateRegressionTables, pSize = 2, Flip = FALSE) {
  lifecycle::deprecate_soft("20.5.0", "plotForestFromTable()", "PlotForestFromTable()")
  PlotForestFromTable(UnivariateRegressionTables, pSize = pSize, Flip = Flip)
}

ScidrForestPlotResolveResults <- function(x) {
  required_cols <- c("OutcomeLabel", "TermLabel", "Estimate", "ConfLow", "ConfHigh", "PValue")

  if (is.data.frame(x)) {
    missing_cols <- setdiff(required_cols, names(x))
    if (length(missing_cols) > 0) {
      stop(
        "The results dataframe is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ". Pass the Results dataframe from MakeUnivariateRegressionTable()."
      )
    }
    return(x)
  }

  if (is.list(x)) {
    if (is.data.frame(x$Results)) {
      return(x$Results)
    }
    if (length(x) && !is.null(names(x)) && all(nzchar(names(x))) &&
        all(vapply(x, function(item) is.data.frame(item) || (is.list(item) && is.data.frame(item$Results)), logical(1)))) {
      return(dplyr::bind_rows(lapply(names(x), function(group) {
        df_Group <- ScidrForestPlotResolveResults(x[[group]])
        df_Group$Group <- group
        df_Group
      })))
    }
    if (!is.null(x$FormattedTable)) {
      # Objects saved by versions before 20.5.0 have no Results element;
      # rebuild it from the gtsummary tables.
      return(ScidrForestPlotLegacyResults(x))
    }
  }

  stop(
    "UnivariateRegressionTables must be the output of MakeUnivariateRegressionTable() ",
    "or a dataframe with its Results columns."
  )
}

ScidrForestPlotLegacyResults <- function(UnivariateRegressionTables) {
  list_Tables <- UnivariateRegressionTables$FormattedTable$tbls
  title_Tables <- names(list_Tables)
  if (is.null(title_Tables) || any(is.na(title_Tables)) || any(title_Tables == "")) {
    title_Tables <- paste0("Outcome ", seq_along(list_Tables))
  }
  label_Tables <- title_Tables
  spanner_labels <- ScidrForestPlotSpannerLabels(
    UnivariateRegressionTables$FormattedTable,
    length(list_Tables)
  )
  if (all(!is.na(spanner_labels)) && all(spanner_labels != "")) {
    label_Tables <- spanner_labels
  }
  outcome_metadata <- UnivariateRegressionTables$Metadata$Outcomes
  if (!is.null(outcome_metadata) &&
      all(c("Outcome", "OutcomeLabel") %in% names(outcome_metadata))) {
    outcome_label_lookup <- stats::setNames(
      as.character(outcome_metadata$OutcomeLabel),
      as.character(outcome_metadata$Outcome)
    )
    matched_labels <- outcome_label_lookup[title_Tables]
    use_metadata_label <- !is.na(matched_labels) &
      matched_labels != "" &
      (matched_labels != title_Tables | label_Tables == title_Tables)
    label_Tables <- ifelse(
      use_metadata_label,
      matched_labels,
      label_Tables
    )
  }

  # Combine all tables into a single dataframe
  for (i in seq_along(list_Tables)) {
    d <- list_Tables[[i]]$table_body
    d$Header <- title_Tables[i]
    d$HeaderLabel <- label_Tables[i]
    if (i == 1) {
      df_Combined <- d
    } else {
      df_Combined <- rbind(df_Combined, d)
    }
  }

  df_Combined <- df_Combined %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(
      TermLabel = if_else(
        !is.na(label) & label != var_label,
        paste0(var_label, " : ", label), # Combine var_label and label
        var_label # Keep var_label unchanged for other cases
      ),
      OutcomeLabel = HeaderLabel
    )
  reference_value <- ifelse(
    "effect_type" %in% names(df_Combined) & df_Combined$effect_type == "Odds ratio",
    1,
    0
  )

  data.frame(
    Outcome = df_Combined$Header,
    OutcomeLabel = df_Combined$OutcomeLabel,
    TermLabel = df_Combined$TermLabel,
    Estimate = df_Combined$estimate,
    ConfLow = df_Combined$conf.low,
    ConfHigh = df_Combined$conf.high,
    PValue = df_Combined$p.value,
    Significant = df_Combined$p.value < 0.05,
    ReferenceValue = reference_value,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

ScidrForestPlotSpannerLabels <- function(gtsummary_table, n_tables) {
  out <- rep(NA_character_, n_tables)
  spanning_header <- gtsummary_table$table_styling$spanning_header
  if (is.null(spanning_header) ||
      !all(c("column", "spanning_header") %in% names(spanning_header))) {
    return(out)
  }

  table_index <- suppressWarnings(as.integer(sub(".*_([0-9]+)$", "\\1", spanning_header$column)))
  valid_rows <- !is.na(table_index) &
    table_index >= 1 &
    table_index <= n_tables &
    !is.na(spanning_header$spanning_header) &
    spanning_header$spanning_header != ""

  if (!any(valid_rows)) {
    return(out)
  }

  label_data <- data.frame(
    TableIndex = table_index[valid_rows],
    Label = as.character(spanning_header$spanning_header[valid_rows]),
    stringsAsFactors = FALSE
  )
  label_data <- label_data[!duplicated(label_data$TableIndex), , drop = FALSE]
  out[label_data$TableIndex] <- label_data$Label
  out
}
