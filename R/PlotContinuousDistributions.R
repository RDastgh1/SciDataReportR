#' Plot Continuous Distributions
#'
#' Creates rain-cloud plots (half-violin + box/median + scatter) for one or
#' more continuous variables, with optional group-wise colouring.
#'
#' @param data A data frame containing the variables to be plotted.
#' @param variables Character vector of column names to plot.
#' @param Fill Optional column name for grouping.
#' @param Relabel Logical; use variable labels when available.
#' @param FacetLabelStyle One of "both", "label_only", "variable_only", "auto".
#' @param ncol Number of columns in the facet grid.
#' @param TreatOrdinalAs How ordinal variables are handled. `"Continuous"`
#' includes their numeric score; `"Exclude"` omits them. `"Both"` is not
#' meaningful for this plot and errors.
#' @param Ordinal Deprecated logical compatibility option; use
#' `TreatOrdinalAs` instead.
#'
#' @return A ggplot object.
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # Eight variables across three columns show wrapped labels and multi-row facets.
#' PlotContinuousDistributions(
#'   data = Labelled,
#'   variables = c("AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Ferritin",
#'                 "Gamma_Interferon_induced_Monokin", "MMP7", "tau", "p_tau"),
#'   ncol = 3
#' )
#'
#' # Grouped rain-clouds use the Diagnosis fill to compare distributions.
#' PlotContinuousDistributions(
#'   data = Labelled,
#'   variables = c("Ab_42", "p_tau", "tau", "GRO_alpha", "MMP10", "TRAIL_R3"),
#'   Fill = "Diagnosis",
#'   ncol = 3
#' )
#' @details
#' Categorical input safeguards apply to `Fill`.
#' Numeric columns are converted locally to factors, with a warning, only when
#' all observed values are finite whole-number codes and there are at most 10
#' distinct non-missing values. Other numeric inputs error: explicitly convert
#' the column to a factor first if categorical treatment is intended (for
#' example, `df_All$cohort <- factor(df_All$cohort)`). Few unique values alone
#' do not establish categorical meaning. Factors, character columns, and logical
#' columns are accepted; more than 10 observed categories warn and continue.
#' Entirely missing or unsupported columns error. Caller data are unchanged,
#' and existing restrictions on the required number of groups still apply.
#' @export
PlotContinuousDistributions <- function(data,
    variables = NULL,
    Fill = NULL,
    Relabel = TRUE,
    FacetLabelStyle = c("both", "label_only", "variable_only", "auto"),
    ncol = 3,
    Ordinal = lifecycle::deprecated(),
    TreatOrdinalAs = "Categorical",
    DataFrame = lifecycle::deprecated(),
    Variables = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "PlotContinuousDistributions(DataFrame)", "PlotContinuousDistributions(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "PlotContinuousDistributions(Variables)", "PlotContinuousDistributions(variables)")
    variables <- Variables
  }
  Variables <- variables

  if (lifecycle::is_present(Ordinal)) {
    lifecycle::deprecate_warn("20.20.0", "PlotContinuousDistributions(Ordinal)", "PlotContinuousDistributions(TreatOrdinalAs)")
    TreatOrdinalAs <- if (isTRUE(Ordinal)) "Continuous" else "Exclude"
  }
  TreatOrdinalAs <- match.arg(TreatOrdinalAs, c("Categorical", "Continuous", "Both", "Exclude"))
  if (TreatOrdinalAs %in% c("Categorical", "Both")) {
    if (TreatOrdinalAs == "Both") {
      stop("TreatOrdinalAs = 'Both' is not meaningful for PlotContinuousDistributions().", call. = FALSE)
    }
    TreatOrdinalAs <- "Exclude"
  }


  FacetLabelStyle <- match.arg(FacetLabelStyle)

  ScidrValidateVariable(DataFrame, Fill, "Fill", allow_null = TRUE)
  DataFrame <- ScidrPrepareCategorical(DataFrame, Fill, "Fill")
  if (!is.null(Fill) && Fill %in% Variables) {
    stop("`Fill` must not also be included in `variables`.", call. = FALSE)
  }

  # Select Variables
  if (is.null(Variables)) {
    Variables <- getNumVars(DataFrame, Ordinal = FALSE)
    if (TreatOrdinalAs == "Continuous") Variables <- getNumVars(DataFrame, Ordinal = TRUE)
  }

  Variables <- setdiff(Variables, Fill)

  ordinal <- ConvertOrdinalToNumeric(
    DataFrame, Variables, TreatOrdinalAs = TreatOrdinalAs,
    Relabel = Relabel, ReturnMetadata = TRUE
  )
  DataFrame <- ordinal$data
  Variables <- ordinal$variables

  # Build facet labels
  var_labels <- sjlabelled::get_label(DataFrame[Variables], def.value = Variables)

  build_label <- function(var, label) {

    if (!Relabel) return(var)

    if (is.null(label) || is.na(label) || label == "") {
      return(var)
    }

    switch(
      FacetLabelStyle,
      both = paste0(var, "\n", label),
      label_only = label,
      variable_only = var,
      auto = if (identical(var, label)) var else label
    )
  }

  facetlabels <- mapply(build_label, Variables, var_labels, USE.NAMES = FALSE)

  # Pivot Longer
  long_vars <- c(Variables, Fill)

  ContData  <- DataFrame %>%
    dplyr::select(dplyr::all_of(long_vars)) %>%
    # Labelled numeric columns can carry different variable labels, which
    # vctrs correctly refuses to combine in pivot_longer(). Plot values are
    # numeric; labels have already been captured above for the facet text.
    dplyr::mutate(dplyr::across(dplyr::all_of(Variables), as.numeric)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(Variables)) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(Mean = mean(value, na.rm = TRUE)) %>%
    dplyr::ungroup()

  ContData$name <- factor(
    ContData$name,
    levels = Variables,
    labels = facetlabels
  )

  # Fill handling
  if (is.null(Fill)) {
    ContData[[".fill_const"]] <- "1"
    fill_var <- ".fill_const"
    legend   <- FALSE
  } else {
    fill_var <- Fill
    legend   <- TRUE
  }

  # Plot
  p <- ggplot2::ggplot(
    ContData,
    ggplot2::aes(
      y = value,
      x = 1,
      fill   = .data[[fill_var]],
      colour = .data[[fill_var]]
    )
  ) +
    ggrain::geom_rain(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ name, scales = "free", ncol = ncol) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  # Legend handling. With no grouping there is a single series, which keeps its
  # own fixed colour; with a Fill grouping the groups take the package palette.
  if (!legend) {
    p <- p +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::scale_fill_manual(values = "#6EC259") +
      ggplot2::scale_colour_manual(values = "#6EC259")
  } else {
    fill_label <- ScidrDisplayLabels(DataFrame, Fill, Relabel)[[Fill]]
    p <- p + .SciDataFillScale(name = fill_label) + .SciDataColourScale(name = fill_label)
  }

  return(p)
}
