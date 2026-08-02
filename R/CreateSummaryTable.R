
#' Create Summary Table
#'
#' Generate a descriptive summary table for specified variables in a dataset.
#'
#' @param data The dataset containing the variables of interest.
#' @param variables A character vector specifying the variables for which summary statistics will be calculated.
#' @param digits Number of decimal places to round the summary statistics.
#' @param Relabel Logical, indicating whether to use variable labels as column headers.
#' @param TreatOrdinalAs How ordinal variables are handled. This numeric
#' descriptive table accepts `"Continuous"` or `"Exclude"`.
#' @param Ordinal Deprecated logical compatibility option; use
#' `TreatOrdinalAs` instead.
#' @param ScrollBoxHeight Height of the scroll box for displaying the table.
#' @return A formatted HTML table displaying summary statistics.
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate if_else
#' @importFrom tibble rownames_to_column
#' @importFrom sjlabelled get_label set_label
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @param numdecimals \strong{Deprecated} (since 19.15.0). Use \code{digits} instead.
#' @export


CreateSummaryTable <- function(data,
    variables = NULL,
    digits = 2,
    Relabel = TRUE,
    Ordinal = lifecycle::deprecated(),
    TreatOrdinalAs = "Categorical",
    ScrollBoxHeight = "700px",
    Data = lifecycle::deprecated(),
    Variables = lifecycle::deprecated(),
    numdecimals = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "CreateSummaryTable(Data)", "CreateSummaryTable(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "CreateSummaryTable(Variables)", "CreateSummaryTable(variables)")
    variables <- Variables
  }
  Variables <- variables
  if (lifecycle::is_present(numdecimals)) {
    lifecycle::deprecate_warn("19.15.0", "CreateSummaryTable(numdecimals)", "CreateSummaryTable(digits)")
    digits <- numdecimals
  }
  numdecimals <- digits

  if (lifecycle::is_present(Ordinal)) {
    lifecycle::deprecate_warn("20.20.0", "CreateSummaryTable(Ordinal)", "CreateSummaryTable(TreatOrdinalAs)")
    TreatOrdinalAs <- if (isTRUE(Ordinal)) "Continuous" else "Exclude"
  }
  TreatOrdinalAs <- ScidrMatchOrdinalTreatment(TreatOrdinalAs)
  check_vars <- if (is.null(Variables)) names(Data) else intersect(Variables, names(Data))
  has_ordinal <- any(vapply(Data[check_vars], ScidrIsOrdinal, logical(1)))
  if (has_ordinal && TreatOrdinalAs %in% c("Categorical", "Both")) {
    stop("CreateSummaryTable() is numeric-only; use TreatOrdinalAs = 'Continuous' or 'Exclude'. Use MakeTable1() for categorical or both ordinal summaries.", call. = FALSE)
  }

  for (pkg in c("summarytools", "kableExtra")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required by CreateSummaryTable(). ",
        "Install it with install.packages('", pkg, "')."
      )
    }
  }

  if (is.null(Variables)) {
    Variables <- colnames(Data)
  }

  suppressWarnings({
    Data <- dplyr::select(Data, dplyr::all_of(Variables))
    prep <- ScidrPrepareOrdinal(Data, Variables, TreatOrdinalAs)
    Data <- prep$data
    Variables <- prep$variables

    d <- summarytools::descr(Data)
    statVars <- c("Mean", "Std.Dev", "Median", "IQR", "Min", "Max", "Skewness", "Kurtosis", "N.Valid", "Pct.Valid")
    d2 <- as.data.frame(t(as.data.frame(d)))
    d2 <- dplyr::select(d2, dplyr::all_of(statVars))
    d2 <- d2[Variables, ]
    d2[statVars] <- lapply(d2[statVars], round, numdecimals)

    if (Relabel) {
      labels <- ScidrDisplayLabels(Data, Variables, Relabel) %>%
        as.data.frame() %>% tibble::rownames_to_column()
      colnames(labels) <- c("Variable", "label")
      d2$label <- labels$label
      d2 <- dplyr::select(d2, label, dplyr::all_of(statVars))
    }

    SummaryTable <- d2 %>%
      tibble::rownames_to_column("Variable") %>%
      dplyr::mutate(
        Skewness = kableExtra::cell_spec(
          round(Skewness, numdecimals), "html",
          background = dplyr::if_else(abs(Skewness) > 10, "yellow", "", missing = "grey")
        ),
        Kurtosis = kableExtra::cell_spec(
          round(Kurtosis, numdecimals), "html",
          background = dplyr::if_else(abs(Kurtosis) > 10, "yellow", "", missing = "grey")
        ),
        IQR = kableExtra::cell_spec(
          round(IQR, numdecimals), "html",
          background = dplyr::if_else(abs((IQR/Std.Dev/1.34) - 1) > 0.5, "yellow", "", missing = "grey")
        ),
        Pct.Valid = kableExtra::cell_spec(
          round(Pct.Valid, numdecimals), "html",
          background = dplyr::if_else(Pct.Valid < 70, "red", "", missing = "grey")
        )
      ) %>%
      kableExtra::kable(
        format = "html", escape = FALSE, digits = numdecimals,
        row.names = TRUE,
        caption = "Descriptive Summary Table. IQR, Skewness, and Kurtosis are highlighted in yellow if they are indicative of a non-normal distribution. Pct.Valid is highlighted in red if over 30% of data is missing"
      ) %>%
      kableExtra::kable_styling(
        bootstrap_options = c("striped", "hover", "condensed", "responsive")
      ) %>%
      kableExtra::scroll_box(width = "100%", height = ScrollBoxHeight)
  })

  return(SummaryTable)
}
