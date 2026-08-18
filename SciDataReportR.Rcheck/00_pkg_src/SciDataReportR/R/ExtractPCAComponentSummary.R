#' Extract PCA component summaries
#'
#' Extract variables contributing to each PCA component based on an
#' absolute loading threshold. Returns both a tidy long-format table
#' and compact summary tables suitable for reporting.
#'
#' Negative contributors can optionally be formatted in red HTML text
#' for improved readability in HTML tables.
#'
#' @details
#' Use this after [CreatePCAObject()], at the point where the components exist
#' but do not yet mean anything. PCA hands back axes named `PC1`, `PC2`, `PC3`
#' - each a weighted combination of every input variable - and the analysis
#' cannot be written up until those names are replaced by something
#' interpretable.
#'
#' The full loadings matrix is the wrong tool for that job: with 40 variables
#' and 8 components it is 320 numbers, most of them near zero and irrelevant.
#' This function keeps only the loadings whose absolute value clears
#' `loading_threshold` (0.4 by default), which are the variables actually
#' driving each component, and reports them per component, sorted by strength.
#' Reading down the retained variables is what tells you that PC1 is "an
#' inflammatory axis" or that PC2 separates volume from thickness.
#'
#' The sign matters as much as the magnitude. A component with some variables
#' loading positively and others negatively is a *contrast* - it scores the
#' balance between two sets of measures rather than their overall level - and
#' `html_format = TRUE` prints the negative contributors in red so that
#' structure is visible at a glance instead of having to be hunted for in a
#' column of numbers.
#'
#' Once a component has a name, that name is what belongs on the axis of every
#' downstream plot and in every table of PCA scores.
#'
#' @section Choosing a threshold:
#' `loading_threshold` trades completeness against readability. Raise it
#' (0.5-0.6) when a component retains so many variables that no theme is
#' visible; lower it (0.3) when a component comes back nearly empty, which
#' usually means its variance is spread thinly across many measures rather
#' than concentrated in a few. `top_n` caps the list per component instead,
#' which is the better control when what you want is a compact table for a
#' manuscript rather than a different scientific claim.
#'
#' @seealso [CreatePCAObject()] to fit the PCA, [CreatePCATable()] for the
#'   variance-explained table, and [ProjectPCA()] to score new data on the
#'   components once they have been interpreted.
#'
#' @param PCAObject Output object from CreatePCAObject().
#' @param loading_threshold Minimum absolute loading required for inclusion.
#'   Default is 0.4.
#' @param top_n Optional maximum number of contributors per component.
#'   If NULL, all contributors above threshold are retained.
#' @param use_labels Logical indicating whether variable labels should
#'   be used when available. Default TRUE.
#' @param html_format Logical indicating whether negative contributors
#'   should be formatted using red HTML text. Default TRUE.
#'
#' @return A list containing:
#' \item{LongTable}{
#' A tidy tibble with one row per contributor.
#' }
#' \item{SummaryTable}{
#' A compact tibble with one row per component and comma-separated
#' contributor summaries.
#' }
#' \item{SummaryTableLines}{
#' A compact tibble with one row per component and line-separated
#' contributor summaries.
#' }
#' \item{FormattedSummaryTable}{
#' A formatted gt table with comma-separated contributors.
#' }
#' \item{FormattedSummaryTableLines}{
#' A formatted gt table with line-separated contributors.
#' }
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' vars_Biomarkers <- c(
#'   "AXL", "Adiponectin", "Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
#'   "Apolipoprotein_A1", "Apolipoprotein_B", "C_Reactive_Protein",
#'   "Cortisol", "Cystatin_C", "Ferritin", "Insulin", "Leptin", "p_tau"
#' )
#'
#' # Thirteen correlated biomarkers reduced to a handful of components
#' pca_obj <- CreatePCAObject(
#'   data = df_Labelled,
#'   VarsToReduce = vars_Biomarkers
#' )
#'
#' # At this point the components are called RC1, RC2, RC3
#' summary_obj <- ExtractPCAComponentSummary(pca_obj)
#'
#' # One row per component, listing only the variables that drive it
#' summary_obj$FormattedSummaryTableLines
#'
#' # The same information, one contributor per row
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(
#'     dplyr::mutate(
#'       summary_obj$LongTable,
#'       dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3))
#'     ),
#'     height = "320px", full_width = TRUE
#'   )
#' )))
#'
#' # A stricter threshold
#' ExtractPCAComponentSummary(pca_obj, loading_threshold = 0.6)$SummaryTable
#'
#' # Capping the list per component instead
#' ExtractPCAComponentSummary(pca_obj, top_n = 3)$SummaryTable
#' }
#'
#' @export
ExtractPCAComponentSummary <- function(
    PCAObject,
    loading_threshold = 0.4,
    top_n = NULL,
    use_labels = TRUE,
    html_format = TRUE
) {

  # Validate inputs

  if (!is.list(PCAObject)) {
    stop("PCAObject must be a list returned from CreatePCAObject().")
  }

  if (!"LoadingTable" %in% names(PCAObject)) {
    stop("PCAObject does not contain LoadingTable.")
  }

  LoadingTable <- PCAObject$LoadingTable

  rc_cols <- names(LoadingTable)[grepl("^RC", names(LoadingTable))]

  if (length(rc_cols) == 0) {
    stop("No PCA component columns beginning with 'RC' were found.")
  }

  if (!"Variable" %in% names(LoadingTable)) {
    stop("LoadingTable must contain a Variable column.")
  }

  # Prepare data

  LongTable <- LoadingTable %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(rc_cols),
      names_to = "Component",
      values_to = "Loading"
    ) %>%
    dplyr::mutate(
      AbsLoading = abs(Loading),
      Direction = dplyr::if_else(
        Loading >= 0,
        "Positive",
        "Negative"
      )
    ) %>%
    dplyr::filter(
      AbsLoading >= loading_threshold
    )

  # Apply labels

  if (use_labels && "Labels" %in% names(LongTable)) {

    LongTable <- LongTable %>%
      dplyr::mutate(
        Label = dplyr::if_else(
          is.na(Labels) | Labels == "",
          Variable,
          Labels
        )
      )

  } else {

    LongTable <- LongTable %>%
      dplyr::mutate(
        Label = Variable
      )

  }

  # Remove old Labels column if present

  if ("Labels" %in% names(LongTable)) {

    LongTable <- LongTable %>%
      dplyr::select(-Labels)

  }

  # Rank contributors

  LongTable <- LongTable %>%
    dplyr::group_by(Component) %>%
    dplyr::arrange(
      dplyr::desc(AbsLoading),
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      Rank = dplyr::row_number()
    )

  # Apply optional top_n filter

  if (!is.null(top_n)) {

    LongTable <- LongTable %>%
      dplyr::filter(Rank <= top_n)

  }

  LongTable <- LongTable %>%
    dplyr::ungroup()

  # Build formatted contributor strings

  ContributorTable <- LongTable %>%
    dplyr::mutate(

      ContributorText = dplyr::case_when(

        Direction == "Negative" & html_format ~
          paste0(
            "<span style='color:red'>",
            Label,
            "</span>"
          ),

        Direction == "Negative" & !html_format ~
          paste0("-", Label),

        TRUE ~
          Label
      )

    )

  # Build comma-separated summary table

  SummaryTable <- ContributorTable %>%
    dplyr::group_by(Component) %>%
    dplyr::summarise(
      Contributors = paste(
        ContributorText,
        collapse = ", "
      ),
      .groups = "drop"
    )

  # Build line-separated summary table

  SummaryTableLines <- ContributorTable %>%
    dplyr::group_by(Component) %>%
    dplyr::summarise(
      Contributors = paste(
        ContributorText,
        collapse = "<br>"
      ),
      .groups = "drop"
    )

  # Build formatted gt table (comma separated)

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop(
      "Package 'gt' is required by ExtractPCAComponentSummary(). ",
      "Install it with install.packages('gt')."
    )
  }

  FormattedSummaryTable <- SummaryTable %>%
    gt::gt() %>%
    gt::fmt_markdown(
      columns = "Contributors"
    ) %>%
    gt::cols_label(
      Component = "Component",
      Contributors = "Top Contributors"
    ) %>%
    gt::tab_options(
      table.width = gt::pct(100),
      data_row.padding = gt::px(4)
    )

  # Build formatted gt table (line separated)

  FormattedSummaryTableLines <- SummaryTableLines %>%
    gt::gt() %>%
    gt::fmt_markdown(
      columns = "Contributors"
    ) %>%
    gt::cols_label(
      Component = "Component",
      Contributors = "Top Contributors"
    ) %>%
    gt::tab_options(
      table.width = gt::pct(100),
      data_row.padding = gt::px(4)
    )

  # Return result

  return(list(
    LongTable = LongTable,
    SummaryTable = SummaryTable,
    SummaryTableLines = SummaryTableLines,
    FormattedSummaryTable = FormattedSummaryTable,
    FormattedSummaryTableLines = FormattedSummaryTableLines
  ))

}
