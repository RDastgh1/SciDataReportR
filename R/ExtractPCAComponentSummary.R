#' Extract PCA component summaries
#'
#' Extract variables contributing to each PCA component based on an
#' absolute loading threshold. Returns both a tidy long-format table
#' and compact summary tables suitable for reporting.
#'
#' Negative contributors can optionally be formatted in red HTML text
#' for improved readability in HTML tables.
#'
#' @param PCAObject Output object from CreatePCATable().
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
#' \dontrun{
#' pca_obj <- CreatePCATable(
#'   Data = mtcars,
#'   VarsToReduce = colnames(mtcars)
#' )
#'
#' summary_obj <- ExtractPCAComponentSummary(pca_obj)
#'
#' summary_obj$LongTable
#' summary_obj$FormattedSummaryTable
#' summary_obj$FormattedSummaryTableLines
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
    stop("PCAObject must be a list returned from CreatePCATable().")
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
