#' Create facet labels for ggplot2 based on variable labels in a data frame
#'
#' This function takes a data frame containing variable labels and creates facet labels
#' suitable for use with ggplot2 facet functions.
#'
#' @param data A data frame containing variable labels.
#'
#' @return A character vector containing facet labels.
#'
#' @details
#' Faceted plots label their strips with the column name, which in study data
#' is usually the least readable thing available: `p_tau`, `ACE_CD143`,
#' `Weight_kg_v2`. This builds strip labels that carry the column name *and*
#' the variable label attached by [RevalueData()], one per line, so a reader
#' can see what a panel measures without also losing the name needed to find
#' that column again.
#'
#' The returned vector is in column order, which is what makes it usable
#' directly as the `labels` argument of `factor()` when reshaping to long
#' format for faceting.
#'
#' Variables with no attached label fall back to their name, so the vector is
#' always the right length and mixed labelled/unlabelled frames are safe.
#'
#' @seealso [RevalueData()], which attaches the labels this reads, and
#'   [PlotContinuousDistributions()], which applies the same idea internally
#'   through its `FacetLabelStyle` argument.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' vars_Show <- c("AXL", "Adiponectin", "tau", "p_tau")
#'
#' # Column name on the first line, attached label on the second
#' labels_Facet <- createFacetLabels(Labelled[vars_Show])
#' labels_Facet
#'
#' \donttest{
#' df_Long <- tidyr::pivot_longer(
#'   dplyr::mutate(
#'     Labelled[vars_Show],
#'     dplyr::across(dplyr::everything(), as.numeric)
#'   ),
#'   cols = dplyr::everything()
#' )
#'
#' # Without the labels: strips carry bare column names
#' ggplot2::ggplot(df_Long, ggplot2::aes(x = value)) +
#'   ggplot2::geom_histogram(bins = 30, na.rm = TRUE) +
#'   ggplot2::facet_wrap(~ name, scales = "free") +
#'   ggplot2::theme_bw()
#'
#' # With them: each panel says what it measures
#' df_Long$name <- factor(
#'   df_Long$name,
#'   levels = vars_Show,
#'   labels = labels_Facet
#' )
#'
#' ggplot2::ggplot(df_Long, ggplot2::aes(x = value)) +
#'   ggplot2::geom_histogram(bins = 30, na.rm = TRUE) +
#'   ggplot2::facet_wrap(~ name, scales = "free") +
#'   ggplot2::theme_bw()
#' }
#'
#' @importFrom sjlabelled get_label
#' @import dplyr
#' @import tidyr
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
createFacetLabels <- function(data,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "createFacetLabels(DataFrame)", "createFacetLabels(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data

  # Extract labels from DataFrame
  l <- sjlabelled::get_label(DataFrame) %>% as.data.frame() %>% rownames_to_column()

  # Concatenate labels with new line separator
  facetlabels <- do.call(paste0, c(l[1], sep = "\n", l[2]))

  return(facetlabels)
}
