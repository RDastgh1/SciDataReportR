#' Create facet labels for ggplot2 based on variable labels in a data frame
#'
#' This function takes a data frame containing variable labels and creates facet labels
#' suitable for use with ggplot2 facet functions.
#'
#' @param data A data frame containing variable labels.
#'
#' @return A character vector containing facet labels.
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
