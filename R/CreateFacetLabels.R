#' Create facet labels for ggplot2 based on variable labels in a data frame
#'
#' This function takes a data frame containing variable labels and creates facet labels
#' suitable for use with ggplot2 facet functions.
#'
#' @param DataFrame A data frame containing variable labels.
#'
#' @return A character vector containing facet labels.
#'
#' @importFrom sjlabelled get_label
#' @import dplyr
#' @import tidyr
#'

createFacetLabels <- function(DataFrame) {
  # Extract labels from DataFrame
  l <- sjlabelled::get_label(DataFrame) %>% as.data.frame() %>% rownames_to_column()

  # Concatenate labels with new line separator
  facetlabels <- do.call(paste0, c(l[1], sep = "\n", l[2]))

  return(facetlabels)
}
