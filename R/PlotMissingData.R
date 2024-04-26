#' Plot Missing Data
#'
#' A wrapper function to plot missing data using the naniar package and flip coordinates.
#'
#' @param DataFrame The dataframe containing missing data.
#' @param Variables A character vector specifying the variables to be included in the plot. If NULL (default), all variables in the dataframe will be used.
#' @param Relabel Logical indicating whether to relabel the columns based on their labels. Default is TRUE.
#' @return A ggplot object visualizing missing data with flipped coordinates.
#'
#' @import naniar
#' @import ggplot2
#' @importFrom ggplot2 coord_flip
#'
PlotMissingData <- function(DataFrame, Variables = NULL, Relabel = TRUE){
  # If Variables argument is not provided, use all columns of the DataFrame
  if(is.null(Variables)){
    Variables = colnames(DataFrame)
  }

  # Subset DataFrame to include only specified Variables
  DataFrame <- DataFrame[Variables]

  # Relabel columns based on their labels if Relabel is TRUE
  if(Relabel){
    DataFrame <- ReplaceMissingLabels(DataFrame)
    labels <- get_label(DataFrame)
    names(labels) <- NULL
    labels <- make.unique(labels)
    colnames(DataFrame) <- labels
  }

  # Plot missing data and flip coordinates
  pMissing <- naniar::vis_miss(DataFrame) + coord_flip()

  return(pMissing)
}
