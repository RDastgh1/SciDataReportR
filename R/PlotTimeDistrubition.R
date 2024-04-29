#' Plot Time Distribution
#'
#' This function plots the distribution of time-based data.
#'
#' @param Data The data frame containing the time-based data.
#' @param DateVariable The name of the column in the data frame containing the date information. Default is "Date".
#' @return A ggplot object displaying the distribution of time-based data.
#' @importFrom ggplot2 ggplot geom_boxplot theme_linedraw
#' @importFrom ggdist stat_halfeye stat_dots
#' @importFrom lubridate decimal_date
#' @importFrom scales scale_y_reverse
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @import ggdist
#' @export
PlotTimeDistribution <- function(Data, DateVariable = "Date") {

  # Convert date variable to numeric year
  Data$Year <- lubridate::decimal_date(.data[[DateVariable]])

  # Create ggplot object
  pTimeline <- Data %>%
    ggplot(aes(x = 1, y = Year)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, point_color = NA) +
    geom_boxplot(width = 0.12) +
    ggdist::stat_dots(side = "left", justification = 1.1) +
    scale_y_reverse() +
    theme_linedraw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  return(pTimeline)
}
