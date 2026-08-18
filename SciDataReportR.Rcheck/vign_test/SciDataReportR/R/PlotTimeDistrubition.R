#' Plot Time Distribution
#'
#' This function plots the distribution of time-based data.
#'
#' @param data The data frame containing the time-based data.
#' @param DateVariable The name of the column in the data frame containing the date information. Default is "Date".
#' @return A ggplot object displaying the distribution of time-based data.
#' @importFrom ggplot2 ggplot geom_boxplot scale_y_reverse theme_linedraw
#' @importFrom dplyr %>%
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   Date = as.Date("2024-01-01") + sample(0:364, 200, replace = TRUE)
#' )
#'
#' PlotTimeDistribution(df, DateVariable = "Date")
#' @export
PlotTimeDistribution <- function(data,
    DateVariable = "Date",
    Data = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotTimeDistribution(Data)", "PlotTimeDistribution(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data


  # Convert date variable to numeric year
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop(
      "Package 'lubridate' is required by PlotTimeDistribution(). ",
      "Install it with install.packages('lubridate')."
    )
  }
  Data$Year <- lubridate::decimal_date(Data[[DateVariable]])

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
