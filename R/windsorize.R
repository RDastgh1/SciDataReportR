#' Windsorize
#'
#' This function performs windsorization on a numeric vector, replacing values below or above a certain threshold
#' with the corresponding threshold value.
#'
#' @param Data A numeric vector to be windsorized.
#' @param sdlim The number of standard deviations to use for the windsorization threshold. Values below
#'              (mean - sdlim * sd) or above (mean + sdlim * sd) will be replaced with the threshold values.
#' @return The windsorized numeric vector.
#' @export
windsorize <- function(Data, sdlim = 2.5) {
  # Calculate mean and standard deviation
  mean_d <- mean(Data, na.rm = TRUE)
  sd_d <- sd(Data, na.rm = TRUE)

  # Calculate windsorization thresholds
  m_min <- mean_d - sdlim * sd_d
  m_max <- mean_d + sdlim * sd_d

  # Replace values below the lower threshold with the lower threshold
  Data[Data < m_min] <- m_min

  # Replace values above the upper threshold with the upper threshold
  Data[Data > m_max] <- m_max

  return(Data)
}
