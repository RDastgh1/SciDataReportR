#' Plot Bland-Altman Agreement Plot
#'
#' Generates a Bland-Altman plot to visualize the agreement between two variables.
#'
#' @param data A data frame containing the variables to compare.
#' @param Variable1 The name of the first variable (as a string) to compare.
#' @param Variable2 The name of the second variable (as a string) to compare.
#'
#' @return A list containing:
#'   \item{plot}{A ggplot2 object of the Bland-Altman plot.}
#'   \item{stats}{A list of Bland-Altman statistics from the BlandAltmanLeh package.}
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble rownames_to_column
#'
#'
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @examples
#' # Bland-Altman compares two measurements of the SAME quantity on the same
#' # scale. Here two devices measure the same underlying value, with device B
#' # carrying a small constant bias plus noise.
#' set.seed(101)
#' n <- 80
#' truth <- rnorm(n, mean = 100, sd = 15)
#' method_data <- data.frame(
#'   SampleID = paste0("S", 1:n),
#'   DeviceA  = truth + rnorm(n, 0, 3),
#'   DeviceB  = truth + 2 + rnorm(n, 0, 3)
#' )
#'
#' result <- PlotBlandAltman(method_data, "DeviceA", "DeviceB")
#'
#' # Agreement plot: mean difference (bias) and 95% limits of agreement
#' result$plot
#'
#' # Underlying statistics
#' result$stats$mean.diffs
#' @export
#'
#' @note This function is adapted from code written by Eran Shorer.
PlotBlandAltman <- function(data,
    Variable1,
    Variable2,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "PlotBlandAltman(DataFrame)", "PlotBlandAltman(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data


  # Check if BlandAltmanLeh package is available
  if (!requireNamespace("BlandAltmanLeh", quietly = TRUE)) {
    stop("The package 'BlandAltmanLeh' is required for this function. Please install it using install.packages('BlandAltmanLeh').")
  }

  # Remove rows with missing values
  DataFrame <- na.omit(DataFrame)

  # Calculate Bland-Altman statistics
  BAstats <- BlandAltmanLeh::bland.altman.stats(DataFrame[[Variable1]], DataFrame[[Variable2]])

  # Create a new data frame for plotting
  df <- data.frame(
    Average = BAstats$means,
    lowerlimit = BAstats$lower.limit,
    upperlimit = BAstats$upper.limit,
    Diffs = BAstats$diffs
  )

  # Create caption for the plot
  caption <- paste(
    "Mean Difference =", round(BAstats$lines[2], 2), "\n",
    "Limits of Agreement: [", round(BAstats$lower.limit, 2), ",", round(BAstats$upper.limit, 2), "]"
  )

  # Generate Bland-Altman plot
  p_BA <- ggplot2::ggplot(df, ggplot2::aes(x = Average, y = Diffs)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = BAstats$CI.lines[[1]], ymax = BAstats$CI.lines[[2]]),
      fill = "lightgrey", alpha = 0.3, color = 'white'
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = BAstats$CI.lines[[5]], ymax = BAstats$CI.lines[[6]]),
      fill = "lightgrey", alpha = 0.3, color = 'white'
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = BAstats$CI.lines[[3]], ymax = BAstats$CI.lines[[4]]),
      fill = "lightgray", alpha = 0.3, color = 'white'
    ) +
    ggplot2::geom_point(shape = 16, size = 1.5, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = BAstats$mean.diffs, linetype = "dashed", color = "blue", size = 0.8) +
    ggplot2::geom_hline(yintercept = BAstats$upper.limit, linetype = "dashed", color = "blue", size = 0.7) +
    ggplot2::geom_hline(yintercept = BAstats$lower.limit, linetype = "dashed", color = "blue", size = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(subtitle = caption) +
    ggplot2::xlab("Averages") +
    ggplot2::ylab(paste0(Variable1, "-", Variable2))

  return(list(plot = p_BA, stats = BAstats))
}
