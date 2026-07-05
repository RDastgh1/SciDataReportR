#' Winsorize a numeric vector using SD or IQR thresholds
#'
#' This function performs winsorization on a numeric vector by capping extreme
#' values at calculated lower and upper thresholds. Thresholds can be based on
#' either standard deviation (assuming approximate normality) or interquartile
#' range (robust to skewed distributions).
#'
#' @param data A numeric vector to be winsorized.
#' @param method Character string specifying the method: "sd" (default) or "iqr".
#' @param sdlim Numeric. Number of standard deviations for the "sd" method.
#' @param iqrlim Numeric. Multiplier for the IQR when method = "iqr" (default 1.5).
#'
#' @return A numeric vector with values winsorized to the specified thresholds.
#'
#' @details
#' Winsorization replaces extreme values with the nearest threshold value rather
#' than discarding them, limiting the influence of outliers while retaining the
#' sample size. The term originates with the statistician Charles P. Winsor; see
#' Dixon (1960) and Tukey (1962) for early treatments of the estimator.
#'
#' @references
#' Tukey, J. W. (1962). The future of data analysis. \emph{The Annals of
#' Mathematical Statistics}, 33(1), 1-67. \doi{10.1214/aoms/1177704711}
#'
#' Dixon, W. J. (1960). Simplified estimation from censored normal samples.
#' \emph{The Annals of Mathematical Statistics}, 31(2), 385-391.
#' \doi{10.1214/aoms/1177705900}
#'
#' @examples
#' x <- c(rnorm(100), 10, 15, -12)
#'
#' # SD-based winsorization
#' windsorize(x, method = "sd", sdlim = 2.5)
#'
#' # IQR-based winsorization
#' windsorize(x, method = "iqr", iqrlim = 1.5)
#'
#' # Compare the distribution before and after winsorization
#' set.seed(42)
#' x <- c(rnorm(200, mean = 10, sd = 2), 30, 32, -8, -10)
#' compare_df <- data.frame(
#'   raw        = x,
#'   winsorized = windsorize(x, method = "iqr", iqrlim = 1.5)
#' )
#'
#' PlotContinuousDistributions(compare_df, variables = c("raw", "winsorized"))
#'
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
windsorize <- function(data,
    sdlim = 2.5,
    iqrlim = 1.5,
    method = "sd",
    Data = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "windsorize(Data)", "windsorize(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data


  # Validate inputs
  if (!is.numeric(Data)) {
    stop("Data must be a numeric vector.")
  }

  if (!method %in% c("sd", "iqr")) {
    stop("method must be either 'sd' or 'iqr'.")
  }

  # Remove NA for calculations (but preserve positions)
  valid_data <- Data[!is.na(Data)]

  # Calculate thresholds
  if (method == "sd") {
    mean_d <- mean(valid_data)
    sd_d <- sd(valid_data)

    lower <- mean_d - sdlim * sd_d
    upper <- mean_d + sdlim * sd_d

  } else if (method == "iqr") {
    q1 <- stats::quantile(valid_data, 0.25, names = FALSE)
    q3 <- stats::quantile(valid_data, 0.75, names = FALSE)
    iqr <- q3 - q1

    lower <- q1 - iqrlim * iqr
    upper <- q3 + iqrlim * iqr
  }

  # Apply winsorization
  Data[Data < lower] <- lower
  Data[Data > upper] <- upper

  return(Data)
}
