#' Winsorize a numeric vector using SD or IQR thresholds
#'
#' This function performs winsorization on a numeric vector by capping extreme
#' values at calculated lower and upper thresholds. Thresholds can be based on
#' either standard deviation (assuming approximate normality) or interquartile
#' range (robust to skewed distributions).
#'
#' @param Data A numeric vector to be winsorized.
#' @param method Character string specifying the method: "sd" (default) or "iqr".
#' @param sdlim Numeric. Number of standard deviations for the "sd" method.
#' @param iqrlim Numeric. Multiplier for the IQR when method = "iqr" (default 1.5).
#'
#' @return A numeric vector with values winsorized to the specified thresholds.
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
#' @export
windsorize <- function(Data,
                       method = "sd",
                       sdlim = 2.5,
                       iqrlim = 1.5) {

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
