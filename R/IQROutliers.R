#' Detect outliers using the Tukey IQR rule and visualize results
#'
#' This function identifies potential outliers in a numeric variable using the
#' Tukey interquartile range (IQR) rule (Tukey, 1977). It returns the full dataset
#' with an outlier flag, a tibble of outlier rows only, a list of thresholds used
#' for classification, and a ggplot boxplot showing the distribution and detected
#' outliers.
#'
#' The outlier rule is defined as:
#'   \deqn{value < Q1 - 1.5 * IQR  \; \textrm{or} \;  value > Q3 + 1.5 * IQR}
#'
#' This approach is widely used in exploratory data analysis for identifying
#' extreme observations that may reflect measurement error, batch effects, or
#' genuine biological variability.
#'
#' @param df A data frame or tibble containing the variable to evaluate.
#' @param variable A string specifying the name of the numeric variable to test.
#' @param id A string specifying the unique identifier column in \code{df}.
#'   Defaults to \code{"ID"}.
#' @param group A string specifying the grouping or batch column to appear
#'   on the x-axis of the diagnostic plot. Defaults to \code{"SourceFile"}.
#'
#' @return A list with four elements:
#'   \itemize{
#'     \item \code{data}: the input dataframe with an added \code{outlier} column.
#'     \item \code{outliers}: a tibble containing only rows flagged as outliers.
#'     \item \code{thresholds}: a named list containing Q1, Q3, IQR, lower bound,
#'       and upper bound used for classification.
#'     \item \code{plot}: a ggplot2 object showing a boxplot and jittered points
#'       colored by outlier status.
#'   }
#'
#' @details
#' The function uses \code{rlang::ensym()} to safely evaluate variable names and
#' generates a ggplot using \code{geom_boxplot()} and \code{geom_jitter()} to allow
#' visual inspection of outlier distribution across groups.
#'
#' Missing values are ignored when computing quartiles. The function checks that
#' the specified \code{variable}, \code{id}, and \code{group} columns exist in the
#' input dataframe.
#'
#' @references
#' Tukey, J. W. (1977). \emph{Exploratory Data Analysis}. Addison-Wesley.
#'
#' @examples
#' \dontrun{
#'   results <- IQROutliers(df = mydata, variable = "RT", id = "SubjectID", group = "Session")
#'   results$plot
#'   results$outliers
#' }
#'
#' @export
IQROutliers <- function(df, variable, id = "ID", group = "SourceFile") {

  if (!variable %in% names(df)) stop("Variable not found in dataframe.")
  if (!id %in% names(df)) stop("ID column not found.")
  if (!group %in% names(df)) stop("Grouping column not found.")

  var <- rlang::ensym(variable)
  grp <- rlang::ensym(group)

  data <- df

  Q1 <- quantile(data[[variable]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[variable]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR

  data$outlier <- data[[variable]] < lower | data[[variable]] > upper

  p <- ggplot2::ggplot(data, ggplot2::aes(x = !!grp, y = !!var, color = outlier)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    ggplot2::labs(
      title = paste("IQR Outlier Detection for", variable),
      y = variable,
      x = "",
      color = "Outlier"
    ) +
    ggplot2::theme_minimal()

  list(
    data = data,
    outliers = dplyr::filter(data, outlier) %>%
      dplyr::select(dplyr::all_of(id), dplyr::all_of(variable),
                    dplyr::all_of(group), outlier),
    thresholds = list(
      Q1 = Q1, Q3 = Q3, IQR = IQR,
      lower = lower, upper = upper
    ),
    plot = p
  )
}
