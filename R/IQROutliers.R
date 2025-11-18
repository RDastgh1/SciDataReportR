#' Detect outliers using the Tukey IQR rule and visualize results
#'
#' This function identifies potential outliers in a numeric variable using the
#' Tukey interquartile range (IQR) rule (Tukey, 1977). It returns a tibble of the
#' detected outlier rows and a ggplot visualization showing the variable across
#' groups with outlier points highlighted.
#'
#' The outlier rule is defined as:
#'   \deqn{value < Q1 - 1.5 * IQR  \; \textrm{or} \;  value > Q3 + 1.5 * IQR}
#'
#' Outliers are visually highlighted using jittered points colored red, while the
#' boxplot remains uncolored to prevent creation of a separate outlier-only box.
#'
#' @param df A data frame or tibble containing the variable to evaluate.
#' @param Variable A string specifying the name of the numeric variable to test.
#' @param id A string specifying the identifier column to include in the returned
#'   outlier table. Defaults to \code{"ID"}.
#' @param group A string specifying the grouping or batch column to use on the
#'   x-axis of the diagnostic plot. Defaults to \code{"SourceFile"}.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{outlierdf}: a tibble containing only rows flagged as outliers,
#'       including the ID, variable, group, and outlier flag.
#'     \item \code{p}: a ggplot2 object showing a boxplot and jittered points
#'       colored by outlier status.
#'   }
#'
#' @details
#' The color aesthetic is mapped only within \code{geom_jitter()}, ensuring the
#' boxplot is drawn once per group rather than once per outlier class. Missing
#' values are ignored when computing quartiles. Columns specified by \code{id}
#' and \code{group} must be present in the input data frame.
#'
#' @references
#' Tukey, J. W. (1977). \emph{Exploratory Data Analysis}. Addison-Wesley.
#'
#' @examples
#' \dontrun{
#'   IQROutliers2(df_Revalued_Data, "PS",
#'                id = "record_id",
#'                group = "Cohort")
#' }
#'
#' @export
IQROutliers <- function(df, Variable, id = "ID", group = "SourceFile") {

  if (!Variable %in% names(df)) stop("Variable not found in dataframe.")
  if (!id %in% names(df)) stop("ID column not found in dataframe.")
  if (!group %in% names(df)) stop("Group column not found in dataframe.")

  Q1 <- quantile(df[[Variable]], 0.25, na.rm = TRUE)
  Q3 <- quantile(df[[Variable]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1

  df$outlier <- (df[[Variable]] < (Q1 - 1.5 * IQR)) |
    (df[[Variable]] > (Q3 + 1.5 * IQR))

  p_outlier <- ggplot2::ggplot(df,
                               ggplot2::aes(y = !!rlang::sym(Variable),
                                            x = !!rlang::sym(group))) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(
      width = 0.1,
      size = 2,
      alpha = 0.7,
      ggplot2::aes(color = outlier)
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red")
    ) +
    ggplot2::labs(
      title = paste("IQR Outlier Detection for", Variable),
      y = Variable,
      x = "",
      color = "Outlier"
    ) +
    ggplot2::theme_minimal()

  out_df <- dplyr::select(
    dplyr::filter(df, outlier),
    dplyr::all_of(id),
    dplyr::all_of(Variable),
    dplyr::all_of(group),
    outlier
  )

  list(
    outlierdf = out_df,
    p = p_outlier
  )
}
