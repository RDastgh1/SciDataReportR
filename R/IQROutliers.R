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
#'   outlier table. If \code{NULL}, no ID column is included in the returned table.
#'   Defaults to \code{NULL}.
#' @param group A string specifying the grouping or batch column to use on the
#'   x-axis of the diagnostic plot. If \code{NULL}, the function will produce a
#'   single combined boxplot across all rows. Defaults to \code{NULL}.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{outlierdf}: a tibble containing only rows flagged as outliers,
#'       including the ID (when requested), variable, group (when requested), and outlier flag.
#'     \item \code{p}: a ggplot2 object showing a boxplot and jittered points
#'       colored by outlier status.
#'   }
#'
#' @details
#' The color aesthetic is mapped only within \code{geom_jitter()}, ensuring the
#' boxplot is drawn once per group rather than once per outlier class. Missing
#' values are ignored when computing quartiles. If \code{id} or \code{group} are
#' provided they must be present in the input data frame.
#'
#' @references
#' Tukey, J. W. (1977). \emph{Exploratory Data Analysis}. Addison-Wesley.
#'
#' @examples
#' \dontrun{
#'   IQROutliers(df_Revalued_Data, "PS",
#'               id = "record_id",
#'               group = "Cohort")
#'
#'   # Without id or group
#'   IQROutliers(df_Revalued_Data, "PS", id = NULL, group = NULL)
#' }
#'
#' @export
IQROutliers <- function(df, Variable, id = NULL, group = NULL) {

  if (!Variable %in% names(df)) stop("Variable not found in dataframe.")
  if (!is.null(id) && !id %in% names(df)) stop("ID column not found in dataframe.")
  if (!is.null(group) && !group %in% names(df)) stop("Group column not found in dataframe.")

  if (!is.numeric(df[[Variable]])) stop("Variable must be numeric.")

  Q1 <- quantile(df[[Variable]], 0.25, na.rm = TRUE)
  Q3 <- quantile(df[[Variable]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1



  ## compute outlier flag robustly using vectorized logicals
  df$.iqr_outlier <- (df[[Variable]] < (Q1 - 1.5 * IQR)) |
    (df[[Variable]] > (Q3 + 1.5 * IQR))

  ## prepare grouping column for plotting: if group is NULL make a single group
  if (is.null(group)) {
    df$.group_for_plot <- factor("All")
    x_aes_sym <- rlang::sym(".group_for_plot")
    include_group_in_output <- FALSE
  } else {
    df$.group_for_plot <- df[[group]]
    x_aes_sym <- rlang::sym(".group_for_plot")
    include_group_in_output <- TRUE
  }

  p_outlier <- ggplot2::ggplot(df,
                               ggplot2::aes(y = !!rlang::sym(Variable),
                                            x = !!x_aes_sym)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(
      width = 0.1,
      size = 2,
      alpha = 0.7,
      ggplot2::aes(color = .iqr_outlier)
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red")
    ) +
    ggplot2::labs(
      title = paste("IQR Outlier Detection for", Variable),
      y = Variable,
      x = if (is.null(group)) "" else group,
      color = "Outlier"
    ) +
    ggplot2::theme_minimal()

  ## build the outlier tibble: include id and group only if requested
  out_select_cols <- c()
  if (!is.null(id)) out_select_cols <- c(out_select_cols, id)
  out_select_cols <- c(out_select_cols, Variable)
  if (include_group_in_output) out_select_cols <- c(out_select_cols, group)
  out_select_cols <- c(out_select_cols, ".iqr_outlier")

  out_df <- dplyr::select(
    dplyr::filter(df, .iqr_outlier),
    dplyr::all_of(out_select_cols)
  )

  ## rename internal flag back to 'outlier' for user-facing output
  out_df <- dplyr::rename(out_df, outlier = .iqr_outlier)

  list(
    outlierdf = out_df,
    p = p_outlier
  )
}
