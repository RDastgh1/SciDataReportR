#' Transform selected features in a wide data frame
#'
#' @param data A data frame with samples in rows.
#' @param variables Character vector of numeric feature columns.
#' @param method One of `"none"`, `"log1p"`, `"log2"`, `"log10"`, or `"glog"`.
#' @param pseudocount Numeric value added before plain-log transformations.
#' @param glog_lambda Positive stabilization value for generalized log.
#'
#' @return A list with transformed `data`, a feature-level `summary`, and a
#'   before/after skewness `plot`.
#' @export
TransformFeatureMatrix <- function(
    data,
    variables = NULL,
    method = c("none", "log1p", "log2", "log10", "glog"),
    pseudocount = 1,
    glog_lambda = 1
) {
  variables <- ValidateFeatureMatrix(data, variables)
  method <- match.arg(method)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 || is.na(pseudocount) || pseudocount <= 0) {
    stop("`pseudocount` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(glog_lambda) || length(glog_lambda) != 1 || is.na(glog_lambda) || glog_lambda <= 0) {
    stop("`glog_lambda` must be a single positive numeric value.", call. = FALSE)
  }

  feature_mat <- as.matrix(data[variables])
  storage.mode(feature_mat) <- "double"
  before_skewness <- vapply(seq_len(ncol(feature_mat)), function(i) ColumnSkewness(feature_mat[, i]), numeric(1))
  transformed <- switch(
    method,
    none = feature_mat,
    log1p = log1p(feature_mat),
    log2 = log2(feature_mat + pseudocount),
    log10 = log10(feature_mat + pseudocount),
    glog = log2((feature_mat + sqrt(feature_mat^2 + glog_lambda^2)) / 2)
  )
  after_skewness <- vapply(seq_len(ncol(transformed)), function(i) ColumnSkewness(transformed[, i]), numeric(1))
  out_data <- data
  out_data[variables] <- as.data.frame(transformed)

  plot_data <- tibble::tibble(
    Skewness = c(before_skewness, after_skewness),
    Stage = rep(c("Before", "After"), each = length(variables))
  )
  plot_data <- dplyr::filter(plot_data, is.finite(Skewness))
  diagnostic_plot <- if (nrow(plot_data) > 0) {
    ggplot2::ggplot(plot_data, ggplot2::aes(x = Stage, y = Skewness, fill = Stage)) +
      ggplot2::geom_boxplot(na.rm = TRUE) +
      ggplot2::theme_minimal() +
      ggplot2::guides(fill = "none")
  } else {
    NULL
  }

  list(
    data = out_data,
    summary = tibble::tibble(Variable = variables, BeforeSkewness = before_skewness, AfterSkewness = after_skewness, Method = method),
    plot = diagnostic_plot
  )
}

ColumnSkewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || stats::sd(x) == 0) return(NA_real_)
  mean((x - mean(x))^3) / stats::sd(x)^3
}
