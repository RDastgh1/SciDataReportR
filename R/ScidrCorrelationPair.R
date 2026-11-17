# Shared pairwise engine for continuous and point-biserial heatmaps.
ScidrCorrelationPair <- function(df_Pair, method, min_n, eps, yname, xname) {
  # Positional internal names avoid collisions with covariates named x or y.
  names(df_Pair) <- c("x", "y", if (ncol(df_Pair) > 2L) paste0("covariate", seq_len(ncol(df_Pair) - 2L)) else character(0))
  df_Pair <- stats::na.omit(df_Pair)
  n_complete <- nrow(df_Pair)
  result <- list(r = NA_real_, p = NA_real_, n = n_complete)
  HasVariance <- function(x) {
    variance <- suppressWarnings(stats::var(x))
    is.finite(variance) && variance > eps
  }
  if (n_complete < min_n || !HasVariance(df_Pair$x) || !HasVariance(df_Pair$y)) {
    return(result)
  }
  TestCorrelation <- function(x, y) {
    tryCatch(suppressWarnings(ScidrObserveInference(
      stats::cor.test(x, y, method = method)
    )), error = function(e) NULL)
  }
  vars_Covariates <- setdiff(names(df_Pair), c("x", "y"))
  x <- df_Pair$x
  y <- df_Pair$y
  residual_df <- NULL
  if (length(vars_Covariates)) {
    df_Covariates <- df_Pair[, vars_Covariates, drop = FALSE]
    df_Covariates[] <- lapply(df_Covariates, function(z) {
      if (is.character(z) || is.logical(z)) factor(z) else z
    })
    design <- tryCatch(stats::model.matrix(~ ., data = df_Covariates),
      error = function(e) NULL)
    reason <- NULL
    if (is.null(design)) {
      reason <- "Covariate model matrix failed"
    } else {
      design <- design[, -1, drop = FALSE]
      keep <- vapply(seq_len(ncol(design)), function(i) HasVariance(design[, i]), logical(1))
      design <- design[, keep, drop = FALSE]
      if (!ncol(design)) {
        reason <- "All requested covariate columns collapsed"
      } else {
        design <- cbind(1, design)
        fit_x <- stats::lm.fit(design, x)
        fit_y <- stats::lm.fit(design, y)
        residual_df <- n_complete - fit_x$rank - 1L
        if (residual_df <= 0) {
          reason <- "Insufficient observations for covariate adjustment"
          residual_df <- NULL
        } else {
          x <- fit_x$residuals
          y <- fit_y$residuals
          if (!HasVariance(x) || !HasVariance(y)) return(result)
        }
      }
    }
    if (!is.null(reason)) {
      ScidrRecordInferenceFailure(yname, xname,
        paste0(reason, "; unadjusted correlation used."), n_complete, stage = "adjustment")
    }
  }
  if (!is.null(residual_df) && method == "pearson") {
    result$r <- unname(stats::cor(x, y))
    statistic <- result$r * sqrt(residual_df / max(0, 1 - result$r^2))
    result$p <- 2 * stats::pt(-abs(statistic), df = residual_df)
  } else {
    test <- TestCorrelation(x, y)
    if (!is.null(test)) {
      result$r <- unname(test$estimate)
      result$p <- test$p.value
    }
  }
  result
}
