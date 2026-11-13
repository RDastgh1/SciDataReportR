#' Analyze differential feature abundance with moderated statistics
#'
#' Fits one linear model per feature with `limma` and moderates the
#' feature-wise variances by empirical Bayes, which stabilizes inference in
#' small samples. Feature abundance is the response:
#' `feature ~ predictor + covariates`. This is the feature-neutral engine
#' shared by the omics packages (proteins, metabolites, and similar wide
#' feature matrices).
#'
#' @param data A wide data frame with samples in rows, containing the feature
#'   columns, the predictor, and optional covariates. Abundances should be on a
#'   log (or otherwise approximately normal) scale.
#' @param variables Character vector of numeric feature columns.
#' @param predictor_var Character. Group column (character, factor, or logical)
#'   or a numeric column. Groups are compared pairwise; a numeric predictor
#'   gives the change in abundance per unit.
#' @param covariates Optional character vector of adjustment columns.
#' @param contrasts Optional group comparisons, as a list of length-two
#'   character vectors `c(numerator, denominator)`; for example
#'   `list(c("Case", "Control"))` estimates Case minus Control. Defaults to
#'   every other group versus `reference`. Ignored for a numeric predictor.
#' @param reference Reference group for the default contrasts. Defaults to the
#'   first factor level.
#' @param block_var Optional column identifying repeated or paired samples from
#'   the same participant. The within-participant correlation is estimated with
#'   `limma::duplicateCorrelation()` (requires the `statmod` package).
#' @param trend Logical. Allow the prior variance to depend on mean abundance
#'   (`limma::eBayes(trend = TRUE)`).
#' @param robust Logical. Make the variance moderation robust to outlier
#'   features (`limma::eBayes(robust = TRUE)`).
#' @param adjust_method P-value adjustment method, applied across features
#'   within each contrast.
#' @param alpha Significance threshold for `Significant` (P-value) and
#'   `FDRSignificant`.
#' @param min_samples Minimum non-missing samples required to test a feature.
#' @param feature_label Singular noun used in the methods text and plot labels,
#'   such as `"protein"` or `"metabolite"`.
#'
#' @return A `feature_differential_abundance` object with:
#' \describe{
#'   \item{`ResultsTable`}{One row per contrast and feature: `Contrast`,
#'     `Variable`, `Effect`/`Estimate` (difference in mean abundance, or slope),
#'     `EffectType`, `ConfLow`, `ConfHigh`, `StdError`, `TStatistic`, `PValue`,
#'     `FDR`, `AveAbundance`, `N`, `Significant`, `FDRSignificant`, and
#'     `Direction`.}
#'   \item{`SampleAudit`}{Samples supplied, dropped for an incomplete design,
#'     and analyzed.}
#'   \item{`ExcludedFeatures`}{Features not tested, with the reason.}
#'   \item{`Design`, `Fit`}{The design matrix and the `limma` fit.}
#'   \item{`Metadata`}{Settings used, contrast names, group sizes, the
#'     within-block correlation, and a `Methods` sentence.}
#' }
#' @details Samples missing the predictor, a covariate, or the block are
#'   dropped and counted in `SampleAudit`. Missing feature values are handled
#'   per feature.
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   df_wide <- data.frame(
#'     Group = rep(c("Control", "Case"), each = 10),
#'     Age = round(runif(20, 30, 70)),
#'     matrix(rnorm(20 * 30), nrow = 20, dimnames = list(NULL, paste0("F", 1:30)))
#'   )
#'   df_wide$F1 <- df_wide$F1 + 2 * (df_wide$Group == "Case")
#'   da <- AnalyzeFeatureDifferentialAbundance(
#'     df_wide, paste0("F", 1:30), "Group",
#'     covariates = "Age", reference = "Control"
#'   )
#'   da
#'   PlotFeatureDifferentialAbundance(da)
#' }
#' @export
AnalyzeFeatureDifferentialAbundance <- function(
    data,
    variables,
    predictor_var,
    covariates = NULL,
    contrasts = NULL,
    reference = NULL,
    block_var = NULL,
    trend = FALSE,
    robust = FALSE,
    adjust_method = "BH",
    alpha = 0.05,
    min_samples = 3,
    feature_label = "feature"
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Moderated differential abundance requires the optional `limma` package (Bioconductor).", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  is_name <- function(x) is.character(x) && length(x) == 1 && !is.na(x)
  if (!is.character(variables) || length(variables) < 2 || anyNA(variables) || anyDuplicated(variables)) {
    stop("The feature columns must name at least two unique columns.", call. = FALSE)
  }
  if (!is_name(predictor_var)) stop("The predictor must be a single column name.", call. = FALSE)
  if (!is.null(block_var) && !is_name(block_var)) stop("The block variable must be NULL or a single column name.", call. = FALSE)
  if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates))) stop("The covariates must be NULL or a character vector of column names.", call. = FALSE)
  if (!is_name(feature_label)) stop("`feature_label` must be a single string.", call. = FALSE)
  design_vars <- c(predictor_var, covariates, block_var)
  missing_cols <- setdiff(c(variables, design_vars), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in `data`: ", paste(utils::head(missing_cols, 10), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(design_vars) || any(design_vars %in% variables)) {
    stop("The predictor, covariates, block variable, and feature columns must not overlap.", call. = FALSE)
  }
  if (!all(vapply(data[variables], is.numeric, logical(1)))) stop("The feature columns must be numeric.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(min_samples) || length(min_samples) != 1 || is.na(min_samples) || min_samples < 3) {
    stop("The minimum number of samples must be at least 3.", call. = FALSE)
  }

  # Samples with a complete design
  complete <- stats::complete.cases(data[design_vars])
  sample_audit <- tibble::tibble(
    SamplesSupplied = nrow(data),
    SamplesDroppedIncompleteDesign = sum(!complete),
    SamplesAnalyzed = sum(complete)
  )
  df_design <- as.data.frame(data[complete, design_vars, drop = FALSE])
  expression <- t(as.matrix(data[complete, variables, drop = FALSE]))
  rownames(expression) <- variables
  colnames(expression) <- NULL

  predictor <- df_design[[predictor_var]]
  is_group <- !is.numeric(predictor)
  # Syntactic stand-in names keep model.matrix safe with any column names.
  model_data <- data.frame(row.names = seq_len(nrow(df_design)))
  if (is_group) {
    predictor <- droplevels(as.factor(predictor))
    if (nlevels(predictor) < 2) {
      stop("The predictor must have at least two groups among samples with a complete design.", call. = FALSE)
    }
    if (!is.null(reference)) {
      if (!is_name(reference) || !reference %in% levels(predictor)) {
        stop("The reference must be one of the predictor groups: ", paste(levels(predictor), collapse = ", "), call. = FALSE)
      }
      predictor <- stats::relevel(predictor, ref = reference)
    }
    group_levels <- levels(predictor)
    model_data$.Group <- factor(paste0("G", as.integer(predictor)), levels = paste0("G", seq_along(group_levels)))
  } else {
    if (length(unique(predictor)) < 2) stop("The predictor must vary across samples.", call. = FALSE)
    if (!is.null(contrasts) || !is.null(reference)) {
      warning("Contrasts and the reference group are ignored for a numeric predictor.", call. = FALSE)
    }
    model_data$.Predictor <- predictor
  }
  for (i in seq_along(covariates)) {
    covariate <- df_design[[covariates[[i]]]]
    model_data[[paste0(".Cov", i)]] <- if (is.numeric(covariate)) covariate else droplevels(as.factor(covariate))
  }
  covariate_terms <- if (length(covariates) > 0) paste0(".Cov", seq_along(covariates)) else character()
  formula <- stats::as.formula(paste("~", paste(c(if (is_group) "0 + .Group" else ".Predictor", covariate_terms), collapse = " + ")))
  design <- stats::model.matrix(formula, data = model_data)
  if (qr(design)$rank < ncol(design)) {
    stop("The design matrix is not full rank: the predictor and covariates are confounded or a group has no samples.", call. = FALSE)
  }
  if (nrow(design) <= ncol(design)) {
    stop("Too few samples (", nrow(design), ") for a model with ", ncol(design), " coefficients.", call. = FALSE)
  }

  # Contrast matrix
  if (is_group) {
    if (is.null(contrasts)) {
      contrasts <- lapply(group_levels[-1], function(level) c(level, group_levels[[1]]))
    } else {
      if (is.character(contrasts) && length(contrasts) == 2) contrasts <- list(contrasts)
      valid <- is.list(contrasts) && length(contrasts) > 0 &&
        all(vapply(contrasts, function(pair) is.character(pair) && length(pair) == 2 && !anyNA(pair) && pair[[1]] != pair[[2]], logical(1)))
      if (!valid) stop("Contrasts must be a list of `c(numerator, denominator)` group pairs.", call. = FALSE)
      unknown <- setdiff(unlist(contrasts), group_levels)
      if (length(unknown) > 0) {
        stop("Unknown group(s) in the contrasts: ", paste(unknown, collapse = ", "), ". Available: ", paste(group_levels, collapse = ", "), call. = FALSE)
      }
    }
    contrast_names <- vapply(contrasts, function(pair) paste(pair[[1]], "vs", pair[[2]]), character(1))
    if (anyDuplicated(contrast_names)) stop("Contrasts must not repeat a comparison.", call. = FALSE)
    contrast_matrix <- matrix(0, nrow = ncol(design), ncol = length(contrasts), dimnames = list(colnames(design), contrast_names))
    for (j in seq_along(contrasts)) {
      contrast_matrix[paste0(".GroupG", match(contrasts[[j]][[1]], group_levels)), j] <- 1
      contrast_matrix[paste0(".GroupG", match(contrasts[[j]][[2]], group_levels)), j] <- -1
    }
    effect_type <- "Mean difference"
  } else {
    contrast_names <- paste0(predictor_var, " (per unit)")
    contrast_matrix <- matrix(as.numeric(colnames(design) == ".Predictor"), ncol = 1, dimnames = list(colnames(design), contrast_names))
    effect_type <- "Slope"
  }

  # Features with enough data
  n_observed <- rowSums(!is.na(expression))
  variance <- apply(expression, 1, stats::var, na.rm = TRUE)
  exclude_reason <- unname(ifelse(
    n_observed < min_samples, "too_few_samples",
    ifelse(is.na(variance) | variance == 0, "no_variance", NA_character_)
  ))
  excluded <- tibble::tibble(
    Variable = variables[!is.na(exclude_reason)],
    N = unname(n_observed[!is.na(exclude_reason)]),
    ExcludeReason = exclude_reason[!is.na(exclude_reason)]
  )
  expression <- expression[is.na(exclude_reason), , drop = FALSE]
  if (nrow(expression) < 2) {
    stop("Fewer than two ", feature_label, "s have enough non-missing, varying values to test.", call. = FALSE)
  }

  correlation <- NULL
  block <- NULL
  if (!is.null(block_var)) {
    if (!requireNamespace("statmod", quietly = TRUE)) stop("A block variable requires the optional `statmod` package.", call. = FALSE)
    block <- as.factor(df_design[[block_var]])
    if (nlevels(block) == length(block)) {
      stop("The block variable has no repeated samples; omit it for independent samples.", call. = FALSE)
    }
    correlation <- limma::duplicateCorrelation(expression, design, block = block)$consensus.correlation
    if (!is.finite(correlation)) stop("The within-block correlation could not be estimated.", call. = FALSE)
  }
  fit <- if (is.null(block)) {
    limma::lmFit(expression, design)
  } else {
    limma::lmFit(expression, design, block = block, correlation = correlation)
  }
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), trend = trend, robust = robust)

  n_tested <- unname(rowSums(!is.na(expression)))
  df_results <- dplyr::bind_rows(lapply(seq_along(contrast_names), function(j) {
    table <- limma::topTable(fit, coef = j, number = Inf, sort.by = "none", confint = TRUE, adjust.method = adjust_method)
    df_contrast <- tibble::tibble(
      Contrast = contrast_names[[j]],
      Variable = rownames(expression),
      Effect = table$logFC,
      EffectType = effect_type,
      Estimate = table$logFC,
      ConfLow = table$CI.L,
      ConfHigh = table$CI.R,
      StdError = table$logFC / table$t,
      TStatistic = table$t,
      PValue = table$P.Value,
      FDR = table$adj.P.Val,
      AveAbundance = table$AveExpr,
      N = n_tested
    )
    df_contrast$Significant <- !is.na(df_contrast$PValue) & df_contrast$PValue < alpha
    df_contrast$FDRSignificant <- !is.na(df_contrast$FDR) & df_contrast$FDR < alpha
    df_contrast$Direction <- ifelse(df_contrast$Effect > 0, "Higher", ifelse(df_contrast$Effect < 0, "Lower", "None"))
    df_contrast
  }))

  structure(
    list(
      ResultsTable = df_results,
      SampleAudit = sample_audit,
      ExcludedFeatures = excluded,
      Design = design,
      Fit = fit,
      Metadata = list(
        Predictor = predictor_var,
        PredictorType = if (is_group) "group" else "numeric",
        Groups = if (is_group) table(Group = predictor) else NULL,
        Covariates = covariates,
        Block = block_var,
        BlockCorrelation = correlation,
        Contrasts = contrast_names,
        Trend = trend,
        Robust = robust,
        AdjustMethod = adjust_method,
        Alpha = alpha,
        FeatureLabel = feature_label,
        NFeatures = nrow(expression),
        NSamples = ncol(expression),
        Methods = paste0(
          "Differential ", feature_label, " abundance was tested with limma: a linear model was fitted per ", feature_label,
          " (abundance ~ ", paste(c(predictor_var, covariates), collapse = " + "), ")",
          if (!is.null(block_var)) paste0(", with within-", block_var, " correlation estimated by duplicateCorrelation") else "",
          ", and variances were moderated by empirical Bayes",
          if (trend) " with an abundance trend" else "",
          if (robust) " using robust estimation" else "",
          ". P-values were adjusted across ", feature_label, "s within each contrast using ", adjust_method, "."
        )
      )
    ),
    class = c("feature_differential_abundance", "list")
  )
}

#' Print a feature differential abundance object
#'
#' @param x A `feature_differential_abundance` object, as returned by
#'   [AnalyzeFeatureDifferentialAbundance()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.feature_differential_abundance <- function(x, ...) {
  metadata <- x$Metadata
  label <- if (is.null(metadata$FeatureLabel)) "feature" else metadata$FeatureLabel
  excluded <- x$ExcludedFeatures
  cat("<", class(x)[[1]], ">\n", sep = "")
  cat(sprintf("  Predictor: %s (%s)\n", metadata$Predictor, metadata$PredictorType))
  if (length(metadata$Covariates) > 0) cat(sprintf("  Adjusted:  %s\n", paste(metadata$Covariates, collapse = ", ")))
  if (!is.null(metadata$Block)) cat(sprintf("  Block:     %s (correlation %.2f)\n", metadata$Block, metadata$BlockCorrelation))
  cat(sprintf("  Samples:   %d (%d dropped for incomplete design)\n", metadata$NSamples, x$SampleAudit$SamplesDroppedIncompleteDesign))
  cat(sprintf("  Tested:    %d %ss, %d excluded\n", metadata$NFeatures, label, if (is.null(excluded)) 0L else nrow(excluded)))
  for (contrast in metadata$Contrasts) {
    df_results <- x$ResultsTable[x$ResultsTable$Contrast == contrast, , drop = FALSE]
    cat(sprintf(
      "  %s: %d with FDR < %s (%d with P < %s)\n",
      contrast, sum(df_results$FDRSignificant), metadata$Alpha, sum(df_results$Significant), metadata$Alpha
    ))
  }
  invisible(x)
}

#' Plot differential feature abundance as a volcano plot
#'
#' @param x A `feature_differential_abundance` object, as returned by
#'   [AnalyzeFeatureDifferentialAbundance()] or one of its omics wrappers.
#' @param contrast Contrast to plot. Defaults to the first.
#' @param significance_metric Metric used for coloring and the threshold line:
#'   `"FDR"` (default) or `"PValue"`. The y-axis always shows the P-value.
#' @param alpha Significance threshold. Defaults to the value used in the
#'   analysis.
#' @param min_effect Optional minimum absolute effect for a feature to be
#'   colored as significant.
#' @param label_top_n Number of most significant features to label.
#' @param label_by Column of the results used for labels. Defaults to
#'   `"Variable"`.
#'
#' @return A ggplot object.
#' @export
PlotFeatureDifferentialAbundance <- function(
    x,
    contrast = NULL,
    significance_metric = c("FDR", "PValue"),
    alpha = NULL,
    min_effect = 0,
    label_top_n = 10,
    label_by = "Variable"
) {
  if (!inherits(x, "feature_differential_abundance")) {
    stop("`x` must be a differential abundance object.", call. = FALSE)
  }
  significance_metric <- match.arg(significance_metric)
  contrast_names <- x$Metadata$Contrasts
  if (is.null(contrast)) contrast <- contrast_names[[1]]
  if (!is.character(contrast) || length(contrast) != 1 || !contrast %in% contrast_names) {
    stop("The contrast must be one of: ", paste(contrast_names, collapse = ", "), call. = FALSE)
  }
  if (is.null(alpha)) alpha <- x$Metadata$Alpha
  df_plot <- x$ResultsTable[x$ResultsTable$Contrast == contrast, , drop = FALSE]
  if (!is.character(label_by) || length(label_by) != 1 || !label_by %in% names(df_plot)) {
    stop("The label column must name a column of the results.", call. = FALSE)
  }
  if (!is.numeric(min_effect) || length(min_effect) != 1 || is.na(min_effect) || min_effect < 0) {
    stop("The minimum effect must be a non-negative number.", call. = FALSE)
  }
  below <- !is.na(df_plot[[significance_metric]]) & df_plot[[significance_metric]] < alpha
  passes <- below & abs(df_plot$Effect) >= min_effect
  df_plot$Status <- factor(
    ifelse(!passes, "Not significant", ifelse(df_plot$Effect > 0, "Higher", "Lower")),
    levels = c("Lower", "Not significant", "Higher")
  )
  df_plot$NegLog10P <- -log10(df_plot$PValue)
  df_plot$Label <- dplyr::coalesce(as.character(df_plot[[label_by]]), df_plot$Variable)
  .data <- rlang::.data
  plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$Effect, y = .data$NegLog10P, color = .data$Status)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70") +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_color_manual(
      values = c(Lower = "#2166AC", `Not significant` = "grey65", Higher = "#B2182B"),
      drop = FALSE, name = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = if (identical(x$Metadata$PredictorType, "group")) "Difference in mean abundance" else "Change in abundance per unit",
      y = expression(-log[10] ~ italic(P)),
      title = contrast,
      caption = paste0(
        "Colored: ", significance_metric, " < ", alpha,
        if (min_effect > 0) paste0(" and |effect| >= ", min_effect) else "",
        "; limma moderated t-test"
      )
    )
  # Threshold line: the largest P-value still passing the chosen metric.
  if (any(below)) {
    plot <- plot + ggplot2::geom_hline(yintercept = min(df_plot$NegLog10P[below]), linetype = 2, color = "grey40")
  }
  if (min_effect > 0) {
    plot <- plot + ggplot2::geom_vline(xintercept = c(-min_effect, min_effect), linetype = 2, color = "grey40")
  }
  if (label_top_n > 0 && any(passes)) {
    df_labels <- df_plot[passes, , drop = FALSE]
    df_labels <- utils::head(df_labels[order(df_labels$PValue), , drop = FALSE], label_top_n)
    label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = df_labels, ggplot2::aes(label = .data$Label), size = 3, show.legend = FALSE, max.overlaps = Inf)
    } else {
      ggplot2::geom_text(data = df_labels, ggplot2::aes(label = .data$Label), size = 3, vjust = -0.7, show.legend = FALSE, check_overlap = TRUE)
    }
    plot <- plot + label_layer
  }
  plot
}
