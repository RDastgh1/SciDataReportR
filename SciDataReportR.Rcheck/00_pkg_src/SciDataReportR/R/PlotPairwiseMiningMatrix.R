#' Plot a mixed-type pairwise mining matrix
#'
#' Builds referent-centred phenotype contrasts for continuous and categorical
#' measures. It follows [MakePairwiseHeatmap()] conventions while retaining
#' categorical category-level contrasts instead of silently discarding them.
#'
#' @param data A data frame with labelled variables.
#' @param group_var Character scalar naming the grouping variable.
#' @param variables Character vector of continuous or categorical variables.
#' @param Referent Character scalar naming the referent level of `group_var`.
#' @param covariates Optional character vector of covariates.
#' @param adjust_scope Multiple-comparison correction scope: `"per_group"`,
#'   `"per_variable"`, `"matrix"`, or `"none"`.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param star_p Which p-values drive cell stars: `"raw"`, `"adjusted"`, or
#'   `"none"`.
#' @param variable_metadata Optional data frame with `Variable` and optional
#'   `AnchorN`, `AlignedN`, `FilledPrior`, `FilledFuture`, and `TimeExtended`
#'   columns. These fields are joined to the returned audit table.
#' @param max_levels Maximum categorical levels allowed for one variable.
#' @param continuous_fill_limits,categorical_fill_limits Optional symmetric
#'   plotting limits for the continuous and categorical matrices.
#' @param x_axis_text_angle Numeric angle for phenotype labels.
#' @param row_label_width Approximate character width used to wrap displayed
#'   variable and category labels.
#'
#' @return An object of class `"SciDataReportRPairwiseMiningMatrix"` with
#'   `Plots`, `Results`, `OmnibusResults`, `Settings`, `Models`, and `Warnings`.
#'   `Plots$Continuous` is a reference-SD mean-difference matrix;
#'   `Plots$Categorical` is a percentage-point prevalence-difference matrix;
#'   and `Plots$Omnibus` summarizes the overall phenotype association for every
#'   source variable.
#'
#' @details
#' Continuous cells are `Group - Referent` contrasts after scaling each outcome
#' to the referent mean and standard deviation. Categorical variables expand to
#' one row per observed level, with cells equal to `Group - Referent` prevalence
#' in percentage points. The two plot types intentionally have separate scales.
#' Do not compare their colour magnitudes as if they were the same effect size.
#'
#' @export
PlotPairwiseMiningMatrix <- function(data,
                                     group_var,
                                     variables,
                                     Referent,
                                     covariates = NULL,
                                     adjust_scope = c("per_group", "per_variable", "matrix", "none"),
                                     p_adjust_method = c("fdr", "bonferroni", "holm", "none"),
                                     star_p = c("raw", "adjusted", "none"),
                                     variable_metadata = NULL,
                                     max_levels = 30,
                                     continuous_fill_limits = NULL,
                                     categorical_fill_limits = NULL,
                                     x_axis_text_angle = 0,
                                     row_label_width = 58) {
  adjust_scope <- match.arg(adjust_scope)
  p_adjust_method <- match.arg(p_adjust_method)
  star_p <- match.arg(star_p)

  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is.character(group_var) || length(group_var) != 1 || !group_var %in% names(data)) {
    stop("group_var must name one column in data.")
  }
  if (!is.character(variables) || !length(variables) || !all(variables %in% names(data))) {
    stop("variables must be a non-empty character vector of columns in data.")
  }
  if (!is.character(Referent) || length(Referent) != 1) {
    stop("Referent must be supplied as a single character level name.")
  }
  if (!is.null(covariates) && !all(covariates %in% names(data))) {
    stop("Covariate(s) not found: ", paste(setdiff(covariates, names(data)), collapse = ", "))
  }
  if (!is.data.frame(variable_metadata) && !is.null(variable_metadata)) {
    stop("variable_metadata must be NULL or a data frame.")
  }
  if (!is.null(variable_metadata) && !all(c("Variable") %in% names(variable_metadata))) {
    stop("variable_metadata must include a Variable column.")
  }

  variables <- unique(variables)
  df_Data <- data
  df_Data$.PMM_Group <- if (is.factor(df_Data[[group_var]])) droplevels(df_Data[[group_var]]) else factor(df_Data[[group_var]])
  if (!Referent %in% levels(df_Data$.PMM_Group)) stop("Referent level not found: ", Referent)
  group_levels <- setdiff(levels(df_Data$.PMM_Group), Referent)
  if (!length(group_levels)) stop("At least one non-referent group is required.")

  MakeContrast <- function(level) {
    contrast <- stats::setNames(rep(0, length(levels(df_Data$.PMM_Group))), levels(df_Data$.PMM_Group))
    contrast[Referent] <- -1
    contrast[level] <- 1
    out <- list(contrast)
    names(out) <- level
    out
  }
  AddStars <- function(p) {
    dplyr::case_when(is.na(p) ~ "", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "")
  }
  AdjustP <- function(df_Results) {
    if (!nrow(df_Results) || p_adjust_method == "none" || adjust_scope == "none") {
      df_Results$AdjustedPValue <- df_Results$PValue
      return(df_Results)
    }
    df_Results$AdjustedPValue <- NA_real_
    split_key <- switch(adjust_scope,
      matrix = rep("Matrix", nrow(df_Results)),
      per_group = df_Results$Group,
      per_variable = paste(df_Results$Variable, df_Results$Level, sep = "::")
    )
    for (idx in split(seq_len(nrow(df_Results)), split_key)) {
      df_Results$AdjustedPValue[idx] <- stats::p.adjust(df_Results$PValue[idx], method = p_adjust_method)
    }
    df_Results
  }
  MakeMatrixPlot <- function(df_Results, value_col, limits, scale_name) {
    if (!nrow(df_Results)) return(ggplot2::ggplot() + ggplot2::theme_void())
    row_labels <- unique(df_Results$DisplayRowLabel)
    df_Results$DisplayRowLabel <- factor(df_Results$DisplayRowLabel, levels = rev(row_labels))
    df_Results$ColumnLabel <- factor(df_Results$Group, levels = group_levels)
    if (is.null(limits)) {
      max_abs <- max(abs(df_Results[[value_col]]), na.rm = TRUE)
      limits <- c(-max_abs, max_abs)
    }
    plot <- ggplot2::ggplot(df_Results, ggplot2::aes(x = .data$ColumnLabel, y = .data$DisplayRowLabel, fill = .data[[value_col]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::scale_fill_gradient2(low = "#52BCA3FF", mid = "white", high = "#E58606FF", midpoint = 0, limits = limits, oob = scales::squish, name = scale_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, hjust = if (x_axis_text_angle == 0) 0.5 else 1))
    if (star_p != "none") {
      plot <- plot + ggplot2::geom_text(ggplot2::aes(label = .data$SignificanceLabel), color = "black", size = 3)
    }
    plot
  }

  df_ContinuousResults <- tibble::tibble()
  df_CategoricalResults <- tibble::tibble()
  df_OmnibusResults <- tibble::tibble()
  models <- list()
  warnings <- character()

  for (v in variables) {
    label <- .SdrLabelOrName(df_Data, v)
    values <- df_Data[[v]]
    if (is.numeric(values)) {
      ref_values <- values[df_Data$.PMM_Group == Referent]
      ref_mean <- mean(ref_values, na.rm = TRUE)
      ref_sd <- stats::sd(ref_values, na.rm = TRUE)
      if (!is.finite(ref_sd) || ref_sd == 0) {
        warnings <- c(warnings, paste0(v, ": excluded because the referent SD is zero or missing."))
        next
      }
      scaled_name <- paste0(".PMM_", v)
      df_Model <- df_Data %>% dplyr::mutate("{scaled_name}" := (.data[[v]] - ref_mean) / ref_sd)
      formula <- stats::reformulate(c(".PMM_Group", covariates), response = scaled_name)
      fit <- tryCatch(stats::lm(formula, data = df_Model), error = function(e) e)
      if (inherits(fit, "error")) {
        warnings <- c(warnings, paste0(v, ": ", conditionMessage(fit)))
        next
      }
      models[[v]] <- fit
      emm <- emmeans::emmeans(fit, specs = ".PMM_Group")
      contrast_list <- lapply(group_levels, MakeContrast)
      contrast_list <- unlist(contrast_list, recursive = FALSE)
      contrast_df <- as.data.frame(summary(emmeans::contrast(emm, method = contrast_list, adjust = "none"), infer = c(TRUE, TRUE)))
      df_ContinuousResults <- dplyr::bind_rows(df_ContinuousResults, tibble::tibble(
        Variable = v, VariableLabel = label, Level = NA_character_, RowLabel = label,
        Type = "Continuous", Group = contrast_df$contrast, Referent = Referent,
        N = stats::nobs(fit), NGroup = vapply(contrast_df$contrast, function(g) sum(df_Model$.PMM_Group == g & !is.na(df_Model[[scaled_name]])), numeric(1)),
        NReferent = sum(df_Model$.PMM_Group == Referent & !is.na(df_Model[[scaled_name]])),
        Missing = sum(is.na(values)), EstimatedMeanDifference = contrast_df$estimate,
        PercentagePointDifference = NA_real_, PValue = contrast_df$p.value,
        Test = "Linear model + emmeans referent contrast", ModelFormula = paste(deparse(formula), collapse = " ")
      ))
      anova_table <- stats::anova(fit)
      group_row <- which(rownames(anova_table) == ".PMM_Group")
      ss_group <- anova_table$`Sum Sq`[group_row]
      ss_error <- anova_table$`Sum Sq`[nrow(anova_table)]
      df_OmnibusResults <- dplyr::bind_rows(df_OmnibusResults, tibble::tibble(Variable = v, VariableLabel = label, Type = "Continuous", N = stats::nobs(fit), EffectSize = ss_group / (ss_group + ss_error), PValue = anova_table$`Pr(>F)`[group_row], Test = "Linear-model partial eta squared", ModelFormula = paste(deparse(formula), collapse = " ")))
    } else {
      factor_values <- if (is.logical(values)) factor(values, levels = c(FALSE, TRUE)) else droplevels(factor(values))
      levels_values <- levels(factor_values)
      if (length(levels_values) < 2 || length(levels_values) > max_levels) {
        warnings <- c(warnings, paste0(v, ": excluded because it has ", length(levels_values), " observed categorical levels."))
        next
      }
      for (level in levels_values) {
        binary_name <- paste0(".PMM_", v, "_", make.names(level))
        df_Model <- df_Data %>% dplyr::mutate("{binary_name}" := as.integer(factor_values == level))
        formula <- stats::reformulate(c(".PMM_Group", covariates), response = binary_name)
        fit <- tryCatch(stats::glm(formula, data = df_Model, family = stats::binomial()), error = function(e) e)
        if (inherits(fit, "error") || !isTRUE(fit$converged)) {
          warnings <- c(warnings, paste0(v, " / ", level, ": categorical model unavailable."))
          next
        }
        models[[paste(v, level, sep = "::")]] <- fit
        emm <- emmeans::regrid(emmeans::emmeans(fit, specs = ".PMM_Group"), transform = "response")
        contrast_list <- lapply(group_levels, MakeContrast)
        contrast_list <- unlist(contrast_list, recursive = FALSE)
        contrast_df <- as.data.frame(summary(emmeans::contrast(emm, method = contrast_list, adjust = "none"), infer = c(TRUE, TRUE)))
        df_CategoricalResults <- dplyr::bind_rows(df_CategoricalResults, tibble::tibble(
          Variable = v, VariableLabel = label, Level = level, RowLabel = paste0(label, ": ", level),
          Type = "Categorical", Group = contrast_df$contrast, Referent = Referent,
          N = stats::nobs(fit), NGroup = vapply(contrast_df$contrast, function(g) sum(df_Model$.PMM_Group == g & !is.na(factor_values)), numeric(1)),
          NReferent = sum(df_Model$.PMM_Group == Referent & !is.na(factor_values)),
          Missing = sum(is.na(values)), EstimatedMeanDifference = NA_real_,
          PercentagePointDifference = 100 * contrast_df$estimate, PValue = contrast_df$p.value,
          Test = "Binomial model + emmeans prevalence contrast", ModelFormula = paste(deparse(formula), collapse = " ")
        ))
      }
      table_values <- table(df_Data$.PMM_Group, factor_values, useNA = "no")
      chi <- suppressWarnings(tryCatch(stats::chisq.test(table_values), error = function(e) NULL))
      if (!is.null(chi)) {
        n_table <- sum(table_values)
        cramers_v <- sqrt(as.numeric(chi$statistic) / (n_table * min(nrow(table_values) - 1, ncol(table_values) - 1)))
        df_OmnibusResults <- dplyr::bind_rows(df_OmnibusResults, tibble::tibble(Variable = v, VariableLabel = label, Type = "Categorical", N = n_table, EffectSize = cramers_v, PValue = chi$p.value, Test = "Chi-square Cramer's V", ModelFormula = if (is.null(covariates)) "Unadjusted contingency table" else "Unadjusted contingency table; directed contrasts include covariates"))
      }
    }
  }

  results <- dplyr::bind_rows(df_ContinuousResults, df_CategoricalResults) %>% AdjustP()
  results$SignificanceLabel <- AddStars(if (star_p == "adjusted") results$AdjustedPValue else if (star_p == "raw") results$PValue else NA_real_)
  if (!is.null(variable_metadata) && nrow(results)) results <- dplyr::left_join(results, variable_metadata, by = "Variable")
  if (!is.null(variable_metadata) && nrow(df_OmnibusResults)) df_OmnibusResults <- dplyr::left_join(df_OmnibusResults, variable_metadata, by = "Variable")
  if (nrow(results) && all(c("AnchorN", "AlignedN", "FilledPrior", "FilledFuture", "TimeExtended") %in% names(results))) {
    results <- results %>% dplyr::mutate(
      AlignmentLabel = paste0(" [A ", .data$AnchorN, "/", .data$AlignedN, "; +/- ", .data$FilledPrior, "/", .data$FilledFuture, if_else(.data$TimeExtended, " dagger", ""), "]"),
      DisplayRowLabel = paste0(.data$RowLabel, .data$AlignmentLabel)
    )
  } else if (nrow(results)) {
    results$DisplayRowLabel <- results$RowLabel
  }
  if (nrow(results)) results$DisplayRowLabel <- stringr::str_wrap(results$DisplayRowLabel, width = row_label_width)
  df_OmnibusResults <- AdjustP(df_OmnibusResults)
  df_OmnibusResults$SignificanceLabel <- AddStars(if (star_p == "adjusted") df_OmnibusResults$AdjustedPValue else if (star_p == "raw") df_OmnibusResults$PValue else NA_real_)

  continuous_plot <- MakeMatrixPlot(dplyr::filter(results, Type == "Continuous"), "EstimatedMeanDifference", continuous_fill_limits, "Referent SD units")
  categorical_plot <- MakeMatrixPlot(dplyr::filter(results, Type == "Categorical"), "PercentagePointDifference", categorical_fill_limits, "Percentage-point difference")
  df_OmnibusResults$DisplayVariableLabel <- stringr::str_wrap(df_OmnibusResults$VariableLabel, width = row_label_width)
  omnibus_plot <- ggplot2::ggplot(df_OmnibusResults, ggplot2::aes(x = .data$EffectSize, y = reorder(.data$DisplayVariableLabel, .data$EffectSize), color = .data$Type, shape = .data$SignificanceLabel)) +
    ggplot2::geom_point(size = 3) + ggplot2::theme_bw() + ggplot2::theme(axis.title.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) + ggplot2::labs(x = "Overall phenotype-association magnitude", color = "Measure type", shape = "Raw p value")

  out <- list(
    Plots = list(Continuous = continuous_plot, Categorical = categorical_plot, Omnibus = omnibus_plot),
    Results = results,
    OmnibusResults = df_OmnibusResults,
    Settings = list(group_var = group_var, variables = variables, Referent = Referent, covariates = covariates, adjust_scope = adjust_scope, p_adjust_method = p_adjust_method, star_p = star_p, continuous_fill_limits = continuous_fill_limits, categorical_fill_limits = categorical_fill_limits, row_label_width = row_label_width),
    Models = models,
    Warnings = unique(warnings)
  )
  class(out) <- c("SciDataReportRPairwiseMiningMatrix", class(out))
  out
}
