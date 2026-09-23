#' Plot non-directional mixed-type pairwise mining matrices
#'
#' Builds referent-centred phenotype screening matrices that show the magnitude
#' of pairwise separation only. Signed continuous contrasts and category-level
#' prevalence contrasts are retained in returned audit objects for follow-up,
#' but are deliberately not displayed in the mining matrices.
#'
#' @param data A data frame with labelled variables.
#' @param group_var Character scalar naming the grouping variable.
#' @param variables Character vector of continuous or categorical variables.
#' @param Referent Character scalar naming the referent level of `group_var`.
#' @param covariates Optional character vector of covariates. Continuous
#'   contrasts use these covariates; Cramer's V screening remains unadjusted.
#' @param adjust_scope Multiple-comparison correction scope: `"per_group"`,
#'   `"per_variable"`, `"matrix"`, or `"none"`.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param star_p Which p-values drive cell stars: `"raw"`, `"adjusted"`, or
#'   `"none"`.
#' @param adjusted_outline Logical; outline cells significant after adjustment.
#' @param adjusted_significance_threshold Threshold for adjusted-significant
#'   outlines.
#' @param adjusted_outline_color,adjusted_outline_linewidth Appearance of the
#'   adjusted-significant outline.
#' @param variable_metadata Optional data frame with `Variable` and optional
#'   `AnchorN`, `AlignedN`, `FilledPrior`, `FilledFuture`, and `TimeExtended`.
#' @param max_levels Maximum categorical levels allowed for one variable.
#' @param continuous_fill_limits,categorical_fill_limits Optional non-negative
#'   plotting limits for the continuous and categorical magnitude matrices.
#' @param x_axis_text_angle Numeric angle for phenotype labels.
#' @param row_label_width Approximate character width used to wrap row labels.
#'
#' @return An object of class `"SciDataReportRPairwiseMiningMatrix"` with
#'   `Plots`, `Results`, `ContinuousAudit`, `CategoryLevelResults`, `Settings`,
#'   `Models`, and `Warnings`. `Plots$Continuous` displays absolute
#'   referent-SD mean differences and `Plots$Categorical` displays pairwise
#'   Cramer's V. Both use an independent pale-gray-to-navy magnitude scale.
#'
#' @details
#' The required Referent determines the phenotype comparison columns.
#' Continuous cells are absolute `Group - Referent` contrasts in reference-SD
#' units. Categorical cells are pairwise Cramer's V values from contingency
#' tables, one row per variable. These metrics intentionally have separate
#' legends and should not be compared as interchangeable effect sizes.
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
                                     adjusted_outline = TRUE,
                                     adjusted_significance_threshold = 0.05,
                                     adjusted_outline_color = "black",
                                     adjusted_outline_linewidth = 1.0,
                                     variable_metadata = NULL,
                                     max_levels = 30,
                                     continuous_fill_limits = NULL,
                                     categorical_fill_limits = NULL,
                                     x_axis_text_angle = 0,
                                     row_label_width = 58) {
  adjust_scope <- match.arg(adjust_scope)
  p_adjust_method <- match.arg(p_adjust_method)
  star_p <- match.arg(star_p)
  if (!is.logical(adjusted_outline) || length(adjusted_outline) != 1 || is.na(adjusted_outline)) stop("adjusted_outline must be TRUE or FALSE.")
  if (!is.numeric(adjusted_significance_threshold) || length(adjusted_significance_threshold) != 1 || is.na(adjusted_significance_threshold) || adjusted_significance_threshold <= 0 || adjusted_significance_threshold >= 1) stop("adjusted_significance_threshold must be a single number between 0 and 1.")
  if (!is.numeric(adjusted_outline_linewidth) || length(adjusted_outline_linewidth) != 1 || is.na(adjusted_outline_linewidth) || adjusted_outline_linewidth <= 0) stop("adjusted_outline_linewidth must be a positive number.")
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!is.character(group_var) || length(group_var) != 1 || !group_var %in% names(data)) stop("group_var must name one column in data.")
  if (!is.character(variables) || !length(variables) || !all(variables %in% names(data))) stop("variables must be a non-empty character vector of columns in data.")
  if (!is.character(Referent) || length(Referent) != 1) stop("Referent must be supplied as a single character level name.")
  if (!is.null(covariates) && !all(covariates %in% names(data))) stop("Covariate(s) not found: ", paste(setdiff(covariates, names(data)), collapse = ", "))
  if (!is.null(variable_metadata) && (!is.data.frame(variable_metadata) || !"Variable" %in% names(variable_metadata))) stop("variable_metadata must be NULL or a data frame containing Variable.")

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
    stats::setNames(list(contrast), level)
  }
  AddStars <- function(p) dplyr::case_when(is.na(p) ~ "", p <= 0.001 ~ "***", p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "")
  AdjustP <- function(df_Results) {
    if (!nrow(df_Results)) {
      df_Results$AdjustedPValue <- numeric(0)
      return(df_Results)
    }
    if (p_adjust_method == "none" || adjust_scope == "none") {
      df_Results$AdjustedPValue <- df_Results$PValue
      return(df_Results)
    }
    df_Results$AdjustedPValue <- NA_real_
    split_key <- switch(adjust_scope, matrix = rep("Matrix", nrow(df_Results)), per_group = df_Results$Group, per_variable = df_Results$Variable)
    for (idx in split(seq_len(nrow(df_Results)), split_key)) df_Results$AdjustedPValue[idx] <- stats::p.adjust(df_Results$PValue[idx], method = p_adjust_method)
    df_Results
  }
  AddMetadata <- function(df_Results) {
    if (!nrow(df_Results)) return(df_Results)
    if (!is.null(variable_metadata)) df_Results <- dplyr::left_join(df_Results, variable_metadata, by = "Variable")
    if (!"TimeExtended" %in% names(df_Results)) df_Results$TimeExtended <- FALSE
    df_Results <- df_Results %>% dplyr::mutate(TimeExtended = dplyr::coalesce(.data$TimeExtended, FALSE), DisplayRowLabel = paste0(.data$RowLabel, if_else(.data$TimeExtended, "\u2020", "")))
    df_Results$DisplayRowLabel <- stringr::str_wrap(df_Results$DisplayRowLabel, width = row_label_width)
    df_Results
  }
  MakeMatrixPlot <- function(df_Results, limits, scale_name) {
    if (!nrow(df_Results)) return(ggplot2::ggplot() + ggplot2::theme_void())
    row_labels <- unique(df_Results$DisplayRowLabel)
    df_Results$DisplayRowLabel <- factor(df_Results$DisplayRowLabel, levels = rev(row_labels))
    df_Results$ColumnLabel <- factor(df_Results$Group, levels = group_levels)
    if (is.null(limits)) {
      max_value <- max(df_Results$EffectMagnitude, na.rm = TRUE)
      limits <- c(0, if (is.finite(max_value) && max_value > 0) max_value else 1)
    }
    if (length(limits) != 2 || limits[1] < 0 || limits[2] <= limits[1]) stop("Magnitude plotting limits must be two increasing non-negative values.")
    plot <- ggplot2::ggplot(df_Results, ggplot2::aes(x = .data$ColumnLabel, y = .data$DisplayRowLabel, fill = .data$EffectMagnitude)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::scale_fill_gradient(low = "#F1F5F9", high = "#2166AC", limits = limits, oob = scales::squish, name = scale_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, hjust = if (x_axis_text_angle == 0) 0.5 else 1))
    .AddHeatmapSignificanceLayers(
      plot = plot,
      data = df_Results,
      x_col = "ColumnLabel",
      y_col = "DisplayRowLabel",
      label_col = "SignificanceLabel",
      outline_col = "IsAdjustedSignificant",
      adjusted_outline = adjusted_outline,
      adjusted_outline_color = adjusted_outline_color,
      adjusted_outline_linewidth = adjusted_outline_linewidth,
      star_color = "black",
      star_size = 3
    )
  }

  df_ContinuousResults <- tibble::tibble()
  df_CategoricalResults <- tibble::tibble()
  df_CategoryLevelResults <- tibble::tibble()
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
      contrast_df <- as.data.frame(summary(emmeans::contrast(emm, method = unlist(lapply(group_levels, MakeContrast), recursive = FALSE), adjust = "none"), infer = c(TRUE, TRUE)))
      df_ContinuousResults <- dplyr::bind_rows(df_ContinuousResults, tibble::tibble(
        Variable = v, VariableLabel = label, Level = NA_character_, RowLabel = label, Type = "Continuous", Group = contrast_df$contrast, Referent = Referent,
        N = stats::nobs(fit), NGroup = vapply(contrast_df$contrast, function(g) sum(df_Model$.PMM_Group == g & !is.na(df_Model[[scaled_name]])), numeric(1)), NReferent = sum(df_Model$.PMM_Group == Referent & !is.na(df_Model[[scaled_name]])), Missing = sum(is.na(values)),
        SignedEstimatedMeanDifference = contrast_df$estimate, EstimatedMeanDifference = contrast_df$estimate, PercentagePointDifference = NA_real_, CramersV = NA_real_, EffectMagnitude = abs(contrast_df$estimate), PValue = contrast_df$p.value,
        Test = "Linear model + emmeans referent contrast", ModelFormula = paste(deparse(formula), collapse = " ")
      ))
    } else {
      factor_values <- if (is.logical(values)) factor(values, levels = c(FALSE, TRUE)) else droplevels(factor(values))
      levels_values <- levels(factor_values)
      if (length(levels_values) < 2 || length(levels_values) > max_levels) {
        warnings <- c(warnings, paste0(v, ": excluded because it has ", length(levels_values), " observed categorical levels."))
        next
      }
      if (!is.null(covariates)) warnings <- c(warnings, paste0(v, ": Cramer's V screening is unadjusted; covariates are not applicable to contingency-table magnitude."))
      for (g in group_levels) {
        include <- df_Data$.PMM_Group %in% c(Referent, g) & !is.na(factor_values)
        pair_groups <- droplevels(factor(df_Data$.PMM_Group[include], levels = c(Referent, g)))
        pair_values <- droplevels(factor(factor_values[include]))
        table_values <- table(pair_groups, pair_values)
        chi <- suppressWarnings(tryCatch(stats::chisq.test(table_values), error = function(e) NULL))
        denominator <- sum(table_values) * min(nrow(table_values) - 1, ncol(table_values) - 1)
        if (is.null(chi) || denominator <= 0) {
          warnings <- c(warnings, paste0(v, " / ", g, ": pairwise Cramer's V unavailable."))
          next
        }
        cramers_v <- sqrt(as.numeric(chi$statistic) / denominator)
        df_CategoricalResults <- dplyr::bind_rows(df_CategoricalResults, tibble::tibble(
          Variable = v, VariableLabel = label, Level = NA_character_, RowLabel = label, Type = "Categorical", Group = g, Referent = Referent, N = sum(table_values), NGroup = sum(pair_groups == g), NReferent = sum(pair_groups == Referent), Missing = sum(is.na(values)),
          SignedEstimatedMeanDifference = NA_real_, EstimatedMeanDifference = NA_real_, PercentagePointDifference = NA_real_, CramersV = cramers_v, EffectMagnitude = cramers_v, PValue = chi$p.value,
          Test = "Pairwise chi-square Cramer's V", ModelFormula = "Unadjusted Referent-versus-group contingency table"
        ))
        for (level in levels_values) {
          prevalence <- tapply(pair_values == level, pair_groups, mean)
          df_CategoryLevelResults <- dplyr::bind_rows(df_CategoryLevelResults, tibble::tibble(
            Variable = v, VariableLabel = label, Level = level, RowLabel = paste0(label, ": ", level), Type = "Categorical level", Group = g, Referent = Referent, N = length(pair_values), NGroup = sum(pair_groups == g), NReferent = sum(pair_groups == Referent), Missing = sum(is.na(values)),
            PercentagePointDifference = 100 * (prevalence[[g]] - prevalence[[Referent]]), PValue = NA_real_, Test = "Observed pairwise prevalence difference", ModelFormula = "Unadjusted Referent-versus-group prevalence"
          ))
        }
      }
    }
  }

  results <- dplyr::bind_rows(df_ContinuousResults, df_CategoricalResults) %>% AdjustP() %>% AddMetadata()
  if (nrow(results)) results$SignificanceLabel <- AddStars(if (star_p == "adjusted") results$AdjustedPValue else if (star_p == "raw") results$PValue else NA_real_)
  if (nrow(results)) results$IsAdjustedSignificant <- is.finite(results$AdjustedPValue) & results$AdjustedPValue <= adjusted_significance_threshold
  df_CategoryLevelResults <- df_CategoryLevelResults %>% AdjustP() %>% AddMetadata()
  if (nrow(df_CategoryLevelResults)) df_CategoryLevelResults$SignificanceLabel <- AddStars(if (star_p == "adjusted") df_CategoryLevelResults$AdjustedPValue else if (star_p == "raw") df_CategoryLevelResults$PValue else NA_real_)
  continuous_results <- dplyr::filter(results, .data$Type == "Continuous")
  categorical_results <- dplyr::filter(results, .data$Type == "Categorical")

  out <- list(
    Plots = list(Continuous = MakeMatrixPlot(continuous_results, continuous_fill_limits, "Absolute Referent-scaled difference"), Categorical = MakeMatrixPlot(categorical_results, categorical_fill_limits, "Pairwise Cramer's V")),
    Results = results, ContinuousAudit = continuous_results, CategoryLevelResults = df_CategoryLevelResults,
    Settings = list(group_var = group_var, variables = variables, Referent = Referent, covariates = covariates, adjust_scope = adjust_scope, p_adjust_method = p_adjust_method, star_p = star_p, adjusted_outline = adjusted_outline, adjusted_significance_threshold = adjusted_significance_threshold, adjusted_outline_color = adjusted_outline_color, adjusted_outline_linewidth = adjusted_outline_linewidth, continuous_fill_limits = continuous_fill_limits, categorical_fill_limits = categorical_fill_limits, colors = c(low = "#F1F5F9", high = "#2166AC"), row_label_width = row_label_width),
    Models = models, Warnings = unique(warnings)
  )
  class(out) <- c("SciDataReportRPairwiseMiningMatrix", class(out))
  out
}
