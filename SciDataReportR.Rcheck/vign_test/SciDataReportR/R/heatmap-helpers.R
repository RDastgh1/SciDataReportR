.FormatPValueStars <- function(p_values, thresholds = ScidrStarTiersThree) {
  ScidrPValueStars(p_values, tiers = thresholds, ns_label = "", na_label = "")
}

.ApplyTidyPAdjustment <- function(data,
                                  p_col,
                                  adjust_scope = c("per_group", "per_variable", "matrix", "none"),
                                  group_col = "Group",
                                  variable_col = "Variable",
                                  method = "fdr") {
  adjust_scope <- match.arg(adjust_scope)

  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!p_col %in% names(data)) {
    stop("p_col was not found in data: ", p_col)
  }

  p_values <- suppressWarnings(as.numeric(data[[p_col]]))

  if (identical(adjust_scope, "none") || identical(method, "none")) {
    return(p_values)
  }

  if (identical(adjust_scope, "matrix")) {
    return(ApplyFDRCorrection(p_values, fdr_scope = "matrix", method = method))
  }

  if (identical(adjust_scope, "per_group")) {
    if (!group_col %in% names(data)) {
      stop("group_col was not found in data: ", group_col)
    }
    return(ApplyFDRCorrection(
      p_values,
      fdr_scope = "per_outcome",
      outcome_ids = data[[group_col]],
      method = method
    ))
  }

  if (!variable_col %in% names(data)) {
    stop("variable_col was not found in data: ", variable_col)
  }

  ApplyFDRCorrection(
    p_values,
    fdr_scope = "per_outcome",
    outcome_ids = data[[variable_col]],
    method = method
  )
}

.BuildHeatmapSignificanceLabels <- function(data,
                                            raw_p_col = "PValue",
                                            adjusted_p_col = "AdjustedPValue",
                                            star_p = c("raw", "adjusted", "none"),
                                            adjusted_significance_threshold = 0.05) {
  star_p <- match.arg(star_p)

  if (!raw_p_col %in% names(data)) {
    stop("raw_p_col was not found in data: ", raw_p_col)
  }
  if (!adjusted_p_col %in% names(data)) {
    stop("adjusted_p_col was not found in data: ", adjusted_p_col)
  }

  star_values <- switch(
    star_p,
    raw = data[[raw_p_col]],
    adjusted = data[[adjusted_p_col]],
    none = rep(NA_real_, nrow(data))
  )

  data$SignificanceLabel <- if (identical(star_p, "none")) {
    rep("", nrow(data))
  } else {
    .FormatPValueStars(star_values)
  }

  adjusted_p <- suppressWarnings(as.numeric(data[[adjusted_p_col]]))
  data$IsAdjustedSignificant <- is.finite(adjusted_p) &
    adjusted_p <= adjusted_significance_threshold

  data
}

.GetHeatmapColorScale <- function(low_color = "#52BCA3FF",
                                  mid_color = "white",
                                  high_color = "#E58606FF",
                                  fill_midpoint = 0,
                                  fill_limits = NULL,
                                  fill_oob = scales::squish,
                                  name = "Estimated mean difference") {
  ggplot2::scale_fill_gradient2(
    low = low_color,
    mid = mid_color,
    high = high_color,
    midpoint = fill_midpoint,
    limits = fill_limits,
    oob = fill_oob,
    na.value = "grey90",
    name = name
  )
}

.AddHeatmapSignificanceLayers <- function(plot,
                                          data,
                                          x_col,
                                          y_col,
                                          label_col = "SignificanceLabel",
                                          outline_col = "IsAdjustedSignificant",
                                          adjusted_outline = TRUE,
                                          adjusted_outline_color = "black",
                                          adjusted_outline_linewidth = 1.0,
                                          star_color = "black",
                                          star_size = 4) {
  if (!label_col %in% names(data)) {
    stop("label_col was not found in data: ", label_col)
  }
  if (!outline_col %in% names(data)) {
    stop("outline_col was not found in data: ", outline_col)
  }

  out <- plot

  if (isTRUE(adjusted_outline)) {
    outline_data <- data[data[[outline_col]] %in% TRUE, , drop = FALSE]
    if (nrow(outline_data) > 0) {
      out <- out +
        ggplot2::geom_tile(
          data = outline_data,
          ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
          inherit.aes = FALSE,
          fill = NA,
          color = adjusted_outline_color,
          linewidth = adjusted_outline_linewidth
        )
    }
  }

  out +
    ggplot2::geom_text(
      data = data,
      ggplot2::aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        label = .data[[label_col]]
      ),
      inherit.aes = FALSE,
      color = star_color,
      size = star_size,
      na.rm = TRUE
    )
}

.BuildHeatmapCaption <- function(score_type,
                                 star_p = c("raw", "adjusted", "none"),
                                 adjusted_outline = TRUE) {
  star_p <- match.arg(star_p)

  score_text <- if (identical(score_type, "MScore")) "M-score" else "Z-score"
  star_text <- switch(
    star_p,
    raw = paste0("Stars show raw p-values (", ScidrStarCaptionText(), ")."),
    adjusted = paste0("Stars show adjusted p-values (", ScidrStarCaptionText(), ")."),
    none = "P-value stars are not shown."
  )
  outline_text <- if (isTRUE(adjusted_outline)) {
    "Thick black outlines mark comparisons significant after multiple-comparison correction."
  } else {
    NULL
  }

  paste(
    c(
      paste0("Fill shows estimated marginal mean difference versus referent (Group - Referent) on the ", score_text, " scale."),
      star_text,
      outline_text
    ),
    collapse = " "
  )
}

.DescribePairwiseAdjustment <- function(adjust_scope, method) {
  method_label <- dplyr::case_when(
    identical(method, "none") ~ "No multiple-comparison correction",
    identical(method, "fdr") ~ "FDR correction",
    identical(method, "bonferroni") ~ "Bonferroni correction",
    identical(method, "holm") ~ "Holm correction",
    TRUE ~ paste0(method, " correction")
  )

  if (identical(method, "none") || identical(adjust_scope, "none")) {
    return("No multiple-comparison correction")
  }

  scope_label <- dplyr::case_when(
    identical(adjust_scope, "per_group") ~ "within each group-vs-referent contrast across variables",
    identical(adjust_scope, "per_variable") ~ "within each variable across group-vs-referent contrasts",
    identical(adjust_scope, "matrix") ~ "across all displayed cells",
    TRUE ~ adjust_scope
  )

  paste(method_label, scope_label)
}

.ComputeSymmetricFillLimits <- function(values, fallback = c(-1, 1)) {
  values <- suppressWarnings(as.numeric(values))
  max_abs <- suppressWarnings(max(abs(values[is.finite(values)]), na.rm = TRUE))

  if (!is.finite(max_abs) || max_abs <= 0) {
    return(fallback)
  }

  c(-max_abs, max_abs)
}

.ProjectMScore <- function(data,
                           variables,
                           parameters,
                           names_prefix = "M_",
                           center = TRUE,
                           scale = TRUE) {
  required_cols <- c("Variable", "N", "Median", "MAD")
  missing_cols <- setdiff(required_cols, names(parameters))
  if (length(missing_cols) > 0) {
    stop(
      "M-score parameter data frame must contain columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are not in data: ", paste(missing_vars, collapse = ", "))
  }

  missing_params <- setdiff(variables, parameters$Variable)
  if (length(missing_params) > 0) {
    stop("The following variables have no M-score parameters: ", paste(missing_params, collapse = ", "))
  }

  param_df <- parameters[match(variables, parameters$Variable), , drop = FALSE]

  transform_fun <- function(x, med, mad_value) {
    out <- x
    if (isTRUE(center)) {
      out <- out - med
    }
    if (isTRUE(scale)) {
      if (is.na(mad_value) || mad_value == 0) {
        out[] <- NA_real_
      } else {
        out <- out / mad_value
      }
    }
    out
  }

  mscore_list <- mapply(
    FUN = transform_fun,
    data[variables],
    param_df$Median,
    param_df$MAD,
    SIMPLIFY = FALSE
  )

  mscore_df <- as.data.frame(mscore_list, stringsAsFactors = FALSE)
  names(mscore_df) <- paste0(names_prefix, variables)

  list(
    MScores = mscore_df,
    DataWithM = cbind(data, mscore_df),
    Parameters = param_df,
    Center = center,
    Scale = scale,
    Constant = if ("Constant" %in% names(parameters)) unique(parameters$Constant)[1] else NA_real_
  )
}

.SdrBacktick <- function(x) {
  paste0("`", gsub("`", "\\\\`", as.character(x)), "`")
}

.SdrFormula <- function(lhs, rhs_terms) {
  rhs_terms <- as.character(rhs_terms)
  rhs <- paste(.SdrBacktick(rhs_terms), collapse = " + ")
  stats::as.formula(paste(.SdrBacktick(lhs), "~", rhs))
}

.SdrLabelOrName <- function(data, variable) {
  lab <- tryCatch(
    sjlabelled::get_label(data[[variable]]),
    error = function(e) NULL
  )

  if (is.null(lab) || length(lab) == 0 || is.na(lab) || !nzchar(lab)) {
    return(variable)
  }

  as.character(lab)
}

.ComputePairwiseReferentContrasts <- function(data,
                                              group_var,
                                              variables,
                                              Referent,
                                              covariates = NULL,
                                              transformed_variables = variables,
                                              score_type = "ZScore",
                                              Parametric = TRUE,
                                              adjust_scope = c("per_group", "per_variable", "matrix", "none"),
                                              p_adjust_method = "fdr",
                                              adjusted_significance_threshold = 0.05,
                                              star_p = c("raw", "adjusted", "none"),
                                              return_models = FALSE) {
  adjust_scope <- match.arg(adjust_scope)
  star_p <- match.arg(star_p)

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for MakePairwiseHeatmap().")
  }
  if (!isTRUE(Parametric) && !requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' is required when Parametric = FALSE.")
  }

  if (!group_var %in% names(data)) {
    stop("group_var was not found in data: ", group_var)
  }
  if (!all(transformed_variables %in% names(data))) {
    stop(
      "Transformed variable(s) not found in data: ",
      paste(setdiff(transformed_variables, names(data)), collapse = ", ")
    )
  }
  if (!is.null(covariates) && !all(covariates %in% names(data))) {
    stop(
      "Covariate(s) not found in data: ",
      paste(setdiff(covariates, names(data)), collapse = ", ")
    )
  }

  if (length(variables) != length(transformed_variables)) {
    stop("variables and transformed_variables must have the same length.")
  }

  df <- data
  if (is.factor(df[[group_var]])) {
    df[[group_var]] <- droplevels(df[[group_var]])
  } else if (is.logical(df[[group_var]])) {
    df[[group_var]] <- factor(df[[group_var]], levels = c(FALSE, TRUE))
  } else {
    df[[group_var]] <- factor(df[[group_var]])
  }

  group_levels <- levels(df[[group_var]])
  if (!Referent %in% group_levels) {
    stop("Referent level not found: ", Referent)
  }

  comparison_groups <- setdiff(group_levels, Referent)
  if (!length(comparison_groups)) {
    stop("No nonreferent group levels are available.")
  }

  variable_map <- tibble::tibble(
    Variable = variables,
    TransformedVariable = transformed_variables
  )

  models <- list()
  warnings <- character(0)

  results <- purrr::map_dfr(seq_len(nrow(variable_map)), function(i) {
    variable <- variable_map$Variable[[i]]
    transformed_variable <- variable_map$TransformedVariable[[i]]
    model_terms <- c(group_var, covariates)
    model_formula <- .SdrFormula(transformed_variable, model_terms)
    model_formula_internal <- paste(deparse(model_formula), collapse = " ")
    score_label <- if (identical(score_type, "MScore")) {
      "Referent-scaled M-score"
    } else {
      "Referent-scaled Z-score"
    }
    model_formula_text <- paste(
      score_label,
      "~",
      paste(model_terms, collapse = " + ")
    )
    test_label <- if (isTRUE(Parametric)) {
      "Linear model + emmeans referent contrast"
    } else {
      "Robust linear model (HC3) + emmeans referent contrast"
    }
    adjustment_label <- .DescribePairwiseAdjustment(adjust_scope, p_adjust_method)

    needed_cols <- c(transformed_variable, group_var, covariates)
    df_cc <- df[stats::complete.cases(df[, needed_cols, drop = FALSE]), , drop = FALSE]
    df_cc[[group_var]] <- droplevels(df_cc[[group_var]])

    dropped_levels <- setdiff(group_levels, levels(df_cc[[group_var]]))
    note <- if (length(dropped_levels)) {
      paste0(
        "Dropped group level(s) after complete-case filtering: ",
        paste(dropped_levels, collapse = ", ")
      )
    } else {
      NA_character_
    }

    # Counted before the tibble is built: inside tibble(), `Referent` would
    # resolve to the column of that name created two lines down, which is the
    # referent recycled to one row per comparison group rather than the scalar.
    n_referent <- sum(df_cc[[group_var]] == Referent, na.rm = TRUE)

    template <- tibble::tibble(
      Variable = variable,
      VariableLabel = .SdrLabelOrName(data, variable),
      TransformedVariable = transformed_variable,
      ScoreType = score_type,
      Group = comparison_groups,
      GroupLabel = comparison_groups,
      Referent = Referent,
      ReferentLabel = Referent,
      NGroup = vapply(
        comparison_groups,
        function(g) sum(df_cc[[group_var]] == g, na.rm = TRUE),
        integer(1)
      ),
      NReferent = n_referent,
      EstimatedMarginalMean = NA_real_,
      ReferentEstimatedMarginalMean = NA_real_,
      EstimatedMeanDifference = NA_real_,
      StandardError = NA_real_,
      ConfidenceIntervalLower = NA_real_,
      ConfidenceIntervalUpper = NA_real_,
      PValue = NA_real_,
      Test = test_label,
      Contrast = "Group - Referent",
      AdjustmentScope = adjust_scope,
      AdjustmentMethod = p_adjust_method,
      Adjustment = adjustment_label,
      ModelFormula = model_formula_text,
      ModelFormulaInternal = model_formula_internal,
      NCompleteModel = nrow(df_cc),
      DroppedGroupLevels = if (length(dropped_levels)) paste(dropped_levels, collapse = ", ") else NA_character_,
      ModelStatus = "Not fit",
      Notes = note
    )

    if (!Referent %in% levels(df_cc[[group_var]]) ||
        nlevels(df_cc[[group_var]]) < 2 ||
        nrow(df_cc) < 3) {
      warnings <<- c(warnings, paste0(variable, ": insufficient data after complete-case filtering."))
      return(template)
    }

    fit <- tryCatch(
      stats::lm(model_formula, data = df_cc),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      warnings <<- c(warnings, paste0(variable, ": model failed - ", fit$message))
      template$ModelStatus <- "Model failed"
      template$Notes <- paste(stats::na.omit(c(template$Notes[[1]], fit$message)), collapse = "; ")
      return(template)
    }

    if (isTRUE(return_models)) {
      models[[variable]] <<- fit
    }

    vcov_arg <- NULL
    if (!isTRUE(Parametric)) {
      vcov_arg <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
    }

    emm <- tryCatch(
      if (is.null(vcov_arg)) {
        emmeans::emmeans(fit, specs = stats::as.formula(paste0("~", .SdrBacktick(group_var))))
      } else {
        emmeans::emmeans(
          fit,
          specs = stats::as.formula(paste0("~", .SdrBacktick(group_var))),
          vcov. = vcov_arg
        )
      },
      error = function(e) e
    )

    if (inherits(emm, "error")) {
      warnings <<- c(warnings, paste0(variable, ": emmeans failed - ", emm$message))
      template$ModelStatus <- "emmeans failed"
      template$Notes <- paste(stats::na.omit(c(template$Notes[[1]], emm$message)), collapse = "; ")
      return(template)
    }

    emm_df <- tryCatch(as.data.frame(emm), error = function(e) NULL)
    if (is.null(emm_df) || !group_var %in% names(emm_df)) {
      warnings <<- c(warnings, paste0(variable, ": emmeans output did not contain group_var."))
      template$ModelStatus <- "emmeans output failed"
      return(template)
    }

    emm_levels <- as.character(emm_df[[group_var]])
    emm_est_col <- if ("emmean" %in% names(emm_df)) "emmean" else "estimate"
    ref_emm <- emm_df[[emm_est_col]][match(Referent, emm_levels)]

    contrast_list <- list()
    for (group in comparison_groups) {
      if (!group %in% emm_levels || !Referent %in% emm_levels) next
      v <- rep(0, length(emm_levels))
      v[match(group, emm_levels)] <- 1
      v[match(Referent, emm_levels)] <- -1
      contrast_list[[group]] <- v
    }

    if (!length(contrast_list)) {
      warnings <<- c(warnings, paste0(variable, ": no valid group-vs-referent contrasts."))
      template$ModelStatus <- "No valid contrasts"
      return(template)
    }

    ctr <- tryCatch(emmeans::contrast(emm, method = contrast_list), error = function(e) e)
    if (inherits(ctr, "error")) {
      warnings <<- c(warnings, paste0(variable, ": contrast failed - ", ctr$message))
      template$ModelStatus <- "contrast failed"
      template$Notes <- paste(stats::na.omit(c(template$Notes[[1]], ctr$message)), collapse = "; ")
      return(template)
    }

    ctr_df <- tryCatch(
      as.data.frame(summary(ctr, infer = c(TRUE, TRUE), adjust = "none")),
      error = function(e) NULL
    )

    if (is.null(ctr_df) || !"contrast" %in% names(ctr_df)) {
      warnings <<- c(warnings, paste0(variable, ": contrast summary failed."))
      template$ModelStatus <- "contrast summary failed"
      return(template)
    }

    lower_col <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ctr_df))[1]
    upper_col <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ctr_df))[1]

    got <- tibble::tibble(
      Group = as.character(ctr_df$contrast),
      EstimatedMeanDifference = as.numeric(ctr_df$estimate),
      StandardError = as.numeric(ctr_df$SE),
      ConfidenceIntervalLower = if (!is.na(lower_col)) as.numeric(ctr_df[[lower_col]]) else NA_real_,
      ConfidenceIntervalUpper = if (!is.na(upper_col)) as.numeric(ctr_df[[upper_col]]) else NA_real_,
      PValue = as.numeric(ctr_df$p.value)
    )

    template %>%
      dplyr::mutate(
        EstimatedMarginalMean = emm_df[[emm_est_col]][match(.data$Group, emm_levels)],
        ReferentEstimatedMarginalMean = ref_emm,
        EstimatedMeanDifference = got$EstimatedMeanDifference[match(.data$Group, got$Group)],
        StandardError = got$StandardError[match(.data$Group, got$Group)],
        ConfidenceIntervalLower = got$ConfidenceIntervalLower[match(.data$Group, got$Group)],
        ConfidenceIntervalUpper = got$ConfidenceIntervalUpper[match(.data$Group, got$Group)],
        PValue = got$PValue[match(.data$Group, got$Group)],
        ModelStatus = "OK"
      )
  })

  if (!nrow(results)) {
    results$AdjustedPValue <- numeric(0)
  } else {
    results$AdjustedPValue <- .ApplyTidyPAdjustment(
      results,
      p_col = "PValue",
      adjust_scope = adjust_scope,
      group_col = "Group",
      variable_col = "Variable",
      method = p_adjust_method
    )
  }

  results <- .BuildHeatmapSignificanceLabels(
    results,
    raw_p_col = "PValue",
    adjusted_p_col = "AdjustedPValue",
    star_p = star_p,
    adjusted_significance_threshold = adjusted_significance_threshold
  )

  list(
    Results = results,
    Models = if (isTRUE(return_models)) models else NULL,
    Warnings = unique(warnings)
  )
}
