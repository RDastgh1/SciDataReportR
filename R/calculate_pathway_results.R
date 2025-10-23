#' Calculate Pathway Results for Metabolite Comparisons
#'
#' @param df Data frame containing metabolite and comparison data
#' @param comparison_var Character string specifying the comparison variable name
#' @param covariates Character vector of covariate names (optional)
#' @param metabolites Character vector of metabolite names to analyze
#' @param comparison_type Character string: "auto", "binary", or "continuous"
#' @param use_point_correlation Logical, if TRUE uses point correlation for binary comparisons
#'
#' @return Data frame with metabolite results including fold change or correlation values
#' @export
#'
#' @examples
#' \dontrun{
#' results <- calculate_pathway_results(
#'   df = my_data,
#'   comparison_var = "poorSleep",
#'   covariates = c("Age", "BMI"),
#'   metabolites = c("Tryptophan", "Kynurenine"),
#'   comparison_type = "binary"
#' )
#' }
calculate_pathway_results <- function(df,
                                      comparison_var,
                                      covariates = NULL,
                                      metabolites,
                                      comparison_type = "auto",
                                      use_point_correlation = FALSE) {

  # Ensure we're using dplyr functions
  library(dplyr)
  library(tidyr)

  # Auto-detect comparison type if needed
  if (comparison_type == "auto") {
    unique_vals <- df[[comparison_var]] %>% na.omit() %>% unique()
    if (length(unique_vals) == 2) {
      comparison_type <- "binary"
    } else {
      comparison_type <- "continuous"
    }
  }

  # Filter out rows with NA in comparison variable
  df_filtered <- df %>%
    dplyr::filter(!is.na(!!sym(comparison_var)))

  # If covariates provided, remove rows with NA in covariates
  if (!is.null(covariates)) {
    # Check which variables exist in the dataframe
    all_vars <- c(metabolites, covariates)
    missing_vars <- all_vars[!all_vars %in% names(df_filtered)]
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in dataframe:", paste(missing_vars, collapse = ", ")))
    }
    df_filtered <- df_filtered %>%
      tidyr::drop_na(all_of(all_vars))
  } else {
    df_filtered <- df_filtered %>%
      tidyr::drop_na(all_of(metabolites))
  }

  # Binary comparison with fold change
  if (comparison_type == "binary" && !use_point_correlation) {

    # Get the two groups
    groups <- unique(df_filtered[[comparison_var]])
    if (length(groups) != 2) {
      stop("Binary comparison requires exactly 2 groups")
    }

    # Calculate group means - using summarize_at for compatibility
    group_means <- df_filtered %>%
      dplyr::group_by(!!sym(comparison_var)) %>%
      dplyr::summarize_at(vars(all_of(metabolites)), ~mean(.x, na.rm = TRUE)) %>%
      tidyr::pivot_longer(-all_of(comparison_var), names_to = "Metabolite", values_to = "Mean")

    # Reshape for fold change
    mean_wide <- group_means %>%
      tidyr::pivot_wider(names_from = all_of(comparison_var), values_from = Mean)

    # Get column names dynamically
    col_names <- names(mean_wide)
    group1_col <- col_names[2]
    group2_col <- col_names[3]

    # Calculate fold change
    mean_wide$Group1 <- mean_wide[[group1_col]]
    mean_wide$Group2 <- mean_wide[[group2_col]]
    mean_wide$fold_change <- mean_wide$Group2 / mean_wide$Group1
    mean_wide$log2_fc <- log2(mean_wide$fold_change)
    mean_wide$`% Change` <- ((mean_wide$Group2 - mean_wide$Group1) / mean_wide$Group1) * 100

    # Remove original group columns
    mean_wide <- mean_wide[, !names(mean_wide) %in% c(group1_col, group2_col)]

    # Run adjusted models and extract p-values
    adjusted_p <- lapply(metabolites, function(metab) {
      # Create safe metabolite name
      metab_safe <- paste0("`", metab, "`")

      if (!is.null(covariates)) {
        formula_str <- paste(metab_safe, "~", comparison_var, "+", paste(covariates, collapse = " + "))
      } else {
        formula_str <- paste(metab_safe, "~", comparison_var)
      }

      tryCatch({
        formula <- as.formula(formula_str)
        model <- lm(formula, data = df_filtered)

        # Get coefficient name
        coef_names <- names(coef(model))
        comp_coef <- coef_names[grep(comparison_var, coef_names)][1]

        if (!is.na(comp_coef) && comp_coef %in% rownames(summary(model)$coefficients)) {
          p_val <- summary(model)$coefficients[comp_coef, "Pr(>|t|)"]
        } else {
          p_val <- NA
        }

        data.frame(Metabolite = metab, p_value = p_val, stringsAsFactors = FALSE)
      }, error = function(e) {
        warning(paste("Error fitting model for", metab, ":", e$message))
        data.frame(Metabolite = metab, p_value = NA, stringsAsFactors = FALSE)
      })
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(p_adj = p.adjust(p_value, method = "fdr"))

    # Combine results
    results_table <- mean_wide %>%
      dplyr::left_join(adjusted_p, by = "Metabolite")

    return(results_table)
  }

  # Binary comparison with point correlation
  if (comparison_type == "binary" && use_point_correlation) {

    # Convert binary to numeric (0, 1)
    df_filtered <- df_filtered %>%
      dplyr::mutate(comparison_numeric = as.numeric(factor(!!sym(comparison_var))) - 1)

    results_list <- lapply(metabolites, function(metab) {

      tryCatch({
        if (!is.null(covariates)) {
          # Partial correlation
          metab_formula <- as.formula(paste0("`", metab, "` ~", paste(covariates, collapse = " + ")))
          comp_formula <- as.formula(paste("comparison_numeric ~", paste(covariates, collapse = " + ")))

          metab_resid <- residuals(lm(metab_formula, data = df_filtered))
          comp_resid <- residuals(lm(comp_formula, data = df_filtered))

          cor_result <- cor.test(metab_resid, comp_resid, method = "pearson")
          correlation <- cor_result$estimate
          p_value <- cor_result$p.value

        } else {
          # Simple point-biserial correlation
          cor_result <- cor.test(df_filtered[[metab]], df_filtered$comparison_numeric, method = "pearson")
          correlation <- cor_result$estimate
          p_value <- cor_result$p.value
        }

        data.frame(
          Metabolite = metab,
          correlation = correlation,
          p_value = p_value,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        warning(paste("Error calculating correlation for", metab, ":", e$message))
        data.frame(
          Metabolite = metab,
          correlation = NA,
          p_value = NA,
          stringsAsFactors = FALSE
        )
      })
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(p_adj = p.adjust(p_value, method = "fdr"))

    return(results_list)
  }

  # Continuous comparison (correlation)
  if (comparison_type == "continuous") {

    results_list <- lapply(metabolites, function(metab) {

      tryCatch({
        if (!is.null(covariates)) {
          # Partial correlation using residuals
          metab_formula <- as.formula(paste0("`", metab, "` ~", paste(covariates, collapse = " + ")))
          comp_formula <- as.formula(paste0("`", comparison_var, "` ~", paste(covariates, collapse = " + ")))

          metab_resid <- residuals(lm(metab_formula, data = df_filtered))
          comp_resid <- residuals(lm(comp_formula, data = df_filtered))

          cor_result <- cor.test(metab_resid, comp_resid, method = "pearson")
          correlation <- cor_result$estimate
          p_value <- cor_result$p.value

        } else {
          # Simple correlation
          cor_result <- cor.test(df_filtered[[metab]], df_filtered[[comparison_var]], method = "pearson")
          correlation <- cor_result$estimate
          p_value <- cor_result$p.value
        }

        data.frame(
          Metabolite = metab,
          correlation = correlation,
          p_value = p_value,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        warning(paste("Error calculating correlation for", metab, ":", e$message))
        data.frame(
          Metabolite = metab,
          correlation = NA,
          p_value = NA,
          stringsAsFactors = FALSE
        )
      })
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(p_adj = p.adjust(p_value, method = "fdr"))

    return(results_list)
  }
}
