#' Make a Group Comparison Table (gtsummary wrapper)
#'
#' @description
#' This function creates a summary comparison table from a given data frame. It summarizes continuous variables with their mean (and standard deviation) and categorical variables with counts and percentages. Comparisons between groups are performed based on a specified grouping (comparison) variable.
#'
#' @param DataFrame A data frame containing the data for generating the comparison table.
#' @param Variables An optional character vector of variable names (besides the grouping variable) to include in the comparison table. If \code{NULL}, all variables in \code{DataFrame} are used.
#' @param CompVariable A character string specifying the name of the variable used to divide the data into comparison groups.
#' @param ValueDigits An integer specifying the number of digits to display for continuous variable statistics (mean and standard deviation). Default is 2.
#' @param pDigits An integer specifying the number of digits to display for p-values. Default is 3.
#' @param AddEffectSize A logical value indicating whether to add effect size measures. Default is FALSE.
#' @param EffectSizeDigits An integer specifying the number of digits to display for effect size statistics. Default is 2.
#' @param AddPairwise A logical value indicating whether to add pairwise post-hoc comparisons for multi-group analyses. Default is FALSE.
#' @param PairwiseMethod A character string specifying the p-value adjustment method for pairwise comparisons. Default is "bonferroni". Other options include "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none".
#' @param ReferentGroup An optional character string specifying a reference group for pairwise comparisons. If provided, only comparisons between this group and other groups will be shown. Default is NULL (all pairwise comparisons).
#'
#' @return A gtsummary table object with the comparison results.
#'
#' @details
#' The function first determines the set of variables to analyze. If \code{Variables} is \code{NULL}, all variables from \code{DataFrame} are used; otherwise, it ensures that \code{CompVariable} is included once in the analysis. Next, it filters the data frame to include only the selected variables and excludes any factor variables with more than 15 unique levels. The summary is then created using \code{\link[gtsummary]{tbl_summary}} with continuous variables summarized as "mean (sd)" and categorical variables summarized as "n (p)". Additional information, such as sample sizes, p-values, and optionally effect sizes, is appended to the table.
#'
#' When effect sizes are requested, the function calculates Cohen's d for continuous variables and Cramer's V for categorical variables.
#'
#' For multi-group comparisons, post-hoc pairwise tests can be added when \code{AddPairwise = TRUE}. Significant p-values are displayed in bold. If \code{ReferentGroup} is specified, only comparisons between the reference group and other groups will be shown.
#'
#' @note This wrapper function is adapted from code written by Aparna Bhattacharyya.
#' @export
MakeComparisonTable <- function(DataFrame, Variables = NULL, CompVariable, Covariates = NULL,
                                ValueDigits = 2, pDigits = 3,
                                AddEffectSize = FALSE, EffectSizeDigits = 2,
                                AddPairwise = FALSE, PairwiseMethod = "bonferroni",
                                ReferentGroup = NULL) {

  # Check for required packages
  required_packages <- c("car")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(paste("Package", pkg, "is needed for covariate adjustment. Installing it..."))
      install.packages(pkg)
      if (!requireNamespace(pkg, quietly = TRUE)) {
        warning(paste("Failed to install", pkg, "package. Proceeding without covariate adjustment."))
        Covariates <- NULL
      }
    }
  }

  # Load required packages for effect size calculation
  if (AddEffectSize) {
    if (!requireNamespace("effectsize", quietly = TRUE)) {
      warning("Package 'effectsize' is needed for effect size calculations. Installing it...")
      install.packages("effectsize")
      if (!requireNamespace("effectsize", quietly = TRUE)) {
        warning("Failed to install 'effectsize' package. Effect sizes will not be calculated.")
        AddEffectSize <- FALSE
      }
    }
  }

  # Check if covariates are valid
  if (!is.null(Covariates)) {
    missing_covariates <- Covariates[!Covariates %in% names(DataFrame)]
    if (length(missing_covariates) > 0) {
      warning(paste("The following covariates were not found in the data frame and will be ignored:",
                    paste(missing_covariates, collapse = ", ")))
      Covariates <- Covariates[Covariates %in% names(DataFrame)]
    }

    # Remove CompVariable from covariates if present
    Covariates <- Covariates[Covariates != CompVariable]

    # Check if any covariates remain
    if (length(Covariates) == 0) {
      warning("No valid covariates provided. Proceeding without covariate adjustment.")
      Covariates <- NULL
    }
  }

  # If variables of interest are not provided then use all variables from provided data frame
  if (is.null(Variables)) {
    Variables <- colnames(DataFrame)
  } else {
    # Ensure the comparison variable is included once
    Variables <- unique(c(CompVariable, Variables))
  }

  # Ensure covariates are included in the data
  if (!is.null(Covariates)) {
    Variables <- unique(c(Variables, Covariates))
  }

  # Filter the data frame to only include the selected variables.
  # Exclude factor variables that have more than 15 unique levels.
  Data_filtered <- DataFrame %>%
    dplyr::select(dplyr::all_of(Variables)) %>%
    dplyr::select_if(function(col) !(is.factor(col) && nlevels(col) > 15))

  # Convert logical variables to factors for proper handling
  for (var in names(Data_filtered)) {
    if (is.logical(Data_filtered[[var]])) {
      Data_filtered[[var]] <- factor(Data_filtered[[var]], levels = c(FALSE, TRUE), labels = c("No", "Yes"))
    }
  }

  # Create a clean version of the data frame for analysis
  Data_clean <- Data_filtered
  names(Data_clean) <- make.names(names(Data_filtered))
  clean_CompVariable <- make.names(CompVariable)
  clean_Covariates <- make.names(Covariates)

  # Create mapping between original and clean names
  name_mapping <- data.frame(
    original = names(Data_filtered),
    clean = names(Data_clean),
    stringsAsFactors = FALSE
  )

  # Create the comparison table using gtsummary functions and provided arguments.
  CompTable <- gtsummary::tbl_summary(
    data = Data_filtered,
    by = CompVariable,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(
      all_continuous() ~ ValueDigits
    ),
    type = list(
      where(is.double) ~ "continuous",
      where(~ is.factor(.) & nlevels(.) == 2) ~ "categorical"
    )
  ) %>%
    gtsummary::add_n()

  # First add standard p-values without adjustment
  CompTable <- CompTable %>%
    gtsummary::add_p(
      test = list(
        all_continuous() ~ "oneway.test",
        all_categorical() ~ "fisher.test"
      ),
      test.args = list(
        "oneway.test" = list(var.equal = TRUE)  # This makes oneway.test equivalent to aov
      )
    ) %>%
    gtsummary::bold_p() %>%
    gtsummary::modify_fmt_fun(
      p.value ~ function(x) gtsummary::style_pvalue(x, digits = pDigits)
    )

  # If covariates are specified, manually calculate adjusted p-values
  if (!is.null(Covariates)) {
    # Get the variables in the table
    table_vars <- unique(CompTable$table_body$variable)

    # Create a data frame to store adjusted p-values
    adjusted_p_values <- data.frame(
      variable = character(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )

    # Calculate adjusted p-values for each variable
    for (var in table_vars) {
      # Skip the comparison variable itself and covariates
      if (var == CompVariable || var %in% Covariates) {
        next
      }

      # Get the clean variable name
      clean_var <- name_mapping$clean[name_mapping$original == var]

      # Check if variable exists in clean data
      if (!(clean_var %in% names(Data_clean))) {
        next
      }

      # Check variable type
      if (is.numeric(Data_clean[[clean_var]])) {
        # For continuous variables, use linear regression
        tryCatch({
          # Create formula with clean variable names
          formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                paste(clean_Covariates, collapse = " + "))

          # Fit the model
          model <- stats::lm(stats::as.formula(formula_str), data = Data_clean)

          # Get ANOVA results
          anova_result <- car::Anova(model, type = 2)

          # Extract p-value for the grouping variable
          p_value <- anova_result[clean_CompVariable, "Pr(>F)"]

          # Add to adjusted p-values data frame
          adjusted_p_values <- rbind(adjusted_p_values, data.frame(
            variable = var,
            p_value = p_value,
            stringsAsFactors = FALSE
          ))
        }, error = function(e) {
          warning("Error calculating adjusted p-value for ", var, ": ", e$message)
        })
      } else if (is.factor(Data_clean[[clean_var]])) {
        # For categorical variables
        # Check if it's binary or multi-level
        n_levels <- length(levels(Data_clean[[clean_var]]))

        if (n_levels <= 2) {
          # For binary variables, use logistic regression
          tryCatch({
            # Create formula with clean variable names
            formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                  paste(clean_Covariates, collapse = " + "))

            # Fit the model
            model <- stats::glm(stats::as.formula(formula_str), family = binomial(), data = Data_clean)

            # Get Wald test results
            wald_test <- car::Anova(model, type = 2)

            # Extract p-value for the grouping variable
            p_value <- wald_test[clean_CompVariable, "Pr(>Chisq)"]

            # Add to adjusted p-values data frame
            adjusted_p_values <- rbind(adjusted_p_values, data.frame(
              variable = var,
              p_value = p_value,
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            warning("Error calculating adjusted p-value for ", var, ": ", e$message)
          })
        } else {
          # For multi-level categorical variables, use multinomial logistic regression
          if (requireNamespace("nnet", quietly = TRUE)) {
            tryCatch({
              # Create formula with clean variable names
              formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                    paste(clean_Covariates, collapse = " + "))

              # Fit the model
              model <- nnet::multinom(stats::as.formula(formula_str), data = Data_clean, trace = FALSE)

              # Fit null model (without the grouping variable)
              null_formula_str <- paste0(clean_var, " ~ ", paste(clean_Covariates, collapse = " + "))
              null_model <- nnet::multinom(stats::as.formula(null_formula_str), data = Data_clean, trace = FALSE)

              # Likelihood ratio test
              lr_test <- stats::anova(null_model, model)

              # Extract p-value
              p_value <- lr_test$`P(>|Chi|)`[2]

              # Add to adjusted p-values data frame
              adjusted_p_values <- rbind(adjusted_p_values, data.frame(
                variable = var,
                p_value = p_value,
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              warning("Error calculating adjusted p-value for ", var, ": ", e$message)
            })
          } else {
            warning("Package 'nnet' is needed for multinomial regression. Skipping adjustment for ", var)
          }
        }
      }
    }

    # Replace p-values in the table with adjusted p-values
    if (nrow(adjusted_p_values) > 0) {
      # Create a copy of the table body
      new_table_body <- CompTable$table_body

      # Update p-values
      for (i in 1:nrow(adjusted_p_values)) {
        var <- adjusted_p_values$variable[i]
        p_val <- adjusted_p_values$p_value[i]

        # Find rows for this variable
        var_rows <- which(new_table_body$variable == var)

        if (length(var_rows) > 0) {
          # Update p-value
          new_table_body$p.value[var_rows] <- p_val

          # Update p-value styling (bold if significant)
          if (!is.na(p_val) && p_val < 0.05) {
            new_table_body$p.value_fmt[var_rows] <- paste0("**", gtsummary::style_pvalue(p_val, digits = pDigits), "**")
          } else {
            new_table_body$p.value_fmt[var_rows] <- gtsummary::style_pvalue(p_val, digits = pDigits)
          }
        }
      }

      # Replace the table body
      CompTable$table_body <- new_table_body
    }

    # Add footnote about covariate adjustment
    CompTable <- CompTable %>%
      gtsummary::modify_footnote(
        p.value = paste0("P-values adjusted for: ", paste(Covariates, collapse = ", "))
      )

    # Set caption explaining the statistical tests used
    test_caption <- paste0("Statistical tests: Linear models with covariate adjustment (for continuous variables) and ",
                           "logistic/multinomial regression with covariate adjustment (for categorical variables)")
  } else {
    # Set caption for unadjusted tests
    test_caption <- "Statistical tests: ANOVA (for continuous variables) and Fisher's exact test (for categorical variables)"
  }

  # Add effect size if requested
  if (AddEffectSize) {
    # Create a data frame to store effect sizes
    effect_sizes <- data.frame(
      variable = character(),
      effect_size = numeric(),
      effect_type = character(),
      stringsAsFactors = FALSE
    )

    # Get the variables in the table
    table_vars <- unique(CompTable$table_body$variable)

    # Calculate effect sizes for each variable
    for (var in table_vars) {
      # Skip the comparison variable itself and covariates
      if (var == CompVariable || var %in% Covariates) {
        next
      }

      # Get the clean variable name
      clean_var <- name_mapping$clean[name_mapping$original == var]

      # Check if variable exists in clean data
      if (!(clean_var %in% names(Data_clean))) {
        next
      }

      # Different effect size calculations based on whether covariates are used
      if (!is.null(Covariates)) {
        # With covariate adjustment

        # For continuous variables
        if (is.numeric(Data_clean[[clean_var]])) {
          # Calculate adjusted effect size (partial eta-squared)
          tryCatch({
            formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                  paste(clean_Covariates, collapse = " + "))
            model <- stats::lm(stats::as.formula(formula_str), data = Data_clean)
            anova_result <- car::Anova(model, type = 2)

            # Calculate partial eta-squared
            ss_effect <- anova_result[clean_CompVariable, "Sum Sq"]
            ss_error <- anova_result["Residuals", "Sum Sq"]
            partial_eta_sq <- ss_effect / (ss_effect + ss_error)

            # Add to effect sizes data frame
            effect_sizes <- rbind(effect_sizes, data.frame(
              variable = var,
              effect_size = partial_eta_sq,
              effect_type = "partial_eta_sq",
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            warning("Error calculating adjusted effect size for ", var, ": ", e$message)
          })
        }
        # For categorical variables
        else if (is.factor(Data_clean[[clean_var]])) {
          # For binary or multi-level outcomes, calculate Nagelkerke's pseudo-R²
          tryCatch({
            # Full model with group and covariates
            formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                  paste(clean_Covariates, collapse = " + "))

            # Reduced model with only covariates
            null_formula_str <- paste0(clean_var, " ~ ", paste(clean_Covariates, collapse = " + "))

            # For binary variables
            if (length(levels(Data_clean[[clean_var]])) <= 2) {
              full_model <- stats::glm(stats::as.formula(formula_str), family = binomial(), data = Data_clean)
              reduced_model <- stats::glm(stats::as.formula(null_formula_str), family = binomial(), data = Data_clean)

              # Calculate McFadden's pseudo R²
              pseudo_r2 <- 1 - (stats::logLik(full_model) / stats::logLik(reduced_model))

              # Add to effect sizes data frame
              effect_sizes <- rbind(effect_sizes, data.frame(
                variable = var,
                effect_size = as.numeric(pseudo_r2),
                effect_type = "pseudo_r2",
                stringsAsFactors = FALSE
              ))
            }
            # For multi-level variables
            else {
              if (requireNamespace("nnet", quietly = TRUE)) {
                full_model <- nnet::multinom(stats::as.formula(formula_str), data = Data_clean, trace = FALSE)
                reduced_model <- nnet::multinom(stats::as.formula(null_formula_str), data = Data_clean, trace = FALSE)

                # Calculate McFadden's pseudo R²
                pseudo_r2 <- 1 - (stats::logLik(full_model) / stats::logLik(reduced_model))

                # Add to effect sizes data frame
                effect_sizes <- rbind(effect_sizes, data.frame(
                  variable = var,
                  effect_size = as.numeric(pseudo_r2),
                  effect_type = "pseudo_r2",
                  stringsAsFactors = FALSE
                ))
              }
            }
          }, error = function(e) {
            warning("Error calculating adjusted effect size for ", var, ": ", e$message)
          })
        }
      }
      else {
        # Without covariate adjustment - use standard effect sizes

        # For continuous variables
        if (is.numeric(Data_clean[[clean_var]])) {
          # Calculate eta-squared
          tryCatch({
            formula <- stats::as.formula(paste(clean_var, "~", clean_CompVariable))
            aov_model <- stats::aov(formula, data = Data_clean)
            result <- effectsize::eta_squared(aov_model)
            eta_sq <- result$Eta2[1]  # First row contains the effect size for the group variable

            # Add to effect sizes data frame
            effect_sizes <- rbind(effect_sizes, data.frame(
              variable = var,
              effect_size = eta_sq,
              effect_type = "eta_sq",
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            warning("Error calculating effect size for ", var, ": ", e$message)
          })
        }
        # For categorical variables
        else if (is.factor(Data_clean[[clean_var]])) {
          # Calculate Cramer's V
          tryCatch({
            # Create contingency table
            tbl <- table(Data_filtered[[var]], Data_filtered[[CompVariable]])

            # Calculate Cramer's V
            result <- effectsize::cramers_v(tbl)
            cramers_v <- result$Cramers_v

            # Add to effect sizes data frame
            effect_sizes <- rbind(effect_sizes, data.frame(
              variable = var,
              effect_size = cramers_v,
              effect_type = "cramers_v",
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            warning("Error calculating effect size for ", var, ": ", e$message)
          })
        }
      }
    }

    # Add effect size column to the table
    if (nrow(effect_sizes) > 0) {
      # Round effect sizes to specified digits
      effect_sizes$effect_size_rounded <- round(effect_sizes$effect_size, EffectSizeDigits)

      # Format effect sizes as strings with the specified number of decimal places
      effect_sizes$effect_size_fmt <- sprintf(paste0("%.", EffectSizeDigits, "f"),
                                              effect_sizes$effect_size_rounded)

      # Determine the effect size header and footnote based on the types of effect sizes
      if (!is.null(Covariates)) {
        # With covariate adjustment
        if (any(effect_sizes$effect_type == "partial_eta_sq") && any(effect_sizes$effect_type == "pseudo_r2")) {
          effect_size_header <- "**Effect Size**"
          effect_size_footnote <- "Effect size measures: Partial η² (for continuous variables) and McFadden's Pseudo-R² (for categorical variables)"
        } else if (any(effect_sizes$effect_type == "partial_eta_sq")) {
          effect_size_header <- "**Partial η²**"
          effect_size_footnote <- "Partial η²: Proportion of variance explained by group after accounting for covariates"
        } else {
          effect_size_header <- "**Pseudo-R²**"
          effect_size_footnote <- "McFadden's Pseudo-R²: Measure of improvement in model fit due to group variable"
        }
      } else {
        # Without covariate adjustment
        if (any(effect_sizes$effect_type == "eta_sq") && any(effect_sizes$effect_type == "cramers_v")) {
          effect_size_header <- "**Effect Size**"
          effect_size_footnote <- "Effect size measures: η² (for continuous variables) and Cramer's V (for categorical variables)"
        } else if (any(effect_sizes$effect_type == "eta_sq")) {
          effect_size_header <- "**η²**"
          effect_size_footnote <- "η²: Proportion of variance explained (0=none, 1=all)"
        } else {
          effect_size_header <- "**Cramer's V**"
          effect_size_footnote <- "Cramer's V: Effect size for categorical variables (0=no association, 1=perfect association)"
        }
      }

      # Add effect size column to the table - CORRECTED IMPLEMENTATION
      CompTable <- CompTable %>%
        gtsummary::add_stat(
          fns = everything() ~ function(data, variable, ...) {
            # Find the effect size for this variable
            es_row <- which(effect_sizes$variable == variable)
            if (length(es_row) > 0) {
              return(effect_sizes$effect_size_fmt[es_row])
            } else {
              return(NA_character_)
            }
          }
        ) %>%
        # Modify the header for the new column
        gtsummary::modify_header(
          add_stat_1 = effect_size_header
        ) %>%
        # Add footnote for the effect size column
        gtsummary::modify_footnote(
          add_stat_1 = effect_size_footnote
        )
    }
  }

  # Add pairwise comparisons if requested and if there are more than 2 groups
  if (AddPairwise && length(unique(DataFrame[[CompVariable]])) > 2) {
    # Get all unique groups
    groups <- unique(DataFrame[[CompVariable]])

    # Determine which pairs to compare
    if (!is.null(ReferentGroup)) {
      # Check if the reference group exists in the data
      if (!(ReferentGroup %in% groups)) {
        warning(paste("Reference group", ReferentGroup, "not found in the data. Using all pairwise comparisons instead."))
        pairs <- utils::combn(groups, 2, simplify = FALSE)
      } else {
        # Create pairs with reference group as first element
        other_groups <- groups[groups != ReferentGroup]
        pairs <- lapply(other_groups, function(g) c(ReferentGroup, g))
      }
    } else {
      # All pairwise comparisons
      pairs <- utils::combn(groups, 2, simplify = FALSE)
    }

    # Calculate total number of comparisons for p-value adjustment
    total_comparisons <- length(pairs)

    # Create column names for each pairwise comparison
    pair_names <- sapply(pairs, function(pair) paste(pair[1], "vs", pair[2]))

    # Create a list to store pairwise comparison tables
    pairwise_tables <- list()

    # Create individual comparison tables for each pair
    for (i in 1:length(pairs)) {
      pair <- pairs[[i]]
      pair_name <- pair_names[i]

      # Subset data for this pair
      pair_data <- Data_filtered[Data_filtered[[CompVariable]] %in% pair, ]
      pair_data_clean <- Data_clean[Data_clean[[clean_CompVariable]] %in% pair, ]

      # Skip if insufficient data
      if (nrow(pair_data) < 2) next

      # Create a comparison table for this pair
      suppressWarnings({
        pair_table <- gtsummary::tbl_summary(
          data = pair_data,
          by = CompVariable,
          missing = "no",
          statistic = list(
            all_continuous() ~ "{mean} ({sd})",
            all_categorical() ~ "{n} ({p})"
          ),
          digits = list(
            all_continuous() ~ ValueDigits
          ),
          type = list(
            where(is.double) ~ "continuous",
            where(~ is.factor(.) & nlevels(.) == 2) ~ "categorical"
          )
        )
      })

      # Add p-values
      suppressWarnings({
        pair_table <- pair_table %>%
          gtsummary::add_p(
            test = list(
              all_continuous() ~ "t.test",
              all_categorical() ~ "fisher.test"
            )
          )
      })

      # If covariates are specified, manually calculate adjusted p-values for pairwise comparisons
      if (!is.null(Covariates)) {
        # Get the variables in the table
        pair_vars <- unique(pair_table$table_body$variable)

        # Create a data frame to store adjusted p-values
        pair_adjusted_p_values <- data.frame(
          variable = character(),
          p_value = numeric(),
          stringsAsFactors = FALSE
        )

        # Calculate adjusted p-values for each variable
        for (var in pair_vars) {
          # Skip the comparison variable itself and covariates
          if (var == CompVariable || var %in% Covariates) {
            next
          }

          # Get the clean variable name
          clean_var <- name_mapping$clean[name_mapping$original == var]

          # Check if variable exists in clean data
          if (!(clean_var %in% names(pair_data_clean))) {
            next
          }

          # Check variable type
          if (is.numeric(pair_data_clean[[clean_var]])) {
            # For continuous variables, use linear regression
            tryCatch({
              # Create formula with clean variable names
              formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                    paste(clean_Covariates, collapse = " + "))

              # Fit the model
              model <- stats::lm(stats::as.formula(formula_str), data = pair_data_clean)

              # Get ANOVA results
              anova_result <- car::Anova(model, type = 2)

              # Extract p-value for the grouping variable
              p_value <- anova_result[clean_CompVariable, "Pr(>F)"]

              # Add to adjusted p-values data frame
              pair_adjusted_p_values <- rbind(pair_adjusted_p_values, data.frame(
                variable = var,
                p_value = p_value,
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              warning("Error calculating adjusted p-value for ", var, " in pair ", pair_name, ": ", e$message)
            })
          } else if (is.factor(pair_data_clean[[clean_var]])) {
            # For categorical variables, use logistic regression
            tryCatch({
              # Create formula with clean variable names
              formula_str <- paste0(clean_var, " ~ ", clean_CompVariable, " + ",
                                    paste(clean_Covariates, collapse = " + "))

              # Fit the model
              model <- stats::glm(stats::as.formula(formula_str), family = binomial(), data = pair_data_clean)

              # Get Wald test results
              wald_test <- car::Anova(model, type = 2)

              # Extract p-value for the grouping variable
              p_value <- wald_test[clean_CompVariable, "Pr(>Chisq)"]

              # Add to adjusted p-values data frame
              pair_adjusted_p_values <- rbind(pair_adjusted_p_values, data.frame(
                variable = var,
                p_value = p_value,
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              warning("Error calculating adjusted p-value for ", var, " in pair ", pair_name, ": ", e$message)
            })
          }
        }

        # Replace p-values in the pair table with adjusted p-values
        if (nrow(pair_adjusted_p_values) > 0) {
          # Create a copy of the table body
          new_pair_table_body <- pair_table$table_body

          # Update p-values
          for (j in 1:nrow(pair_adjusted_p_values)) {
            var <- pair_adjusted_p_values$variable[j]
            p_val <- pair_adjusted_p_values$p_value[j]

            # Find rows for this variable
            var_rows <- which(new_pair_table_body$variable == var)

            if (length(var_rows) > 0) {
              # Update p-value
              new_pair_table_body$p.value[var_rows] <- p_val
            }
          }

          # Replace the table body
          pair_table$table_body <- new_pair_table_body
        }
      }

      # Extract p-values from the pair table
      p_values <- pair_table$table_body$p.value
      variables <- pair_table$table_body$variable

      # Store the table with pair name as key
      pairwise_tables[[pair_name]] <- list(
        table = pair_table,
        p_values = p_values,
        variables = variables
      )
    }

    # Adjust p-values if needed
    if (PairwiseMethod != "none" && length(pairwise_tables) > 0) {
      # Collect all p-values with their variables
      all_p_data <- data.frame(
        pair = character(),
        variable = character(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      )

      for (pair_name in names(pairwise_tables)) {
        if (length(pairwise_tables[[pair_name]]$p_values) > 0) {
          pair_df <- data.frame(
            pair = pair_name,
            variable = pairwise_tables[[pair_name]]$variables,
            p_value = pairwise_tables[[pair_name]]$p_values,
            stringsAsFactors = FALSE
          )
          all_p_data <- rbind(all_p_data, pair_df)
        }
      }

      # Remove NA p-values
      all_p_data <- all_p_data[!is.na(all_p_data$p_value), ]

      if (nrow(all_p_data) > 0) {
        # Adjust p-values
        all_p_data$adjusted_p_value <- stats::p.adjust(all_p_data$p_value, method = PairwiseMethod)

        # Update the p-values in pairwise_tables
        for (i in 1:nrow(all_p_data)) {
          pair_name <- all_p_data$pair[i]
          var_name <- all_p_data$variable[i]
          adj_p <- all_p_data$adjusted_p_value[i]

          # Find the index of this variable in the pair's table
          idx <- which(pairwise_tables[[pair_name]]$variables == var_name)
          if (length(idx) > 0) {
            # Initialize adjusted_p_values if it doesn't exist
            if (is.null(pairwise_tables[[pair_name]]$adjusted_p_values)) {
              pairwise_tables[[pair_name]]$adjusted_p_values <- rep(NA, length(pairwise_tables[[ pair_name]]$p_values))
            }
            pairwise_tables[[pair_name]]$adjusted_p_values[idx] <- adj_p
          }
        }
      }
    }

    # Add pairwise comparison columns to the main table
    for (pair_name in names(pairwise_tables)) {
      col_name <- paste0("pairwise_", gsub(" ", "_", pair_name))

      # Get variables from the pairwise table
      pair_vars <- pairwise_tables[[pair_name]]$variables

      # Get p-values (adjusted if available, otherwise original)
      if (PairwiseMethod != "none" && !is.null(pairwise_tables[[pair_name]]$adjusted_p_values)) {
        pair_p_values <- pairwise_tables[[pair_name]]$adjusted_p_values
      } else {
        pair_p_values <- pairwise_tables[[pair_name]]$p_values
      }

      # Create a mapping between variables and p-values
      p_value_mapping <- data.frame(
        variable = pair_vars,
        p_value = pair_p_values,
        stringsAsFactors = FALSE
      )

      # Add this column to the main table
      CompTable <- CompTable %>%
        gtsummary::modify_table_body(
          ~.x %>%
            dplyr::mutate(
              !!col_name := dplyr::case_when(
                !variable %in% p_value_mapping$variable ~ NA_character_,
                is.na(p_value_mapping$p_value[match(variable, p_value_mapping$variable)]) ~ NA_character_,
                p_value_mapping$p_value[match(variable, p_value_mapping$variable)] < 0.05 ~
                  paste0("**", gtsummary::style_pvalue(p_value_mapping$p_value[match(variable, p_value_mapping$variable)], digits = pDigits), "**"),
                TRUE ~
                  gtsummary::style_pvalue(p_value_mapping$p_value[match(variable, p_value_mapping$variable)], digits = pDigits)
              )
            )
        )

      # Add header for this column
      CompTable <- CompTable %>%
        gtsummary::modify_header(
          !!col_name := paste0("**", pair_name, "**")
        )
    }

    # Add footnote for pairwise comparisons
    comparison_type <- if (!is.null(ReferentGroup)) {
      paste("Comparisons with reference group", ReferentGroup)
    } else {
      "All pairwise comparisons"
    }

    pairwise_footnote <- paste0(comparison_type, " (", PairwiseMethod, " adjusted p-values)")
    if (!is.null(Covariates)) {
      pairwise_footnote <- paste0(pairwise_footnote, ", adjusted for: ", paste(Covariates, collapse = ", "))
    }
    pairwise_footnote <- paste0(pairwise_footnote, ". Bold values indicate p < 0.05.")

    CompTable <- CompTable %>%
      gtsummary::modify_footnote(
        starts_with("pairwise_") ~ pairwise_footnote
      )

    # Update test caption to include pairwise tests
    if (is.null(Covariates)) {
      test_caption <- paste0(test_caption, "; Pairwise comparisons: t-tests (for continuous variables) and Fisher's exact test (for categorical variables)")
    } else {
      test_caption <- paste0(test_caption, "; Pairwise comparisons: Linear models (for continuous variables) and logistic regression (for categorical variables) with covariate adjustment")
    }
  }

  # Make sure p-values are properly formatted
  CompTable <- CompTable %>%
    gtsummary::modify_fmt_fun(
      p.value ~ function(x) gtsummary::style_pvalue(x, digits = pDigits)
    )

  # Add the caption about statistical tests used
  CompTable <- CompTable %>%
    gtsummary::modify_caption(test_caption)

  return(CompTable)
}
