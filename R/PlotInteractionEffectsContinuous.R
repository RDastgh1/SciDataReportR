#' Plot Single Interaction Effect
#'
#' Creates a scatter plot with regression lines showing the interaction between a predictor
#' and outcome variable, moderated by either a continuous or categorical variable.
#'
#' @param Data A data frame containing the variables to be analyzed
#' @param interVar Character string specifying the interaction variable (moderator)
#' @param outcomeVar Character string specifying the outcome variable
#' @param predictorVar Character string specifying the predictor variable
#' @param covars Character vector of covariate names to include in the model
#' @param n_lines For continuous moderators, number of lines to plot (default: 3 for low/med/high)
#' @param alpha Transparency level for points (default: 0.6)
#' @param point_size Size of points (default: 2)
#'
#' @return A ggplot object showing the interaction effect
#'
#' @details
#' For categorical moderators, separate regression lines are plotted for each category.
#' For continuous moderators, regression lines are plotted at the mean and ± 1 SD.
#'
#' The subtitle shows the p-value for the interaction term.
#' The caption lists any covariates included in the model.
#'
#' Variable labels are used if available in the data frame.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @importFrom dplyr mutate case_when
#' @importFrom stats lm coef sd quantile
#' @importFrom paletteer paletteer_d paletteer_c

PlotInteractionEffectsContinuous <- function(Data,
                                  interVar = NULL,
                                  outcomeVar = NULL,
                                  predictorVar = NULL,
                                  covars = NULL,
                                  n_lines = 3,
                                  alpha = 0.6,
                                  point_size = 2) {

  # Input validation
  if (is.null(Data) || !is.data.frame(Data)) {
    stop("Data must be a non-null data frame")
  }

  if (is.null(interVar) || is.null(outcomeVar) || is.null(predictorVar)) {
    stop("interVar, outcomeVar, and predictorVar must all be specified")
  }

  if (!all(c(interVar, outcomeVar, predictorVar) %in% names(Data))) {
    stop("All specified variables must exist in the data frame")
  }

  if (!is.null(covars) && !all(covars %in% names(Data))) {
    stop("All covariates must exist in the data frame")
  }

  # Load required packages
  library(ggplot2)
  library(dplyr)

  # Check if paletteer is available
  use_paletteer <- requireNamespace("paletteer", quietly = TRUE)

  # Function to get variable label if it exists
  get_var_label <- function(data, var_name) {
    if (requireNamespace("sjlabelled", quietly = TRUE)) {
      label <- sjlabelled::get_label(data[[var_name]])
      if (!is.null(label) && label != "") {
        return(label)
      }
    }
    return(var_name)
  }

  # Get labels for all variables
  outcome_label <- get_var_label(Data, outcomeVar)
  predictor_label <- get_var_label(Data, predictorVar)
  inter_label <- get_var_label(Data, interVar)

  # Check if interaction variable is continuous
  interVar_is_continuous <- is.numeric(Data[[interVar]])

  # Build formula for the model
  if (is.null(covars)) {
    formula_str <- paste0("`", outcomeVar, "` ~ `", predictorVar, "` * `", interVar, "`")
  } else {
    covar_terms <- paste(paste0("`", covars, "`"), collapse = " + ")
    formula_str <- paste0("`", outcomeVar, "` ~ ", covar_terms, " + `",
                          predictorVar, "` * `", interVar, "`")
  }

  # Fit the model
  model <- lm(as.formula(formula_str), data = Data)
  model_summary <- summary(model)

  # Find the interaction term p-value
  coef_names <- rownames(model_summary$coefficients)
  interaction_term <- NULL

  # Look for interaction term
  for (coef_name in coef_names) {
    if (grepl(predictorVar, coef_name, fixed = TRUE) &&
        grepl(interVar, coef_name, fixed = TRUE) &&
        grepl(":", coef_name, fixed = TRUE)) {
      interaction_term <- coef_name
      break
    }
  }

  if (is.null(interaction_term)) {
    stop("Could not find interaction term in model")
  }

  # Get interaction p-value
  interaction_p <- model_summary$coefficients[interaction_term, "Pr(>|t|)"]

  # Prepare data for plotting
  plot_data <- Data[, c(outcomeVar, predictorVar, interVar)]
  names(plot_data) <- c("outcome", "predictor", "moderator")

  # Remove rows with missing values
  plot_data <- plot_data[complete.cases(plot_data), ]

  # Create the plot based on moderator type
  if (interVar_is_continuous) {
    # For continuous moderator, create categories for visualization
    mod_mean <- mean(plot_data$moderator, na.rm = TRUE)
    mod_sd <- sd(plot_data$moderator, na.rm = TRUE)

    if (n_lines == 3) {
      # Create three groups: Low (-1SD), Medium (mean), High (+1SD)
      plot_data <- plot_data %>%
        mutate(
          mod_group = case_when(
            moderator <= (mod_mean - mod_sd) ~ paste0(inter_label, " = ", round(mod_mean - mod_sd, 1), " (-1 SD)"),
            moderator >= (mod_mean + mod_sd) ~ paste0(inter_label, " = ", round(mod_mean + mod_sd, 1), " (+1 SD)"),
            TRUE ~ paste0(inter_label, " = ", round(mod_mean, 1), " (Mean)")
          ),
          mod_group = factor(mod_group, levels = c(
            paste0(inter_label, " = ", round(mod_mean - mod_sd, 1), " (-1 SD)"),
            paste0(inter_label, " = ", round(mod_mean, 1), " (Mean)"),
            paste0(inter_label, " = ", round(mod_mean + mod_sd, 1), " (+1 SD)")
          ))
        )
    } else {
      # Create quantile-based groups
      quantiles <- quantile(plot_data$moderator, probs = seq(0, 1, length.out = n_lines + 1))
      plot_data <- plot_data %>%
        mutate(
          mod_group = cut(moderator, breaks = quantiles, include.lowest = TRUE,
                          labels = paste0(inter_label, " Q", 1:n_lines))
        )
    }

    # Create plot for continuous moderator with divergent color scheme
    p <- ggplot(plot_data, aes(x = predictor, y = outcome, color = mod_group, fill = mod_group)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2)

    # Apply divergent color scheme for continuous moderator
    if (use_paletteer) {
      if (n_lines == 3) {
        # Use a divergent palette for 3 levels
        colors <- c("#6874A9FF", "black", "#8A2870FF")
      } else {
        colors <- paletteer::paletteer_c("scico::roma", n = n_lines)
      }
      p <- p +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors)
    } else {
      # Fallback colors if paletteer not available
      p <- p +
        scale_color_manual(values = c("#2166ac", "black", "#b2182b")[1:n_lines]) +
        scale_fill_manual(values = c("#2166ac", "black", "#b2182b")[1:n_lines])
    }

  } else {
    # For categorical moderator
    plot_data$mod_group <- as.factor(plot_data$moderator)
    n_groups <- length(unique(plot_data$mod_group))

    # Create plot for categorical moderator
    p <- ggplot(plot_data, aes(x = predictor, y = outcome, color = mod_group, fill = mod_group)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2)

    # Apply categorical color scheme
    if (use_paletteer) {
      if (n_groups <= 4) {
        colors <- paletteer::paletteer_d("colorBlindness::paletteMartin")[1:n_groups]
      } else if (n_groups <= 8) {
        colors <- paletteer::paletteer_d("rcartocolor::Vivid")[1:n_groups]
      } else {
        colors <- paletteer::paletteer_d("pals::polychrome")[1:n_groups]
      }
      p <- p +
        scale_color_manual(values = colors, name = inter_label) +
        scale_fill_manual(values = colors, name = inter_label)
    } else {
      # Fallback to Set1 if paletteer not available
      p <- p +
        scale_color_brewer(palette = "Set1", name = inter_label) +
        scale_fill_brewer(palette = "Set1", name = inter_label)
    }
  }

  # Add labels and theme
  p <- p +
    labs(
      x = predictor_label,
      y = outcome_label,
      title = paste("Interaction with", inter_label),
      subtitle = paste("Interaction p-value:",
                       ifelse(interaction_p < 0.001, "< 0.001",
                              round(interaction_p, 4))),
      color = inter_label,
      fill = inter_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Add caption about covariates
  if (!is.null(covars)) {
    covar_labels <- sapply(covars, function(x) get_var_label(Data, x))
    caption_text <- paste("Model adjusted for:", paste(covar_labels, collapse = ", "))
  } else {
    caption_text <- "No covariates included in the model"
  }

  p <- p + labs(caption = caption_text)

  # Add model statistics as attributes
  attr(p, "model") <- model
  attr(p, "interaction_p") <- interaction_p
  attr(p, "formula") <- formula_str

  return(p)
}
