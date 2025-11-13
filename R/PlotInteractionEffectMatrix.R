#' Plot Interaction Effects Matrix
#'
#' Creates a heatmap visualization of interaction effects between continuous variables,
#' showing whether interactions result in slope reversals or maintain the same direction.
#'
#' @param Data A data frame containing the variables to be analyzed
#' @param interVar Character string specifying the interaction variable (moderator). Can be categorical or continuous.
#' @param outcomeVars Character vector of outcome variable names (displayed on rows)
#' @param predictorVars Character vector of predictor variable names (displayed on columns)
#' @param covars Character vector of covariate names to include in the models
#' @param Relabel Logical indicating whether to use variable labels if available (default: TRUE)
#' @param Ordinal Logical indicating whether to treat ordered factors as numeric (default: FALSE)
#'
#' @return A list containing:
#' \describe{
#'   \item{Unadjusted}{List with unadjusted results including:
#'     \itemize{
#'       \item C: Matrix of interaction coefficients
#'       \item S: Matrix of slope direction indicators (1 = same, -1 = reversed)
#'       \item P: Matrix of p-values
#'       \item D: Matrix of interaction coefficient signs
#'       \item Slope1: Matrix of slopes at low values of interVar (mean - 1SD for continuous, reference group for categorical)
#'       \item Slope2: Matrix of slopes at high values of interVar (mean + 1SD for continuous, comparison group for categorical)
#'       \item plot: ggplot object of the heatmap
#'       \item pvaltable: P-value table in wide format
#'     }
#'   }
#'   \item{FDRCorrected}{List with FDR-corrected results (same structure as Unadjusted)}
#'   \item{Relabel}{Logical indicating whether relabeling was applied}
#'   \item{Covariates}{Character vector of covariates used}
#'   \item{interVar}{The interaction variable name}
#'   \item{raw_data}{Processed data frame with all calculated values}
#' }
#'
#' @details
#' The function fits linear models of the form: outcome ~ predictor * interVar + covariates
#' for each combination of outcome and predictor variables.
#'
#' For continuous interaction variables, slopes are calculated at mean - 1SD and mean + 1SD.
#' For categorical interaction variables, slopes are calculated for each category.
#'
#' The function then determines whether the interaction causes a slope reversal (opposite signs)
#' or maintains the same direction (same signs) for the predictor-outcome relationship.
#'
#' Color coding:
#' \itemize{
#'   \item Blue gradient: Significant interaction with slopes in the same direction
#'   \item Red gradient: Significant interaction with slope reversal
#'   \item White: Non-significant interaction (p > 0.05)
#'   \item Grey: Missing data or model could not be fit
#' }
#'
#' Darker colors indicate higher significance:
#' \itemize{
#'   \item ***: p ≤ 0.001 (darkest)
#'   \item **: p ≤ 0.01 (medium)
#'   \item *: p ≤ 0.05 (light)
#' }
#'
#' @examples
#' \dontrun{
#' # With categorical interaction variable
#' results <- PlotInteractionEffectsMatrix(
#'   Data = mydata,
#'   interVar = "HIV_status",
#'   outcomeVars = sleep_vars,
#'   predictorVars = metabolites,
#'   covars = c("age", "sex")
#' )
#'
#' # With continuous interaction variable
#' results <- PlotInteractionEffectsMatrix(
#'   Data = mydata,
#'   interVar = "age",
#'   outcomeVars = sleep_vars,
#'   predictorVars = metabolites,
#'   covars = c("sex", "BMI")
#' )
#'
#' # Display the plot
#' print(results$Unadjusted$plot)
#' }
#'
#' @export
#' @importFrom dplyr %>% left_join mutate select rename filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_manual theme_minimal theme element_text element_blank labs
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats lm coef p.adjust as.formula sd

PlotInteractionEffectsMatrix <- function(Data, interVar = NULL,
                                         outcomeVars = NULL, predictorVars = NULL, covars = NULL,
                                         Relabel = TRUE, Ordinal = FALSE) {

  # Check for required packages
  required_packages <- c("dplyr", "tidyr", "ggplot2", "tibble")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Missing required packages:", paste(missing_packages, collapse = ", ")))
  }

  # Load required functions
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)

  # Create a custom function to handle p-values with NAs
  get_significance_stars <- function(p_values) {
    stars <- character(length(p_values))
    stars[is.na(p_values)] <- ""
    stars[!is.na(p_values) & p_values <= 0.001] <- "***"
    stars[!is.na(p_values) & p_values > 0.001 & p_values <= 0.01] <- "**"
    stars[!is.na(p_values) & p_values > 0.01 & p_values <= 0.05] <- "*"
    stars[!is.na(p_values) & p_values > 0.05] <- ""
    return(stars)
  }

  # Define helper functions if they don't exist
  if (!exists("ReplaceMissingLabels")) {
    ReplaceMissingLabels <- function(data) {
      return(data)
    }
  }

  if (!exists("getNumVars")) {
    getNumVars <- function(data, Ordinal = FALSE) {
      if (Ordinal) {
        vars <- names(data)[sapply(data, function(x) is.numeric(x) || is.ordered(x))]
      } else {
        vars <- names(data)[sapply(data, is.numeric)]
      }
      return(vars)
    }
  }

  if (!exists("ConvertOrdinalToNumeric")) {
    ConvertOrdinalToNumeric <- function(data, variables) {
      for (var in variables) {
        if (is.ordered(data[[var]])) {
          data[[var]] <- as.numeric(data[[var]])
        }
      }
      return(data)
    }
  }

  # Input validation
  if (is.null(Data) || !is.data.frame(Data)) {
    stop("Data must be a non-null data frame")
  }

  if (is.null(interVar)) {
    stop("interVar must be specified")
  }

  if (!interVar %in% names(Data)) {
    stop(paste("interVar", interVar, "not found in Data"))
  }

  # Check if interaction variable is continuous
  interVar_is_continuous <- is.numeric(Data[[interVar]])

  # Check for missing values and replace with appropriate labels
  removediag <- FALSE
  if(is.null(predictorVars)){
    predictorVars = outcomeVars
    removediag <- TRUE
  }

  # Replace missing labels
  Data <- ReplaceMissingLabels(Data)

  # If outcomeVars is null, get all numeric variables
  if (is.null(outcomeVars)) {
    outcomeVars <- getNumVars(Data, Ordinal = FALSE)
    if (Ordinal) {
      outcomeVars <- getNumVars(Data, Ordinal = TRUE)
    }
  }

  # Combine variables into a single vector
  Variables <- c(interVar, outcomeVars, predictorVars)

  # If Ordinal is TRUE, convert ordinal variables to numeric and update Variables
  if (Ordinal) {
    Variables <- c(outcomeVars, predictorVars)
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # Update outcomeVars and predictorVars to exclude covariates
  if (!is.null(covars)) {
    outcomeVars <- outcomeVars[!outcomeVars %in% covars]
    predictorVars <- predictorVars[!predictorVars %in% covars]
  }

  # Initialize matrices to store results
  r_P <- r_C <- r_S <- r_D <- matrix(NA, nrow = length(outcomeVars), ncol = length(predictorVars))
  r_Slope1 <- r_Slope2 <- matrix(NA, nrow = length(outcomeVars), ncol = length(predictorVars))

  # Loop through each combination of outcome and predictor variables
  for (i in seq_len(length(outcomeVars))) {
    for (j in seq_len(length(predictorVars))) {
      outcomeVar <- outcomeVars[i]
      predictorVar <- predictorVars[j]

      # Skip if outcomeVar equals predictorVar and removediag is TRUE
      if(removediag && outcomeVar == predictorVar){
        next
      }

      # Convert character variables to factor
      if(is.character(Data[[outcomeVar]])){
        Data[[outcomeVar]] <- as.factor(Data[[outcomeVar]])
      }

      if(is.character(Data[[predictorVar]])){
        Data[[predictorVar]] <- as.factor(Data[[predictorVar]])
      }

      # Try to fit linear model, catching errors
      tryCatch({
        # Build formula with outcomeVar as outcome and predictorVar as predictor
        if(is.null(covars)){
          formula_str <- paste0("`", outcomeVar, "` ~ `", predictorVar, "` * `", interVar, "`")
        } else {
          covar_terms <- paste(paste0("`", covars, "`"), collapse = " + ")
          formula_str <- paste0("`", outcomeVar, "` ~ ", covar_terms, " + `", predictorVar, "` * `", interVar, "`")
        }

        formula <- as.formula(formula_str)
        m <- lm(formula, data = Data)

        # Get coefficient names and summary
        coef_names <- names(coef(m))
        model_summary <- summary(m)

        # Find main effect term - look for predictorVar
        main_effect_term <- paste0("`", predictorVar, "`")
        if (!main_effect_term %in% coef_names) {
          # Try without backticks
          main_effect_term <- predictorVar
        }

        # Find interaction term - need to check multiple patterns
        interaction_patterns <- c(
          paste0("`", predictorVar, "`:`", interVar, "`"),
          paste0("`", interVar, "`:`", predictorVar, "`"),
          paste0(predictorVar, ":", interVar),
          paste0(interVar, ":", predictorVar)
        )

        # For categorical interVar, also check for level-specific patterns
        if (!interVar_is_continuous) {
          levels_interVar <- levels(as.factor(Data[[interVar]]))
          if (length(levels_interVar) > 1) {
            for (level in levels_interVar[-1]) {  # Skip reference level
              interaction_patterns <- c(interaction_patterns,
                                        paste0("`", predictorVar, "`:", interVar, level),
                                        paste0(interVar, level, ":`", predictorVar, "`"))
            }
          }
        }

        interaction_term <- NULL
        for (pattern in interaction_patterns) {
          matches <- grep(pattern, coef_names, value = TRUE, fixed = FALSE)
          if (length(matches) > 0) {
            interaction_term <- matches[1]
            break
          }
        }

        if (!is.null(interaction_term) && main_effect_term %in% coef_names) {
          # Get coefficients
          main_effect <- coef(m)[main_effect_term]
          interC <- coef(m)[interaction_term]
          interP <- model_summary$coefficients[interaction_term, "Pr(>|t|)"]

          # Calculate slopes based on whether interVar is continuous or categorical
          if (interVar_is_continuous) {
            # For continuous interVar, calculate slopes at mean ± 1SD
            interVar_mean <- mean(Data[[interVar]], na.rm = TRUE)
            interVar_sd <- sd(Data[[interVar]], na.rm = TRUE)

            # Get the main effect of interVar if it exists
            interVar_main_term <- paste0("`", interVar, "`")
            if (!interVar_main_term %in% coef_names) {
              interVar_main_term <- interVar
            }

            # Slope at low value (mean - 1SD)
            slope1 <- main_effect + interC * (interVar_mean - interVar_sd)
            # Slope at high value (mean + 1SD)
            slope2 <- main_effect + interC * (interVar_mean + interVar_sd)
          } else {
            # For categorical interVar, slopes for reference and comparison groups
            slope1 <- main_effect  # Reference group
            slope2 <- main_effect + interC  # Other group
          }

          # Store slopes
          r_Slope1[i, j] <- slope1
          r_Slope2[i, j] <- slope2

          # Determine if slopes have same or opposite signs
          same_direction <- (sign(slope1) == sign(slope2))

          # Direction based on interaction coefficient sign
          interD <- sign(interC)

          # Slope reversal indicator
          if (same_direction) {
            interS <- 1  # Same direction
          } else {
            interS <- -1  # Reversed direction
          }

        } else {
          interC <- NA
          interP <- NA
          interS <- NA
          interD <- NA
        }

        # Store results in matrices
        r_P[i, j] <- interP
        r_C[i, j] <- interC
        r_S[i, j] <- interS
        r_D[i, j] <- interD

      }, error = function(err) {
        # If an error occurred, set results to NA
        r_P[i, j] <- NA
        r_C[i, j] <- NA
        r_S[i, j] <- NA
        r_D[i, j] <- NA
        r_Slope1[i, j] <- NA
        r_Slope2[i, j] <- NA
      })
    }
  }

  rownames(r_C) <- outcomeVars
  colnames(r_C) <- predictorVars
  rownames(r_P) <- outcomeVars
  colnames(r_P) <- predictorVars
  rownames(r_S) <- outcomeVars
  colnames(r_S) <- predictorVars
  rownames(r_D) <- outcomeVars
  colnames(r_D) <- predictorVars
  rownames(r_Slope1) <- outcomeVars
  colnames(r_Slope1) <- predictorVars
  rownames(r_Slope2) <- outcomeVars
  colnames(r_Slope2) <- predictorVars

  # Convert matrices to data frames and add variable names
  m_r_C <- r_C %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "C")

  m_r_P <- r_P %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "P")

  m_r_S <- r_S %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "S")

  m_r_D <- r_D %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "D")

  m_r_Slope1 <- r_Slope1 %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "Slope1")

  m_r_Slope2 <- r_Slope2 %>%
    as.data.frame() %>%
    rownames_to_column(var = "Outcome") %>%
    pivot_longer(cols = all_of(predictorVars), names_to = "Predictor", values_to = "Slope2")

  # Join all data frames
  m_G <- left_join(m_r_C, m_r_P, by = c("Outcome", "Predictor")) %>%
    left_join(m_r_S, by = c("Outcome", "Predictor")) %>%
    left_join(m_r_D, by = c("Outcome", "Predictor")) %>%
    left_join(m_r_Slope1, by = c("Outcome", "Predictor")) %>%
    left_join(m_r_Slope2, by = c("Outcome", "Predictor"))

  # Create significance and direction categories
  m_G$sign <- "ns"  # Default to non-significant
  m_G$sign[is.na(m_G$P)] <- "na"  # Set NA values

  # Use custom function for significance stars
  m_G$sig <- get_significance_stars(m_G$P)

  # For significant interactions, determine direction and create gradient categories
  sig_mask <- !is.na(m_G$P) & m_G$P <= 0.05

  # Same direction (blue gradient)
  same_sign_mask <- sig_mask & m_G$S == 1
  m_G$sign[same_sign_mask & m_G$P <= 0.001] <- "same ***"
  m_G$sign[same_sign_mask & m_G$P > 0.001 & m_G$P <= 0.01] <- "same **"
  m_G$sign[same_sign_mask & m_G$P > 0.01 & m_G$P <= 0.05] <- "same *"

  # Reversed direction (red gradient)
  diff_sign_mask <- sig_mask & m_G$S == -1
  m_G$sign[diff_sign_mask & m_G$P <= 0.001] <- "reversed ***"
  m_G$sign[diff_sign_mask & m_G$P > 0.001 & m_G$P <= 0.01] <- "reversed **"
  m_G$sign[diff_sign_mask & m_G$P > 0.01 & m_G$P <= 0.05] <- "reversed *"

  # Factor with proper levels
  m_G$sign <- factor(m_G$sign,
                     levels = c("reversed ***", "reversed **", "reversed *",
                                "ns",
                                "same *", "same **", "same ***",
                                "na"))

  # Create combined sign and significance
  m_G$sigsign <- paste(as.character(m_G$sign), as.character(m_G$sig))
  m_G$sigsign <- trimws(m_G$sigsign)

  ## FDR correction for r_P
  p_values_for_adjustment <- m_G$P[!is.na(m_G$P)]
  if (length(p_values_for_adjustment) > 0) {
    adjusted_p_values <- p.adjust(p_values_for_adjustment, method = "fdr")
    m_G$P_FDR <- NA
    m_G$P_FDR[!is.na(m_G$P)] <- adjusted_p_values
  } else {
    m_G$P_FDR <- NA
  }

  # Repeat for FDR corrected values
  m_G$sign_FDR <- "ns"
  m_G$sign_FDR[is.na(m_G$P_FDR)] <- "na"

  sig_mask_fdr <- !is.na(m_G$P_FDR) & m_G$P_FDR <= 0.05

  # Same direction (blue gradient) - FDR
  same_sign_mask_fdr <- sig_mask_fdr & m_G$S == 1
  m_G$sign_FDR[same_sign_mask_fdr & m_G$P_FDR <= 0.001] <- "same ***"
  m_G$sign_FDR[same_sign_mask_fdr & m_G$P_FDR > 0.001 & m_G$P_FDR <= 0.01] <- "same **"
  m_G$sign_FDR[same_sign_mask_fdr & m_G$P_FDR > 0.01 & m_G$P_FDR <= 0.05] <- "same *"

  # Reversed direction (red gradient) - FDR
  diff_sign_mask_fdr <- sig_mask_fdr & m_G$S == -1
  m_G$sign_FDR[diff_sign_mask_fdr & m_G$P_FDR <= 0.001] <- "reversed ***"
  m_G$sign_FDR[diff_sign_mask_fdr & m_G$P_FDR > 0.001 & m_G$P_FDR <= 0.01] <- "reversed **"
  m_G$sign_FDR[diff_sign_mask_fdr & m_G$P_FDR > 0.01 & m_G$P_FDR <= 0.05] <- "reversed *"

  m_G$sign_FDR <- factor(m_G$sign_FDR,
                         levels = c("reversed ***", "reversed **", "reversed *",
                                    "ns",
                                    "same *", "same **", "same ***",
                                    "na"))

  m_G$sig_FDR <- get_significance_stars(m_G$P_FDR)

  ## add labels
  if (Relabel && requireNamespace("sjlabelled", quietly = TRUE)) {
    Data <- ReplaceMissingLabels(Data)

    outcomelabels <- data.frame(
      Variable = unique(m_G$Outcome),
      label = sjlabelled::get_label(Data[unique(m_G$Outcome)], def.value = unique(m_G$Outcome))
    )

    predictorlabels <- data.frame(
      Variable = unique(m_G$Predictor),
      label = sjlabelled::get_label(Data[unique(m_G$Predictor)], def.value = unique(m_G$Predictor))
    )

    m_G <- m_G %>%
      left_join(outcomelabels, by = c("Outcome" = "Variable")) %>%
      rename(OutcomeLabel = label) %>%
      left_join(predictorlabels, by = c("Predictor" = "Variable")) %>%
      rename(PredictorLabel = label)
  } else {
    m_G$OutcomeLabel <- m_G$Outcome
    m_G$PredictorLabel <- m_G$Predictor
  }

  # Clean labels - remove suffixes
  m_G$OutcomeLabel <- gsub("_hr_w|_mean_w|_w", "", m_G$OutcomeLabel)
  m_G$PredictorLabel <- gsub("_log_w|_w", "", m_G$PredictorLabel)

  m_G$OutcomeLabel <- factor(m_G$OutcomeLabel, ordered = FALSE,
                             levels = rev(unique(m_G$OutcomeLabel)))
  m_G$PredictorLabel <- factor(m_G$PredictorLabel, ordered = FALSE,
                               levels = unique(m_G$PredictorLabel))

  # Create plot text
  slope_description <- ifelse(interVar_is_continuous,
                              paste0("Slope at ", interVar, " -1SD: ", round(m_G$Slope1, 3),
                                     "\nSlope at ", interVar, " +1SD: ", round(m_G$Slope2, 3)),
                              paste0("Slope 1: ", round(m_G$Slope1, 3),
                                     "\nSlope 2: ", round(m_G$Slope2, 3)))

  m_G$PlotText <- paste("Outcome:", m_G$Outcome, "\nPredictor:", m_G$Predictor,
                        "\n", slope_description,
                        "\nInteraction Coef:", round(m_G$C, 4),
                        "\nP-Value:", round(m_G$P, 4), m_G$sig,
                        "\nFDR-corrected P:", round(m_G$P_FDR, 4), m_G$sig_FDR)

  # Define color gradients
  colors <- c("same ***" = "#08519c",      # Darkest blue
              "same **" = "#3182bd",       # Medium blue
              "same *" = "#6baed6",        # Light blue
              "ns" = "#FFFFFF",             # White
              "reversed *" = "#fc9272",     # Light red
              "reversed **" = "#de2d26",    # Medium red
              "reversed ***" = "#a50f15",   # Darkest red
              "na" = "#808080")             # Grey for NA

  # Create ggplot objects with gradient color scheme
  p <- m_G %>%
    ggplot(aes(x = PredictorLabel, y = OutcomeLabel, fill = sign)) +
    geom_tile(aes(text = PlotText), show.legend = TRUE) +
    geom_text(aes(label = sig), size = 4) +
    scale_fill_manual(values = colors,
                      drop = FALSE,
                      name = "Interaction\nEffect",
                      labels = c("same ***" = "+ ***",
                                 "same **" = "+ **",
                                 "same *" = "+ *",
                                 "ns" = "ns",
                                 "reversed *" = "- *",
                                 "reversed **" = "- **",
                                 "reversed ***" = "- ***",
                                 "na" = "NA"),
                      breaks = c("same ***", "same **", "same *",
                                 "ns", "na",
                                 "reversed *", "reversed **", "reversed ***")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.caption = element_text(hjust = 0.5, size = 10, face = "italic")) +
    labs(title = paste("Interaction Effects with", interVar),
         caption = paste("No multiple comparison correction.",
                         ifelse(!is.null(covars),
                                paste("Covariates:", paste(covars, collapse = ", ")),
                                "No covariates included.")),
         x = "", y = "")

  p_FDR <- m_G %>%
    ggplot(aes(x = PredictorLabel, y = OutcomeLabel, fill = sign_FDR)) +
    geom_tile(show.legend = TRUE) +
    geom_text(aes(label = sig_FDR), size = 4) +
    scale_fill_manual(values = colors,
                      drop = FALSE,
                      name = "Interaction\nEffect",
                      labels = c("same ***" = "+ ***",
                                 "same **" = "+ **",
                                 "same *" = "+ *",
                                 "ns" = "ns",
                                 "reversed *" = "- *",
                                 "reversed **" = "- **",
                                 "reversed ***" = "- ***",
                                 "na" = "NA"),
                      breaks = c("same ***", "same **", "same *",
                                 "ns", "na",
                                 "reversed *", "reversed **", "reversed ***")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.caption = element_text(hjust = 0.5, size = 10, face = "italic")) +
    labs(title = paste("Interaction Effects with", interVar),
         caption = paste("FDR corrected for multiple comparisons.",
                         ifelse(!is.null(covars),
                                paste("Covariates:", paste(covars, collapse = ", ")),
                                "No covariates included.")),
         x = "", y = "")

  # Create p-value tables
  pvaltable <- m_G %>%
    select(Outcome, Predictor, P) %>%
    pivot_wider(names_from = Outcome, values_from = P)

  pvaltable_FDR <- m_G %>%
    select(Outcome, Predictor, P_FDR) %>%
    pivot_wider(names_from = Outcome, values_from = P_FDR)

  # Create output list
  M <- list()
  M$C <- r_C
  M$S <- r_S
  M$P <- r_P
  M$D <- r_D
  M$Slope1 <- r_Slope1
  M$Slope2 <- r_Slope2
  M$plot <- p
  M$pvaltable <- pvaltable

  M_FDR <- list()
  M_FDR$C <- r_C
  M_FDR$S <- r_S
  M_FDR$P <- r_P
  M_FDR$D <- r_D
  M_FDR$Slope1 <- r_Slope1
  M_FDR$Slope2 <- r_Slope2
  M_FDR$P_adjusted <- m_G %>%
    select(Outcome, Predictor, P_FDR) %>%
    pivot_wider(names_from = Predictor, values_from = P_FDR) %>%
    column_to_rownames("Outcome") %>%
    as.matrix()
  M_FDR$plot <- p_FDR
  M_FDR$pvaltable <- pvaltable_FDR

  # Return output list
  return(list(
    Unadjusted = M,
    FDRCorrected = M_FDR,
    Relabel = Relabel,
    Covariates = covars,
    interVar = interVar,
    interVar_type = ifelse(interVar_is_continuous, "continuous", "categorical"),
    raw_data = m_G
  ))
}
