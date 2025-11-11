#' PlotInteractionEffectsContinuous
#'
#' Plots the interaction effects between continuous variables in a dataset.
#'
#' @param Data The dataset to be analyzed
#' @param interVar The variable of interest for the interaction effect
#' @param xVars The variables to be plotted on the x-axis
#' @param yVars The variables to be plotted on the y-axis
#' @param covars Covariates to be included in the analysis
#' @param Relabel A logical value indicating whether to relabel the variables with their corresponding labels
#' @param Ordinal A logical value indicating whether the variables are ordinal
#'
#' @return A list containing the results of the analysis, including the p-values, significance levels, and the ggplot object.
#' @export

PlotInteractionEffectsContinuous <- function(Data, interVar = NULL,
                                             xVars = NULL, yVars = NULL, covars = NULL,
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
      # Simple implementation - just return the data as is
      # You should replace this with your actual implementation
      return(data)
    }
  }

  if (!exists("getNumVars")) {
    getNumVars <- function(data, Ordinal = FALSE) {
      # Get numeric variables from the dataset
      if (Ordinal) {
        # Include ordered factors as well
        vars <- names(data)[sapply(data, function(x) is.numeric(x) || is.ordered(x))]
      } else {
        vars <- names(data)[sapply(data, is.numeric)]
      }
      return(vars)
    }
  }

  if (!exists("ConvertOrdinalToNumeric")) {
    ConvertOrdinalToNumeric <- function(data, variables) {
      # Convert ordered factors to numeric
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

  # Check for missing values and replace with appropriate labels
  removediag <- FALSE
  if(is.null(yVars)){
    yVars = xVars
    removediag <- TRUE
  }

  # Replace missing labels
  Data <- ReplaceMissingLabels(Data)

  # If xVars is null, get all numeric variables
  if (is.null(xVars)) {
    xVars <- getNumVars(Data, Ordinal = FALSE)
    if (Ordinal) {
      xVars <- getNumVars(Data, Ordinal = TRUE)
    }
  }

  # Combine variables into a single vector
  Variables <- c(interVar, xVars, yVars)

  # If Ordinal is TRUE, convert ordinal variables to numeric and update Variables
  if (Ordinal) {
    Variables <- c(xVars, yVars)
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # Update xVars and yVars to exclude covariates
  if (!is.null(covars)) {
    xVars <- xVars[!xVars %in% covars]
    yVars <- yVars[!yVars %in% covars]
  }

  # Initialize matrices to store results
  r_P <- r_C <- r_S <- r_D <- matrix(NA, nrow = length(xVars), ncol = length(yVars))

  # Loop through each combination of x and y variables
  for (i in seq_len(length(xVars))) {
    for (j in seq_len(length(yVars))) {
      xVar <- xVars[i]
      yVar <- yVars[j]

      # Skip if xVar equals yVar and removediag is TRUE
      if(removediag && xVar == yVar){
        r_C[i, j] <- NA
        r_P[i, j] <- NA
        r_S[i, j] <- NA
        r_D[i, j] <- NA
        next
      }

      # Convert character variables to factor
      if(is.character(Data[[xVar]])){
        Data[[xVar]] <- as.factor(Data[[xVar]])
      }

      if(is.character(Data[[yVar]])){
        Data[[yVar]] <- as.factor(Data[[yVar]])
      }

      # Try to fit linear model, catching errors
      tryCatch({
        # Build formula properly
        if(is.null(covars)){
          formula_str <- paste(yVar, "~", xVar, "*", interVar)
        } else {
          covar_terms <- paste(covars, collapse = " + ")
          formula_str <- paste(yVar, "~", covar_terms, "+", xVar, "*", interVar)
        }

        formula <- as.formula(formula_str)
        m <- lm(formula, data = Data)

        # Summarize model and extract results
        d <- summary(m)
        dd <- as.data.frame(d$coefficients)

        # Get the interaction term (last row)
        interaction_term <- paste0(xVar, ":", interVar)
        interaction_term_alt <- paste0(interVar, ":", xVar)

        # Find the interaction row
        interaction_row <- which(rownames(dd) == interaction_term | rownames(dd) == interaction_term_alt)

        if (length(interaction_row) > 0) {
          interC <- dd$Estimate[interaction_row]
          interP <- dd$`Pr(>|t|)`[interaction_row]

          # Calculate the direction of the interaction effect
          # Positive interaction: the effect of X on Y increases as interVar increases
          # Negative interaction: the effect of X on Y decreases as interVar increases
          interD <- sign(interC)

          # Alternative: Calculate how interaction modifies the main effect
          main_effect_row <- grep(paste0("^", xVar, "$"), rownames(dd))
          if (length(main_effect_row) > 0) {
            main_effect <- dd$Estimate[main_effect_row]
            # If interaction and main effect have opposite signs, interaction dampens the effect
            # If they have same signs, interaction amplifies the effect
            if (sign(main_effect) == sign(interC)) {
              interS <- 1  # Amplifies
            } else {
              interS <- -1  # Dampens/reverses
            }
          } else {
            interS <- NA
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
      })
    }
  }

  rownames(r_C) <- xVars
  colnames(r_C) <- yVars
  rownames(r_P) <- xVars
  colnames(r_P) <- yVars
  rownames(r_S) <- xVars
  colnames(r_S) <- yVars
  rownames(r_D) <- xVars
  colnames(r_D) <- yVars

  # Convert matrices to data frames and add variable names
  m_r_C <- r_C %>%
    as.data.frame() %>%
    rownames_to_column(var = "X") %>%
    pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "C")

  m_r_P <- r_P %>%
    as.data.frame() %>%
    rownames_to_column(var = "X") %>%
    pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "P")

  m_r_S <- r_S %>%
    as.data.frame() %>%
    rownames_to_column(var = "X") %>%
    pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "S")

  m_r_D <- r_D %>%
    as.data.frame() %>%
    rownames_to_column(var = "X") %>%
    pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "D")

  m_G <- left_join(m_r_C, m_r_P, by = c("X", "Y")) %>%
    left_join(m_r_S, by = c("X", "Y")) %>%
    left_join(m_r_D, by = c("X", "Y"))

  # Use the direction of the interaction coefficient for the sign
  m_G$sign <- factor(m_G$D, levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign[is.na(m_G$P) | m_G$P > 0.05] <- "ns"

  # Use custom function instead of gtools::stars.pval
  m_G$sig <- get_significance_stars(m_G$P)
  m_G$sig[m_G$sig == "." | m_G$sig == "+" | m_G$sig == " "] <- ""
  m_G$sigsign <- paste(as.character(m_G$sign), as.character(m_G$sig))

  m_G$sigsign <- factor(m_G$sigsign,
                        levels = c("+ ***", "+ **", "+ *", "ns ",
                                   "- *", "- **", "- ***"))

  ## FDR correction for r_P
  # Only adjust non-NA p-values
  p_values_for_adjustment <- m_G$P[!is.na(m_G$P)]
  if (length(p_values_for_adjustment) > 0) {
    adjusted_p_values <- p.adjust(p_values_for_adjustment, method = "fdr")
    m_G$P_FDR <- NA
    m_G$P_FDR[!is.na(m_G$P)] <- adjusted_p_values
  } else {
    m_G$P_FDR <- NA
  }

  m_G$sign_FDR <- factor(m_G$D, levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign_FDR[is.na(m_G$P_FDR) | m_G$P_FDR > 0.05] <- "ns"

  # Use custom function instead of gtools::stars.pval
  m_G$sig_FDR <- get_significance_stars(m_G$P_FDR)
  m_G$sig_FDR[m_G$sig_FDR == "." | m_G$sig_FDR == "+" | m_G$sig_FDR == " "] <- ""
  m_G$sig_FDR <- paste(m_G$sign_FDR, m_G$sig_FDR)
  m_G$sig_FDR <- factor(m_G$sig_FDR,
                        levels = c("+ ***", "+ **", "+ *", "ns ",
                                   "- *", "- **", "- ***"))

  ## add labels
  if (Relabel && requireNamespace("sjlabelled", quietly = TRUE)) {
    Data <- ReplaceMissingLabels(Data)

    # Get labels for X variables
    xlabels <- data.frame(
      Variable = unique(m_G$X),
      label = sjlabelled::get_label(Data[unique(m_G$X)], def.value = unique(m_G$X))
    )

    # Get labels for Y variables
    ylabels <- data.frame(
      Variable = unique(m_G$Y),
      label = sjlabelled::get_label(Data[unique(m_G$Y)], def.value = unique(m_G$Y))
    )

    m_G <- m_G %>%
      left_join(xlabels, by = c("X" = "Variable")) %>%
      rename(XLabel = label) %>%
      left_join(ylabels, by = c("Y" = "Variable")) %>%
      rename(YLabel = label)
  } else {
    m_G$XLabel <- m_G$X
    m_G$YLabel <- m_G$Y
  }

  m_G$XLabel <- factor(m_G$XLabel, ordered = FALSE,
                       levels = rev(unique(m_G$XLabel)))
  m_G$YLabel <- factor(m_G$YLabel, ordered = FALSE,
                       levels = unique(m_G$YLabel))

  # Create plot text with proper line breaks for ggplot
  m_G$PlotText <- paste("Y Var:", m_G$Y, "\nX Var:", m_G$X,
                        "\nInteraction Coef:", round(m_G$C, 4),
                        "\nP-Value:", round(m_G$P, 4), m_G$sigsign,
                        "\nFDR-corrected P:", round(m_G$P_FDR, 4), m_G$sig_FDR)

  # Create ggplot objects with better title
  p <- m_G %>%
    ggplot(aes(x = XLabel, y = YLabel, fill = sign)) +
    geom_tile(aes(text = PlotText), show.legend = TRUE) +
    scale_fill_manual(values = c("+ ***" = "darkgreen", "+ **" = "green",
                                 "+ *" = "lightgreen", "ns " = "white",
                                 "- *" = "lightcoral", "- **" = "red",
                                 "- ***" = "darkred"),
                      drop = FALSE, na.value = "grey90",
                      name = "Interaction\nDirection") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = paste("Interaction Effects with", interVar),
         subtitle = "Positive: X effect increases with interVar\nNegative: X effect decreases with interVar\n(No Multiple Comparison Correction)")

  p_FDR <- m_G %>%
    ggplot(aes(x = XLabel, y = YLabel, fill = sig_FDR)) +
    geom_tile(show.legend = TRUE) +
    scale_fill_manual(values = c("+ ***" = "darkgreen", "+ **" = "green",
                                 "+ *" = "lightgreen", "ns " = "white",
                                 "- *" = "lightcoral", "- **" = "red",
                                 "- ***" = "darkred"),
                      drop = FALSE, na.value = "grey90",
                      name = "Interaction\nDirection") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = paste("Interaction Effects with", interVar),
         subtitle = "Positive: X effect increases with interVar\nNegative: X effect decreases with interVar\n(FDR Correction)")

  # Create p-value tables
  pvaltable <- m_G %>%
    select(X, Y, P) %>%
    pivot_wider(names_from = X, values_from = P)

  pvaltable_FDR <- m_G %>%
    select(X, Y, P_FDR) %>%
    pivot_wider(names_from = X, values_from = P_FDR)

  # Create output list
  M <- list()
  M$C <- r_C
  M$S <- r_S
  M$P <- r_P
  M$D <- r_D
  M$plot <- p
  M$pvaltable <- pvaltable

  M_FDR <- list()
  M_FDR$C <- r_C
  M_FDR$S <- r_S
  M_FDR$P <- r_P  # Keep original P-values matrix
  M_FDR$D <- r_D
  M_FDR$P_adjusted <- m_G %>%
    select(X, Y, P_FDR) %>%
    pivot_wider(names_from = Y, values_from = P_FDR) %>%
    column_to_rownames("X") %>%
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
    raw_data = m_G  # Include the processed data frame for debugging
  ))
}
