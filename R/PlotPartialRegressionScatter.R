#' Partial Regression Plot
#'
#' Generate a partial regression plot for a specified independent and dependent variable
#' while adjusting for covariates. In addition to the figure, key parameters such as
#' the correlation method, whether relabeling was used, the covariates, R², p-value, and sample size
#' are returned.
#'
#' @param DataFrame The dataset to use.
#' @param IndepVar A string specifying the independent variable.
#' @param DepVar A string specifying the dependent variable.
#' @param Covariates A character vector of covariate names for adjustment. Defaults to NULL.
#' @param Relabel Logical indicating whether to use labelled names from the data (using sjlabelled::get_label). Defaults to TRUE.
#'
#' @return A list containing:
#' \item{plot}{A ggplot2 object representing the partial regression plot.}
#' \item{method}{The correlation method (as provided).}
#' \item{Relabel}{Logical; whether relabeling was applied.}
#' \item{Covariates}{The vector of covariates.}
#' \item{r2}{The R² of the partial regression model.}
#' \item{p_value}{The p-value for the independent variable coefficient.}
#' \item{n}{The sample size (number of complete cases).}
#' \item{equation}{The regression equation string.}
#'
#' @import ggplot2
#' @importFrom sjlabelled get_label
#' @export
PlotPartialRegressionScatter <- function(DataFrame, IndepVar, DepVar, Covariates = NULL, Relabel = TRUE) {

  # Subset complete cases for the relevant variables.
  allVars <- c(IndepVar, DepVar, Covariates)
  DataSub <- DataFrame[complete.cases(DataFrame[, allVars]), ]

  # If covariates are provided, residualize the independent and dependent variables.
  if (!is.null(Covariates)) {
    # Residualize the independent variable.
    f_indep <- as.formula(paste(IndepVar, "~", paste(Covariates, collapse = " + ")))
    model_indep <- lm(f_indep, data = DataSub)
    DataSub[[paste0(IndepVar, "_resid")]] <- resid(model_indep)

    # Residualize the dependent variable.
    f_dep <- as.formula(paste(DepVar, "~", paste(Covariates, collapse = " + ")))
    model_dep <- lm(f_dep, data = DataSub)
    DataSub[[paste0(DepVar, "_resid")]] <- resid(model_dep)
  } else {
    DataSub[[paste0(IndepVar, "_resid")]] <- DataSub[[IndepVar]]
    DataSub[[paste0(DepVar, "_resid")]] <- DataSub[[DepVar]]
  }

  # Fit the partial regression model: regress the dependent residuals on the independent residuals.
  partial_model <- lm(as.formula(paste0(paste0(DepVar, "_resid"), " ~ ", paste0(IndepVar, "_resid"))),
                      data = DataSub)
  summary_model <- summary(partial_model)
  coef_model <- coef(partial_model)
  r2 <- summary_model$r.squared
  n <- nrow(DataSub)
  p_val <- summary_model$coefficients[2, 4]

  # Build the regression equation string for the subtitle.
  eqString <- paste0("Residual ", DepVar, " = ", round(coef_model[1], 2),
                     " + ", round(coef_model[2], 2), " * ", tolower(paste0(IndepVar, "_resid")),
                     " | R² = ", round(r2, 2),
                     ", n = ", n,
                     ", p = ", formatC(p_val, format = "f", digits = 3))

  # Prepare the caption text to indicate the covariates.
  if (!is.null(Covariates)) {
    if (!Relabel) {
      capText <- paste("Adjusted for", paste(Covariates, collapse = ", "))
    } else {
      covarLabels <- sjlabelled::get_label(DataFrame[, Covariates], def.value = Covariates)
      capText <- paste("Adjusted for", paste(covarLabels, collapse = ", "))
    }
  } else {
    capText <- "Not adjusted for any covariates"
  }

  # Prepare axis labels.
  xLabel <- paste0(IndepVar, "_residual")
  yLabel <- paste("Residualized", DepVar)
  if (Relabel) {
    xLabel <- sjlabelled::get_label(DataFrame[[IndepVar]], def.value = IndepVar)
    yLabel <- sjlabelled::get_label(DataFrame[[DepVar]], def.value = DepVar)
  }

  # Create the ggplot object.
  plotObj <- ggplot(DataSub, aes_string(x = paste0(IndepVar, "_resid"), y = paste0(DepVar, "_resid"))) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "Partial Regression Plot",
         subtitle = eqString,
         x = xLabel,
         y = yLabel,
         caption = capText) +
    theme_minimal()

  # Return a list with the plot and key parameters.
  return(list(
    plot = plotObj,
    method = method,
    Relabel = Relabel,
    Covariates = Covariates,
    r2 = r2,
    p_value = p_val,
    n = n,
    equation = eqString
  ))
}
