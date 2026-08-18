
## Todo: Add relabel options
## Todo: Add Covariate Option
#' Plot Categorical Interaction Effects Matrix
#'
#' This function calculates and visualizes the interaction effects between categorical variables.
#'
#' @param data The dataset containing the variables.
#' @param predictor_vars A character vector of the names of the x-axis categorical variables.
#' @param outcome_vars A character vector of the names of the y-axis categorical variables. Defaults to NULL,
#'        in which case it takes the same values as xVars.
#' @param xVarLabels A character vector of labels for the x-axis variables. Defaults to NULL,
#'        in which case it takes the same values as xVars.
#' @param yVarLabels A character vector of labels for the y-axis variables. Defaults to NULL,
#'        in which case it takes the same values as yVars.
#' @param interVar The name of the interaction variable.
#' @return A list containing matrices of interaction coefficients, p-values, ggplot objects for visualizations,
#'         and tables of FDR-corrected p-values.
#' @param fdr_scope Either `"matrix"` (default) or `"per_outcome"`, passed to
#'   [ApplyFDRCorrection()]. `"matrix"` corrects across all interaction
#'   p-values at once (historical behavior). `"per_outcome"` corrects
#'   separately within each y-axis variable (`outcome_vars`).
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param xVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param yVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels for readable axes
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' result <- PlotCatInteractionEffectsMatrix(
#'   Labelled,
#'   predictor_vars = c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
#'                      "Apolipoprotein_A_IV", "Apolipoprotein_A1",
#'                      "Apolipoprotein_A2", "Apolipoprotein_B",
#'                      "Apolipoprotein_CI", "Apolipoprotein_CIII",
#'                      "Apolipoprotein_D", "Apolipoprotein_E"),
#'   outcome_vars = c("age", "ACE_CD143_Angiotensin_Converti",
#'                    "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
#'                    "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
#'                    "Alpha_1_Microglobulin"),
#'   interVar = "Diagnosis"
#' )
#'
#' # Raw p-value interaction matrix
#' result$p
#'
#' # FDR-adjusted interaction matrix
#' result$p_FDR
#' @export
PlotCatInteractionEffectsMatrix <- function(data,
    predictor_vars,
    outcome_vars = NULL,
    xVarLabels = NULL,
    yVarLabels = NULL,
    interVar,
    Data = lifecycle::deprecated(),
    xVars = lifecycle::deprecated(),
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    yVars = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotCatInteractionEffectsMatrix(Data)", "PlotCatInteractionEffectsMatrix(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(xVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotCatInteractionEffectsMatrix(xVars)", "PlotCatInteractionEffectsMatrix(predictor_vars)")
    predictor_vars <- xVars
  }
  if (!missing(predictor_vars)) xVars <- predictor_vars
  if (lifecycle::is_present(yVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotCatInteractionEffectsMatrix(yVars)", "PlotCatInteractionEffectsMatrix(outcome_vars)")
    outcome_vars <- yVars
  }
  yVars <- outcome_vars
  fdr_scope <- match.arg(fdr_scope)

  if(is.null(yVars)){
    yVars = xVars
  }
  if(is.null(xVarLabels)){
    xVarLabels = xVars
  }
  if(is.null(yVarLabels)){
    yVarLabels = yVars
  }

  # Supplied names must exist. `Data[[xVar]]` returns NULL for a name that is not
  # a column, which used to fill the matrix with NA rather than report the typo.
  # Order and length are preserved because xVarLabels/yVarLabels pair positionally.
  ScidrValidateVariables(Data, xVars, "predictor_vars", unique_only = FALSE)
  ScidrValidateVariables(Data, yVars, "outcome_vars", unique_only = FALSE)
  ScidrValidateVariable(Data, interVar, "interVar")

  r_P <- r_C <- r_S <- matrix(0, nrow = length(xVars), ncol = length(yVars))

  for (i in seq_len(length(xVars))) {
    for (j in seq_len(length(yVars))) {
      xVar <- xVars[i]
      yVar <- yVars[j]

      if(is.character(Data[[xVar]])){
        Data[[xVar]] <- as.factor(Data[[xVar]])
      }

      if(is.character(Data[[yVar]])){
        Data[[yVar]] <- as.factor(Data[[yVar]])
      }
      hadError <- FALSE

      tryCatch({
        m <- lm(get(yVar) ~ get(xVar)*get(interVar), data = Data)
      }, error = function(err){
        hadError <<- TRUE
      }, interrupt = function(){
        hadError <<- FALSE
      })

      if(hadError){
        interC <- NA
        interP <- NA
        interS <- NA
      } else {
        m <- lm(get(yVar) ~ get(xVar)*get(interVar), data = Data)
        d <- summary(m)
        dd <- as.data.frame(d$coefficients)
        interC <- dd$Estimate %>% tail(1)
        interP <- dd$`Pr(>|t|)` %>% tail(1)
        interS <- sum(sign(dd$Estimate))
      }

      r_P[i,j] <- interP
      r_C[i,j] <- interC
      r_S[i,j] <- interS
    }
  }

  rownames(r_C) <- xVarLabels
  colnames(r_C) <- yVarLabels
  rownames(r_P) <- xVarLabels
  colnames(r_P) <- yVarLabels
  rownames(r_S) <- xVarLabels
  colnames(r_S) <- yVarLabels

  m_r_C <- r_C %>% as.data.frame %>% rownames_to_column(var = "X") %>% tidyr::pivot_longer(cols = dplyr::all_of(yVarLabels), names_to = "Y", values_to = "C")
  m_r_P <- r_P %>% as.data.frame %>% rownames_to_column(var = "X") %>% tidyr::pivot_longer(cols = dplyr::all_of(yVarLabels), names_to = "Y", values_to = "P")
  m_r_S <- r_S %>% as.data.frame %>% rownames_to_column(var = "X") %>% tidyr::pivot_longer(cols = dplyr::all_of(yVarLabels), names_to = "Y", values_to = "S")

  m_G <- left_join(m_r_C, m_r_P) %>% left_join(m_r_S)

  m_G$sign <- sign(m_G$S) %>% factor(levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign[m_G$P > 0.05] <- "ns"
  m_G$sig <- gtools::stars.pval(m_G$P)
  m_G$sig[m_G$sig %in% c(".", "+", " ")] <- ""
  m_G$sig <- paste(m_G$sign, m_G$sig)
  m_G$sig[is.na(m_G$P) | is.na(m_G$C)] <- "na"
  m_G$sig <- factor(m_G$sig, levels = c("+ ***", "+ **", "+ *", "ns ", "- *", "- **", "- ***", "na"))

  ## FDR correction for r_P
  # Outcomes are the y-axis variables (outcome_vars / Y) for "per_outcome".
  m_G$P_FDR <- ApplyFDRCorrection(
    m_G$P,
    fdr_scope = fdr_scope,
    outcome_ids = m_G$Y,
    predictor_ids = m_G$X
  )
  m_G$sign_FDR <- sign(m_G$S) %>% factor(levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign_FDR[m_G$P_FDR > 0.05] <- "ns"
  m_G$sig_FDR <- gtools::stars.pval(m_G$P_FDR)
  m_G$sig_FDR[m_G$sig_FDR %in% c(".", "+", " ")] <- ""
  m_G$sig_FDR <- paste(m_G$sign_FDR, m_G$sig_FDR)
  m_G$sig_FDR[is.na(m_G$P_FDR) | is.na(m_G$C)] <- "na"
  m_G$sig_FDR <- factor(m_G$sig_FDR, levels = levels(m_G$sig))

  effect_colors <- c(
    "+ ***" = "#0B1F5E", "+ **" = "#1769D2", "+ *" = "#B8D8E8",
    "ns " = "#FFFFFF", "- *" = "#F2B29A", "- **" = "#C2185B",
    "- ***" = "#841B37", "na" = "grey70")

  p <-  m_G %>%
    ggplot(aes(x = X, y = Y, fill = sig)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = effect_colors, drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "No Multiple Comparison Correction")

  p_FDR <-  m_G %>%
    ggplot(aes(x = X, y = Y, fill = sig_FDR)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = effect_colors, drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "FDR Correction")

  pvaltable_FDR <- m_G %>%
    data.frame() %>%
    dplyr::select(X, Y, P_FDR) %>%
    tidyr::pivot_wider(names_from = X, values_from = P_FDR)

  # Standardized p-value element aliases (old names kept)
  return(list(C = r_C, pvals = r_P, p = p, p_FDR = p_FDR, pvals_FDR = pvaltable_FDR,
              p_fdr = p_FDR, pvals_fdr = pvaltable_FDR))
}
