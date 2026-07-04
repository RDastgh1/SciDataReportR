
## Todo: Add relabel option
#' Plot Numerical Interaction Effects Matrix
#'
#' This function calculates interaction effects between numerical variables and plots them as matrices.
#'
#' @param data The dataset containing the variables.
#' @param predictor_vars A character vector of the names of the x-axis numerical variables.
#' @param outcome_vars A character vector of the names of the y-axis numerical variables. Defaults to NULL.
#' @param xVarLabels A character vector of labels for the x-axis variables. Defaults to NULL.
#' @param yVarLabels A character vector of labels for the y-axis variables. Defaults to NULL.
#' @param interVar The interaction variable. Defaults to NULL.
#' @param covariates A character vector of the names of covariate variables. Defaults to NULL.
#' @return A list containing matrices, ggplot objects for visualizations, and tables of p-values.
#' @param fdr_scope Either `"matrix"` (default) or `"per_outcome"`, passed to
#'   [ApplyFDRCorrection()]. `"matrix"` corrects across all interaction
#'   p-values at once (historical behavior). `"per_outcome"` corrects
#'   separately within each y-axis variable (`outcome_vars`).
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param xVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param yVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @param covars \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @export
PlotNumInteractionEffectsMatrix <- function(data,
    predictor_vars,
    outcome_vars = NULL,
    xVarLabels = NULL,
    yVarLabels = NULL,
    interVar = NULL,
    covariates = NULL,
    Data = lifecycle::deprecated(),
    xVars = lifecycle::deprecated(),
    yVars = lifecycle::deprecated(),
    fdr_scope = c("matrix", "per_outcome"),
    covars = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotNumInteractionEffectsMatrix(Data)", "PlotNumInteractionEffectsMatrix(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(xVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotNumInteractionEffectsMatrix(xVars)", "PlotNumInteractionEffectsMatrix(predictor_vars)")
    predictor_vars <- xVars
  }
  if (!missing(predictor_vars)) xVars <- predictor_vars
  if (lifecycle::is_present(yVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotNumInteractionEffectsMatrix(yVars)", "PlotNumInteractionEffectsMatrix(outcome_vars)")
    outcome_vars <- yVars
  }
  yVars <- outcome_vars
  if (lifecycle::is_present(covars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotNumInteractionEffectsMatrix(covars)", "PlotNumInteractionEffectsMatrix(covariates)")
    covariates <- covars
  }
  covars <- covariates
  fdr_scope <- match.arg(fdr_scope)

  if (is.null(yVars)){
    yVars = xVars
  }
  if (is.null(xVarLabels)){
    xVarLabels = xVars
  }
  if (is.null(yVarLabels)){
    yVarLabels = yVars
  }

  if (is.null(interVar)) {
    stop("interVar must be provided: PlotNumInteractionEffectsMatrix() tests ",
         "predictor-by-interVar interaction effects, so there is nothing to ",
         "compute without it.")
  }

  # Update xVars and yVars to exclude covars
  xVarLabels <- xVarLabels[xVars %in% setdiff(xVars, covars)]
  yVarLabels <- yVarLabels[yVars %in% setdiff(yVars, covars)]
  xVars <- xVars[xVars %in% setdiff(xVars, covars)]
  yVars <- yVars[yVars %in% setdiff(yVars, covars)]

  r_P <- r_C <- r_S <- matrix(0, nrow = length(xVars), ncol = length(yVars))

  for (i in seq_len(length(xVars))) {
    for (j in seq_len(length(yVars))) {
      xVar <- xVars[i]
      yVar <- yVars[j]

      if (is.character(Data[[xVar]])){
        Data[[xVar]] <- as.factor(Data[[xVar]])
      }

      if (is.character(Data[[yVar]])){
        Data[[yVar]] <- as.factor(Data[[yVar]])
      }
      # Covariates are optional: with covars = NULL the covariate term is
      # dropped and the plain interaction model is fit. The formula is built
      # from column names so any number of covariates is supported.
      rhs <- paste(
        c(if (length(covars) > 0) paste0("`", covars, "`"),
          paste0("`", xVar, "` * `", interVar, "`")),
        collapse = " + "
      )
      fml <- stats::as.formula(paste0("`", yVar, "` ~ ", rhs))
      m <- tryCatch(lm(fml, data = Data), error = function(err) NULL)

      if (is.null(m)){
        interC <- NA
        interP <- NA
        interS <- NA
      } else {
        d <- summary(m)
        dd <- as.data.frame(d$coefficients)
        interC <- dd$Estimate %>% tail(1)
        interP <- dd$`Pr(>|t|)` %>% tail(1)
        s <- sign(dd$Estimate %>% tail(2))
        interS <- prod(s) * -1
      }

      if (xVar == yVar){
        interC <- NA
        interP <- NA
        interS <- NA
      }

      r_P[i, j] <- interP
      r_C[i, j] <- interC
      r_S[i, j] <- interS

    }
  }

  rownames(r_C) <- xVarLabels
  colnames(r_C) <- yVarLabels
  rownames(r_P) <- xVarLabels
  colnames(r_P) <- yVarLabels
  rownames(r_S) <- xVarLabels
  colnames(r_S) <- yVarLabels

  m_r_C <- r_C %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "C")

  m_r_P <- r_P %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "P")

  m_r_S <- r_S %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVarLabels), names_to = "Y", values_to = "S")

  m_G <- left_join(m_r_C, m_r_P) %>% left_join(m_r_S)
  m_G$sign <- factor(m_G$S, levels = c(-1, 0, 1), labels = c("-", "ns", "+")) # Switched
  m_G$sign[m_G$P > 0.05] <- "ns"
  m_G$sig <- gtools::stars.pval(m_G$P)
  m_G$sig[m_G$sig == "." | m_G$sig == "+"  | m_G$sig == " "] <- ""
  m_G$sig <- paste(m_G$sign, m_G$sig) %>% factor(levels = c("+ ***", "+ **", "+ *","ns ", "- *", "- **", "- ***"))

  ## FDR correction for r_P
  # Outcomes are the y-axis variables (outcome_vars / Y) for "per_outcome".
  m_G$P_FDR <- ApplyFDRCorrection(
    m_G$P,
    fdr_scope = fdr_scope,
    outcome_ids = m_G$Y
  )
  m_G$sign_FDR <- factor(m_G$S, levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign_FDR[m_G$P_FDR > 0.05] <- "ns"

  m_G$sig_FDR <- gtools::stars.pval(m_G$P_FDR)
  m_G$sig_FDR[m_G$sig_FDR == "." | m_G$sig_FDR == "+"  | m_G$sig_FDR == " "] <- ""
  m_G$sig_FDR <- paste(m_G$sign_FDR, m_G$sig_FDR) %>% factor(levels = c("+ ***", "+ **", "+ *","ns ", "- *", "- **", "- ***"))

  p <-  m_G %>% ggplot(aes(x = X, y = Y, fill = sig)) + geom_tile() +
    scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white", "lightblue2", "steelblue3", "blue")), drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "No Multiple Comparison Correction")

  p_FDR <-  m_G %>% ggplot(aes(x = X, y = Y, fill = sig_FDR)) + geom_tile() +
    scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white", "lightblue2", "steelblue3", "blue")), drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = interVar, subtitle = "FDR Correction")

  pvaltable_FDR <- m_G %>% data.frame() %>% select(X, Y, P_FDR) %>% pivot_wider(names_from = X, values_from = P_FDR)

  # Standardized p-value element aliases (old names kept)
  return(list(C = r_C, pvals = r_P, p = p, p_FDR = p_FDR, pvals_FDR = pvaltable_FDR,
              p_fdr = p_FDR, pvals_fdr = pvaltable_FDR))
}
