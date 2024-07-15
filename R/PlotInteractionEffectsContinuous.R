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

  # Check for missing values and replace with appropriate labels
  removediag <- FALSE
  if(is.null(yVars)){
    yVars = xVars
    removediag <- TRUE
  }

  # Combine variables into a single vector
  Variables <- c(interVar, xVars, yVars)

  # Replace missing labels
  Data <- ReplaceMissingLabels(Data)

  # If xVars is null, get all numeric variables
  if (is.null(xVars)) {
    xVars <- getNumVars(Data, Ordinal = F)
    if (Ordinal) {
      xVars <- getNumVars(Data, Ordinal = T)
    }
  }

  # If Ordinal is TRUE, convert ordinal variables to numeric and update Variables
  if (Ordinal) {
    Variables <- c(xVars, yVars)
    Data <- ConvertOrdinalToNumeric(Data, Variables)
    Data[Variables] <- lapply(Data[Variables], as.numeric)
  }

  # Update xVars and yVars to exclude covariates
  xVars <- xVars[xVars %in% setdiff(xVars, covars)]
  yVars <- yVars[yVars %in% setdiff(yVars, covars)]

  # Initialize matrices to store results
  r_P <- r_C <- r_S <- matrix(NA, nrow = length(xVars), ncol = length(yVars))

  # Loop through each combination of x and y variables
  for (i in seq_len(length(xVars))) {
    for (j in seq_len(length(yVars))) {
      xVar <- xVars[i]
      yVar <- yVars[j]

      # Convert character variables to factor
      if(is.character(Data[[xVar]])){
        Data[[xVar]] <- as.factor(Data[[xVar]])
      }

      if(is.character(Data[[yVar]])){
        Data[[yVar]] <- as.factor(Data[[yVar]])
      }

      # Initialize error flag
      hadError <- FALSE

      # Try to fit linear model, catching errors
      tryCatch({
        if(is.null(covars)){
          m <- lm(get(yVar) ~ get(xVar)*get(interVar), data = Data)
        } else {
          m <- lm(get(yVar) ~ get(covars) + get(xVar)*get(interVar), data = Data)
        }
      },
      error = function(err){
        hadError <<- TRUE
      },
      interrupt = function(){
        hadError <<- FALSE
      })

      # If an error occurred, set results to NA
      if(hadError){
        interC <- NA
        interP <- NA
        interS <- NA
      } else {
        # Summarize model and extract results
        d <- summary(m)
        dd <- as.data.frame(d$coefficients)
        interC <- dd$Estimate %>% tail(1)
        interP <- dd$`Pr(>|t|)` %>% tail(1)
        s <- sign(dd$Estimate %>% tail(2))
        interS <- prod(s) * -1
      }

      # If xVar is equal to yVar, set results to NA
      if(xVar == yVar){
        interC <- NA
        interP <- NA
        interS <- NA
      }

      # Store results in matrices
      r_P[i, j] <- interP
      r_C[i, j] <- interC
      r_S[i, j] <- interS
    }
  }

  rownames(r_C) <- xVars
  colnames(r_C) <- yVars
  rownames(r_P) <- xVars
  colnames(r_P) <- yVars
  rownames(r_S) <- xVars
  colnames(r_S) <- yVars
  # Convert matrices to data frames and add variable names
  m_r_C <- r_C %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "C")

  m_r_P<- r_P %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "P")

  m_r_S<- r_S %>% as.data.frame %>% rownames_to_column(var = "X") %>%  pivot_longer(cols = all_of(yVars), names_to = "Y", values_to = "S")

  m_G<- left_join(m_r_C, m_r_P) %>% left_join(m_r_S)
  m_G$sign <- factor(m_G$S, levels = c(-1, 0, 1), labels = c("-", "ns", "+")) # switched
  #m_G$sign <- -sign(m_G$C) %>% factor(levels = c(-1, 0, 1), labels = c("+", "ns", "-")) # switched
  m_G$sign[m_G$P>0.05]<-"ns"
  m_G$sig <- gtools::stars.pval(m_G$P)
  m_G$sig[ m_G$sig == "." | m_G$sig == "+"  | m_G$sig == " " ] <- ""
  m_G$sigsign<- paste(as.character(m_G$sign), as.character(m_G$sig))

  m_G$sigsign<- m_G$sigsign %>% factor(levels = c("+ ***", "+ **", "+ *","ns ",
                                                  "- *", "- **", "- ***"))


  ## FDR correction for r_P
  m_G$P_FDR <- p.adjust(m_G$P, method = "fdr")
  m_G$sign_FDR<- factor(m_G$S, levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  #m_G$sign_FDR <- sign(m_G$C) %>% factor(levels = c(-1, 0, 1), labels = c("-", "ns", "+"))
  m_G$sign_FDR[m_G$P_FDR>0.05]<-"ns"

  m_G$sig_FDR <- gtools::stars.pval(m_G$P_FDR)
  m_G$sig_FDR[ m_G$sig_FDR == "." | m_G$sig_FDR == "+"  | m_G$sig_FDR == " " ] <- ""
  m_G$sig_FDR<- paste(m_G$sign_FDR, m_G$sig_FDR) %>% factor(levels = c("+ ***", "+ **", "+ *","ns ",
                                                                       "- *", "- **", "- ***"))


  ## add labels

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)
    xlabels <- sjlabelled::get_label(Data[m_G$X],
                                     def.value = colnames(Data[m_G$X])) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")
    ylabels <- sjlabelled::get_label(Data[m_G$Y],
                                     def.value = colnames(Data[m_G$Y])) %>% as.data.frame() %>%
      rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")
    m_G$XLabel <- xlabels$label
    m_G$YLabel <- ylabels$label
  }
  else {
    m_G$XLabel <- m_G$X
    m_G$YLabel <- m_G$Y
  }

  m_G$XLabel <- factor(m_G$XLabel, ordered = FALSE,
                       levels = rev(unique(m_G$XLabel)))
  m_G$YLabel <- factor(m_G$YLabel, ordered = FALSE,
                       levels = unique(m_G$YLabel))

  m_G$PlotText <- paste("YVar", m_G$Y, "</br> XVAR: ",
                        m_G$X, "</br> P-Value: ",m_G$P, m_G$sigsign,
                        "</br> FDR-corrected P: ", m_G$P_FDR, m_G$sig_FDR)



                   # Create ggplot objects
                   p <- m_G %>% ggplot(aes(x = XLabel, y = YLabel, fill = sigsign)) +
                   geom_tile(aes(text = PlotText), show.legend = T) +
                   scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white",
                                                    "lightblue2", "steelblue3", "blue")), drop = F) +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         legend.title=element_blank()) +
                   labs(title = interVar, subtitle = "No Multiple Comparison Correction")

                 p_FDR <- m_G %>% ggplot(aes(x = XLabel, y = YLabel, fill = sig_FDR)) +
                   geom_tile(show.legend = T) +
                   scale_fill_manual(values = rev(c("red4", "firebrick3", "pink2", "white",
                                                    "lightblue2", "steelblue3", "blue")), drop = F) +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         legend.title=element_blank()) +
                   labs(title = interVar, subtitle = "FDR Correction")

                 # Create p-value tables
                 pvaltable <- m_G %>% data.frame() %>% select(X, Y, P) %>% pivot_wider(names_from = X, values_from = P)

                 pvaltable_FDR <- m_G %>% data.frame() %>% select(X, Y, P_FDR) %>% pivot_wider(names_from = X, values_from = P_FDR)

                 # Create output list
                 M <- list()
                 M$C <- m_r_C
                 M$S <- m_r_S
                 M$P <- m_r_P
                 padjusted <- p.adjust(M$P$P, method = "fdr")
                 M_FDR <- M
                 M_FDR$P$P <- padjusted
                 M$plot <- p
                 M_FDR$plot <- p_FDR
                 M_FDR$pvaltable <- pvaltable_FDR
                 M$pvaltable <- pvaltable

                 # Return output list
                 return(list(Unadjusted = M, FDRCorrected = M_FDR, Relabel = Relabel, Covariates = covars, interVar = interVar))
                 }
