% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotMiningMatrix.R
\name{PlotMiningMatrix}
\alias{PlotMiningMatrix}
\title{PlotMiningMatrix}
\usage{
PlotMiningMatrix(
  Data,
  OutcomeVars,
  PredictorVars,
  Covariates = NULL,
  Relabel = TRUE,
  Parametric = TRUE
)
}
\arguments{
\item{Data}{A data frame containing the dataset to analyze.}

\item{OutcomeVars}{A vector of outcome variables to be analyzed.}

\item{PredictorVars}{A vector of predictor variables to be analyzed.}

\item{Covariates}{An optional vector of covariates to adjust the analysis (default is NULL).}

\item{Relabel}{Logical flag indicating whether to relabel the variables in the output (default is TRUE).}

\item{Parametric}{Logical flag indicating whether to use parametric methods (default is TRUE).
If FALSE, non-parametric methods will be used.}
}
\value{
A list containing the following elements:
\item{Unadjusted}{A list with the unadjusted p-value table and corresponding plot.}
\item{FDRCorrected}{A list with the FDR-adjusted p-value table and corresponding plot.}
\item{method}{The method used for correlation ("pearson" for parametric, "spearman" for non-parametric).}
\item{Relabel}{The value of the Relabel parameter.}
\item{Covariates}{The covariates used in the analysis, if any.}
}
\description{
This function generates a matrix of statistical relationships between specified outcome and predictor variables
in a dataset. It includes visualizations for correlations, ANOVA results, and FDR corrections.
}
