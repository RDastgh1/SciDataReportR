% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createBinaryMapping.R
\name{createBinaryMapping}
\alias{createBinaryMapping}
\title{Create a Mapping Table for Binary Variables}
\usage{
createBinaryMapping(Data, CatVars)
}
\arguments{
\item{Data}{A dataframe containing the data to analyze.}

\item{CatVars}{A character vector of column names representing categorical variables in the dataframe.}
}
\value{
A dataframe with columns: \code{VariableName}, \code{Label}, and \code{PositiveValue} for binary variables found in \code{CatVars}.
}
\description{
This function identifies binary variables within a given dataframe and generates a mapping table.
The mapping table includes the variable name, its associated label, and the positive value.
}
