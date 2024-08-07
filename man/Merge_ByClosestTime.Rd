% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Merge_ByClosestTime.R
\name{Merge_ByClosestTime}
\alias{Merge_ByClosestTime}
\title{Merge Two Data Frames by Closest Time}
\usage{
Merge_ByClosestTime(
  DataFrame1,
  DataFrame2,
  TimeVar1,
  TimeVar2,
  IDVar1 = NULL,
  IDVar2 = NULL,
  is_date = FALSE
)
}
\arguments{
\item{DataFrame1}{A data frame containing the first set of data.}

\item{DataFrame2}{A data frame containing the second set of data.}

\item{TimeVar1}{The name of the time variable in DataFrame1 (as a string).}

\item{TimeVar2}{The name of the time variable in DataFrame2 (as a string).}

\item{IDVar1}{Optional. The name of the ID variable in DataFrame1 (as a string). If not provided, no ID-based join is performed.}

\item{IDVar2}{Optional. The name of the ID variable in DataFrame2 (as a string). If not provided, no ID-based join is performed.}

\item{is_date}{Logical. Indicates whether the time variables are dates (TRUE) or POSIXct timestamps (FALSE). Defaults to FALSE.}
}
\value{
A list containing two elements:
\item{merged_dataframe}{A data frame with the closest time matches from both data frames and the time differences.}
\item{time_differences}{A numeric vector of time differences for each closest match.}
}
\description{
This function merges two data frames based on the closest time in the specified time columns.
It can also merge data frames using optional ID columns. The resulting merged data frame
contains the closest time matches and the time differences.
}
