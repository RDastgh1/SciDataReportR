#' Merge Two Data Frames by Closest Time
#'
#' This function merges two data frames based on the closest time in the specified time columns.
#' It can also merge data frames using optional ID columns. The resulting merged data frame
#' contains the closest time matches and the time differences.
#'
#' @param DataFrame1 A data frame containing the first set of data.
#' @param DataFrame2 A data frame containing the second set of data.
#' @param TimeVar1 The name of the time variable in DataFrame1 (as a string).
#' @param TimeVar2 The name of the time variable in DataFrame2 (as a string).
#' @param IDVar1 Optional. The name of the ID variable in DataFrame1 (as a string). If not provided, no ID-based join is performed.
#' @param IDVar2 Optional. The name of the ID variable in DataFrame2 (as a string). If not provided, no ID-based join is performed.
#' @param is_date Logical. Indicates whether the time variables are dates (TRUE) or POSIXct timestamps (FALSE). Defaults to FALSE.
#'
#' @return A list containing two elements:
#' \item{merged_dataframe}{A data frame with the closest time matches from both data frames and the time differences.}
#' \item{time_differences}{A numeric vector of time differences for each closest match.}
#'
#' @import dplyr
#' @importFrom lubridate as_date
#' @export
Merge_ByClosestTime <- function(DataFrame1, DataFrame2, TimeVar1, TimeVar2, IDVar1 = NULL, IDVar2 = NULL, is_date = FALSE) {
  library(dplyr)
  library(lubridate)

  # Convert time variables to Date or POSIXct
  if (is_date) {
    DataFrame1[[TimeVar1]] <- as_date(DataFrame1[[TimeVar1]])
    DataFrame2[[TimeVar2]] <- as_date(DataFrame2[[TimeVar2]])
  }

  # Add time columns for merge and row numbers for later reference
  DataFrame1$FunctionTime1 <- DataFrame1[[TimeVar1]]
  DataFrame2$FunctionTime2 <- DataFrame2[[TimeVar2]]
  DataFrame1$DataFrame1Row <- 1:nrow(DataFrame1)

  # If ID variables are provided, join based on ID
  if (!is.null(IDVar1) && !is.null(IDVar2)) {
    DataFrame1$FunctionID <- as.character(DataFrame1[[IDVar1]])
    DataFrame2$FunctionID <- as.character(DataFrame2[[IDVar2]])
    merged_df <- left_join(DataFrame1, DataFrame2, by = "FunctionID", suffix = c(".df1", ".df2"), relationship = "many-to-many")
  } else {
    merged_df <- left_join(DataFrame1, DataFrame2, suffix = c(".df1", ".df2"), relationship = "many-to-many")
  }

  # Calculate the time difference
  merged_df$TimeDiff <- merged_df$FunctionTime2 - merged_df$FunctionTime1

  # Find the closest time match for each row in DataFrame1
  closest_matches <- merged_df %>%
    group_by(DataFrame1Row) %>%
    arrange(abs(TimeDiff)) %>%
    slice(1) %>%
    ungroup()

  # Prepare the results
  time_differences <- closest_matches$TimeDiff
  closest_matches <- closest_matches %>%
    select(-TimeDiff)

  # Remove unnecessary columns
  merged_dataframe <- closest_matches %>%
    select(-FunctionTime1, -FunctionTime2, -DataFrame1Row) %>% select(-any_of("FunctionID"))

  return(list(merged_dataframe = merged_dataframe, time_differences = time_differences))
}
