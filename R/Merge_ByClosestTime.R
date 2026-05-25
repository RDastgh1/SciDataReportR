#' Merge Two Data Frames by Closest Time
#'
#' This function merges two data frames based on the closest time in the specified time columns.
#' It optionally merges using one or more matching variables (e.g., IDs).
#' The resulting merged data frame contains the closest time matches and time differences.
#'
#' @param DataFrame1 A data frame containing the first set of data.
#' @param DataFrame2 A data frame containing the second set of data.
#' @param TimeVar1 The name of the time variable in DataFrame1 (as a string).
#' @param TimeVar2 The name of the time variable in DataFrame2 (as a string).
#' @param MergeBy Optional. Character vector of variable(s) to merge by.
#'                Must exist in BOTH data frames and be in the same order.
#' @param is_date Logical. Indicates whether the time variables are dates (TRUE) or POSIXct (FALSE).
#'
#' @return A list with:
#' \item{merged_dataframe}{Data frame with closest time matches}
#' \item{time_differences}{Vector of time differences}
#'
#' @import dplyr
#' @import tidyr
#' @import lubridate
#' @export

Merge_ByClosestTime <- function(DataFrame1,
                                DataFrame2,
                                TimeVar1,
                                TimeVar2,
                                MergeBy = NULL,
                                is_date = FALSE) {



  # ---- Convert time variables
  if (is_date) {
    DataFrame1[[TimeVar1]] <- lubridate::as_date(DataFrame1[[TimeVar1]])
    DataFrame2[[TimeVar2]] <- lubridate::as_date(DataFrame2[[TimeVar2]])
  }

  # ---- Create working columns
  DataFrame1 <- DataFrame1 %>%
    dplyr::mutate(
      FunctionTime1 = .data[[TimeVar1]],
      DataFrame1Row = dplyr::row_number()
    )

  DataFrame2 <- DataFrame2 %>%
    dplyr::mutate(
      FunctionTime2 = .data[[TimeVar2]]
    )

  # ---- Handle merge keys
  if (!is.null(MergeBy)) {

    # Ensure variables exist
    missing_vars <- setdiff(MergeBy, intersect(names(DataFrame1), names(DataFrame2)))
    if (length(missing_vars) > 0) {
      stop(paste("MergeBy variable(s) not found in both data frames:", paste(missing_vars, collapse = ", ")))
    }

    # Coerce to character to avoid type mismatch hell
    DataFrame1 <- DataFrame1 %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(MergeBy), as.character))

    DataFrame2 <- DataFrame2 %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(MergeBy), as.character))

    # Controlled join ONLY on MergeBy
    merged_df <- dplyr::left_join(
      DataFrame1,
      DataFrame2,
      by = MergeBy,
      relationship = "many-to-many"
    )

  } else {

    # No merge key: full cross join (dangerous but intentional)
    merged_df <- tidyr::crossing(DataFrame1, DataFrame2)
  }

  # ---- Calculate time difference
  merged_df <- merged_df %>%
    dplyr::mutate(
      TimeDiff = as.numeric(FunctionTime2 - FunctionTime1)
    )

  # ---- Select closest match per row
  closest_matches <- merged_df %>%
    dplyr::group_by(DataFrame1Row) %>%
    dplyr::arrange(abs(TimeDiff), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # ---- Output prep
  time_differences <- closest_matches$TimeDiff

  merged_dataframe <- closest_matches %>%
    dplyr::select(
      -FunctionTime1,
      -FunctionTime2,
      -DataFrame1Row,
      -TimeDiff
    )

  return(list(
    merged_dataframe = merged_dataframe,
    time_differences = time_differences
  ))
}
