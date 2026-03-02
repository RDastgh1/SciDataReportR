# Merge Two Data Frames by Closest Time

This function merges two data frames based on the closest time in the
specified time columns. It can also merge data frames using optional ID
columns. The resulting merged data frame contains the closest time
matches and the time differences.

## Usage

``` r
Merge_ByClosestTime(
  DataFrame1,
  DataFrame2,
  TimeVar1,
  TimeVar2,
  IDVar1 = NULL,
  IDVar2 = NULL,
  is_date = FALSE
)
```

## Arguments

- DataFrame1:

  A data frame containing the first set of data.

- DataFrame2:

  A data frame containing the second set of data.

- TimeVar1:

  The name of the time variable in DataFrame1 (as a string).

- TimeVar2:

  The name of the time variable in DataFrame2 (as a string).

- IDVar1:

  Optional. The name of the ID variable in DataFrame1 (as a string). If
  not provided, no ID-based join is performed.

- IDVar2:

  Optional. The name of the ID variable in DataFrame2 (as a string). If
  not provided, no ID-based join is performed.

- is_date:

  Logical. Indicates whether the time variables are dates (TRUE) or
  POSIXct timestamps (FALSE). Defaults to FALSE.

## Value

A list containing two elements:

- merged_dataframe:

  A data frame with the closest time matches from both data frames and
  the time differences.

- time_differences:

  A numeric vector of time differences for each closest match.
