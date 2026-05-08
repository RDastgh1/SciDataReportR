# Merge Two Data Frames by Closest Time

This function merges two data frames based on the closest time in the
specified time columns. It optionally merges using one or more matching
variables (e.g., IDs). The resulting merged data frame contains the
closest time matches and time differences.

## Usage

``` r
Merge_ByClosestTime(
  DataFrame1,
  DataFrame2,
  TimeVar1,
  TimeVar2,
  MergeBy = NULL,
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

- MergeBy:

  Optional. Character vector of variable(s) to merge by. Must exist in
  BOTH data frames and be in the same order.

- is_date:

  Logical. Indicates whether the time variables are dates (TRUE) or
  POSIXct (FALSE).

## Value

A list with:

- merged_dataframe:

  Data frame with closest time matches

- time_differences:

  Vector of time differences
