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
  keys = NULL,
  is_date = FALSE,
  MergeBy = lifecycle::deprecated()
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

- keys:

  Optional. Character vector of variable(s) to merge by. Must exist in
  BOTH data frames and be in the same order.

- is_date:

  Logical. Indicates whether the time variables are dates (TRUE) or
  POSIXct (FALSE).

- MergeBy:

  **Deprecated** (since 19.15.0). Use `keys` instead.

## Value

A list with:

- merged_dataframe:

  Data frame with closest time matches

- time_differences:

  Vector of time differences

## Examples

``` r
# Clinic visits with blood pressure
visits <- data.frame(
  id         = c("A", "B", "C"),
  visit_date = as.Date(c("2024-03-01", "2024-05-20", "2024-02-10")),
  sbp        = c(120, 135, 128)
)

# Lab draws (multiple per participant, on different dates)
labs <- data.frame(
  id         = c("A", "A", "B", "B", "C"),
  lab_date   = as.Date(c("2024-01-05", "2024-03-10", "2024-01-20",
                         "2024-06-01", "2024-02-15")),
  creatinine = c(0.9, 1.1, 0.8, 1.0, 1.2)
)

# For each visit, attach the lab drawn closest in time (within participant)
res <- Merge_ByClosestTime(
  visits, labs,
  TimeVar1 = "visit_date",
  TimeVar2 = "lab_date",
  keys = "id",
  is_date = TRUE
)

res$merged_dataframe
#> # A tibble: 3 × 5
#>   id    visit_date   sbp lab_date   creatinine
#>   <chr> <date>     <dbl> <date>          <dbl>
#> 1 A     2024-03-01   120 2024-03-10        1.1
#> 2 B     2024-05-20   135 2024-06-01        1  
#> 3 C     2024-02-10   128 2024-02-15        1.2
res$time_differences
#> [1]  9 12  5
```
