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

## Checking the match

A merge like this always produces a match, however far away it had to
reach, so the time gaps are what decide whether a match is usable. The
sign gives the direction: positive means the second record came after
the first.

In the example, three participants are matched within two weeks, while
the fourth has only a lab predating the visit by more than seven
months - that row should be dropped rather than analyzed as if the two
measurements were contemporaneous. Always inspect `time_differences`
before using the merged data.

## Examples

``` r
# Clinic visits with blood pressure
visits <- data.frame(
  id         = c("A", "B", "C", "D"),
  visit_date = as.Date(c("2024-03-01", "2024-05-20", "2024-02-10",
                         "2024-04-01")),
  sbp        = c(120, 135, 128, 142)
)

# Lab draws (multiple per participant, on different dates)
labs <- data.frame(
  id         = c("A", "A", "B", "B", "C", "D"),
  lab_date   = as.Date(c("2024-01-05", "2024-03-10", "2024-01-20",
                         "2024-06-01", "2024-02-15", "2023-08-15")),
  creatinine = c(0.9, 1.1, 0.8, 1.0, 1.2, 1.4)
)

ShowTable <- function(x, caption = NULL) {
  htmltools::browsable(htmltools::HTML(as.character(
    kableExtra::kable_styling(
      knitr::kable(x, format = "html", caption = caption),
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
  )))
}

# The two tables do not line up
ShowTable(visits, "Clinic visits")

Clinic visits
 
 id 
```
