# Summarize participant transitions for a binary longitudinal condition

Create participant-level and condition-level summary tables for a binary
condition observed across repeated visits. This function uses the same
transition logic and data preparation workflow as
[`PlotSwimmerTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSwimmerTransitions.md)
so that plotting and summary outputs remain aligned.

## Usage

``` r
SummarizeTransitions(
  data,
  id_var,
  time_var,
  status_var,
  date_var = NULL,
  participant_subset = NULL,
  max_participants = NULL,
  order_participants_by = c("first_positive", "first_transition", "ever_positive",
    "ever_positive_then_burden", "input_order", "n_visits", "n_positive", "pct_positive"),
  x_axis_type = c("visit", "date", "time_from_baseline"),
  time_from_baseline_unit = c("days", "months", "years")
)
```

## Arguments

- data:

  A data frame containing repeated observations per participant.

- id_var:

  Unquoted column name identifying the participant.

- time_var:

  Unquoted column name representing visit order, visit number, or time
  index.

- status_var:

  Unquoted column name representing the binary condition status.

- date_var:

  Optional unquoted visit date column. This is required when
  `x_axis_type = "date"` or `x_axis_type = "time_from_baseline"`.

- participant_subset:

  Optional vector of participant IDs to include.

- max_participants:

  Optional maximum number of participants to retain after ordering is
  applied.

- order_participants_by:

  Character string controlling participant order. Options are
  `"first_positive"`, `"first_transition"`, `"ever_positive"`,
  `"ever_positive_then_burden"`, `"input_order"`, `"n_visits"`,
  `"n_positive"`, and `"pct_positive"`.

- x_axis_type:

  Character string indicating whether longitudinal ordering should
  follow aligned visit number (`"visit"`), actual date (`"date"`), or
  elapsed time from baseline (`"time_from_baseline"`).

- time_from_baseline_unit:

  Character string specifying the unit for
  `x_axis_type = "time_from_baseline"`. Options are `"days"`,
  `"months"`, and `"years"`.

## Value

A list with:

- `participant_summary`: participant-level summary table

- `condition_summary`: one-row tibble with overall counts

## Details

Transition rules are:

- `0 -> 1` = developed condition

- `1 -> 0` = resolved condition

Missing values remain missing and are not recoded to 0.

The returned condition-level summary includes:

- number of participants

- number ever positive

- number who developed the condition

- number who resolved the condition

- number sustained after development

- number sustained after resolution

## Examples

``` r
toy_df <- tibble::tibble(
  ParticipantID = rep(paste0("P", 1:4), each = 4),
  VisitOrder = rep(1:4, times = 4),
  VisitDate = rep(seq.Date(as.Date("2024-01-01"), by = "month", length.out = 4), times = 4),
  MetSBinary = c(
    0, 0, 1, 1,
    1, 1, 0, 0,
    0, 0, 0, 0,
    TRUE, TRUE, TRUE, TRUE
  )
)

SummarizeTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary,
  date_var = VisitDate,
  x_axis_type = "time_from_baseline",
  time_from_baseline_unit = "months"
)
#> $participant_summary
#> # A tibble: 4 × 23
#>   .plot_id n_rows n_visits n_positive pct_positive ever_positive
#>   <chr>     <int>    <int>      <int>        <dbl> <lgl>        
#> 1 P2            4        4          2          0.5 TRUE         
#> 2 P4            4        4          4          1   TRUE         
#> 3 P1            4        4          2          0.5 TRUE         
#> 4 P3            4        4          0          0   FALSE        
#> # ℹ 17 more variables: developed_condition <lgl>, resolved_condition <lgl>,
#> #   first_positive_time <dbl>, first_positive_date <date>,
#> #   first_transition_time <dbl>, first_transition_date <date>,
#> #   first_developed_time <dbl>, first_resolved_time <dbl>,
#> #   first_developed_date <date>, first_resolved_date <date>,
#> #   baseline_date <date>, input_order <int>, sustained_after_development <lgl>,
#> #   sustained_after_resolution <lgl>, .order_group <dbl>, .order_value <dbl>, …
#> 
#> $condition_summary
#> # A tibble: 1 × 6
#>   n_participants n_ever_positive n_developed_condition n_resolved_condition
#>            <int>           <int>                 <int>                <int>
#> 1              4               3                     1                    1
#> # ℹ 2 more variables: n_sustained_after_development <int>,
#> #   n_sustained_after_resolution <int>
#> 
```
