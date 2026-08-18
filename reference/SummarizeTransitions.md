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

- `Plots`: figures for the two summaries, described below

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

## Figures

`condition_summary` is a single row of six counts, which is exactly the
shape a table reads worst: the numbers are related to each other, and
the relationship is the point. `Plots` renders them so the relationship
is visible:

- `Plots$ConditionCascade` shows the counts as a cascade - the whole
  cohort, the part of it ever affected, the transitions that occurred
  within that part, and how many of those persisted - with each bar
  labelled by its count and its share of the cohort.

- `Plots$TransitionPatterns` classifies every participant into one
  mutually exclusive longitudinal pattern (never positive, positive
  throughout, developed only, resolved only, developed and resolved) and
  plots the counts. The indicator columns in `participant_summary`
  overlap, so this is the view that shows what the cohort is actually
  made of.

The two agree by construction: the pattern counts sum to
`n_participants`, and everything except "never positive" sums to
`n_ever_positive`.

For the per-participant trajectories behind these counts, use
[`PlotSwimmerTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSwimmerTransitions.md),
which prepares its data the same way.

## See also

[`PlotSwimmerTransitions()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotSwimmerTransitions.md)
for the participant-level swimmer plot.

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

transitions <- SummarizeTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary,
  date_var = VisitDate,
  x_axis_type = "time_from_baseline",
  time_from_baseline_unit = "months"
)

transitions$condition_summary
#> # A tibble: 1 × 6
#>   n_participants n_ever_positive n_developed_condition n_resolved_condition
#>            <int>           <int>                 <int>                <int>
#> 1              4               3                     1                    1
#> # ℹ 2 more variables: n_sustained_after_development <int>,
#> #   n_sustained_after_resolution <int>

# \donttest{
# A larger cohort, where a positive visit tends to be followed by another
set.seed(9)
n_participants <- 60

df_Longitudinal <- do.call(rbind, lapply(seq_len(n_participants), function(i) {
  n_visits <- sample(3:5, 1)
  status <- numeric(n_visits)
  status[1] <- stats::rbinom(1, 1, 0.3)
  for (j in seq_len(n_visits)[-1]) {
    status[j] <- stats::rbinom(1, 1, if (status[j - 1] == 1) 0.8 else 0.25)
  }
  data.frame(
    ParticipantID = sprintf("P%02d", i),
    VisitOrder = seq_len(n_visits),
    MetSBinary = status
  )
}))

transitions <- SummarizeTransitions(
  data = df_Longitudinal,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary
)

# Six counts in one row
htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(transitions$condition_summary, full_width = TRUE)
)))


 n_participants 
```
