# Plot swimmer-style transitions for a binary condition over repeated visits

Create a swimmer-style longitudinal plot for a binary condition measured
across repeated visits. The plot shows condition status at each visit,
highlights transition points where the condition develops or resolves,
and supports multiple participant ordering strategies for exploratory
quality control, longitudinal debugging, or figure creation.

## Usage

``` r
PlotSwimmerTransitions(
  data,
  id_var,
  time_var,
  status_var,
  date_var = NULL,
  participant_subset = NULL,
  max_participants = NULL,
  order_participants_by = c("first_positive", "first_transition", "ever_positive",
    "ever_positive_then_burden", "input_order", "n_visits", "n_positive", "pct_positive"),
  x_axis_type = c("visit", "date"),
  show_transition_points = TRUE,
  show_lines = TRUE,
  show_y_axis_labels = FALSE,
  make_interactive = FALSE,
  plot_title = NULL,
  x_label = NULL,
  y_label = NULL,
  return_data = FALSE
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
  Accepted encodings include numeric 0/1, logical TRUE/FALSE, factor
  values such as `"Yes"` and `"No"`, and character values such as `"1"`
  and `"0"`.

- date_var:

  Optional unquoted visit date column. This is required when
  `x_axis_type = "date"`.

- participant_subset:

  Optional vector of participant IDs to include.

- max_participants:

  Optional maximum number of participants to display after ordering is
  applied.

- order_participants_by:

  Character string controlling participant order in the plot. Options
  are `"first_positive"`, `"first_transition"`, `"ever_positive"`,
  `"ever_positive_then_burden"`, `"input_order"`, `"n_visits"`,
  `"n_positive"`, and `"pct_positive"`.

- x_axis_type:

  Character string indicating whether the x-axis should use aligned
  visit number (`"visit"`) or actual calendar date (`"date"`).

- show_transition_points:

  Logical. If `TRUE`, transition visits are highlighted.

- show_lines:

  Logical. If `TRUE`, a swimmer line is drawn across visits within each
  participant.

- show_y_axis_labels:

  Logical. If `TRUE`, show participant labels on the y-axis. Defaults to
  `FALSE`.

- make_interactive:

  Logical. If `TRUE`, return an interactive plotly object. Defaults to
  `FALSE`.

- plot_title:

  Optional custom plot title.

- x_label:

  Optional x-axis label.

- y_label:

  Optional y-axis label. Defaults to `NULL`.

- return_data:

  Logical. If `TRUE`, return both the plot and the processed data.

## Value

A `ggplot` object by default.

If `make_interactive = TRUE`, returns a `plotly` object.

If `return_data = TRUE`, returns a list with:

- `plot`: the `ggplot` or `plotly` object

- `plot_data`: the processed visit-level plotting data

- `participant_summary`: the participant-level summary table

## Details

Visits are sorted within participant by `date_var` if provided,
otherwise by `time_var`.

Transition rules are defined as:

- `0 -> 1` = `"Developed"`

- `1 -> 0` = `"Resolved"`

Missing status values remain `NA` and are not forced to 0.

Participant-level summaries are created internally, including:

- whether the participant was ever positive

- first positive visit

- first transition visit

- number of observed visits

- number and proportion of positive visits

- whether positivity is sustained after first onset

- whether resolution is sustained after first resolution

Sustained after development means that all non-missing observations
after the first `0 -> 1` transition remain positive.

Sustained after resolution means that all non-missing observations after
the first `1 -> 0` transition remain negative.

When `x_axis_type = "visit"`, participants are aligned by visit order.
When `x_axis_type = "date"`, visits are shown at actual calendar dates
and do not need to align across participants.

When `make_interactive = TRUE`, the function returns a `plotly` object
and uses the internally prepared tooltip text for hover labels.

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
    "No", "Yes", "Yes", "Yes"
  )
)

PlotSwimmerTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary
)


PlotSwimmerTransitions(
  data = toy_df,
  id_var = ParticipantID,
  time_var = VisitOrder,
  status_var = MetSBinary,
  date_var = VisitDate,
  x_axis_type = "date"
)
```
