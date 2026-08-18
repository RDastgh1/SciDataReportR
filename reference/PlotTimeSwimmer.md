# Plot longitudinal swimmer timelines

Create swimmer-style longitudinal timeline plots from long-format visit
data. Each row in the input data represents a subject visit or
observation. The function can visualize longitudinal state trajectories,
visit timing, and events over time.

## Usage

``` r
PlotTimeSwimmer(
  data,
  id_var,
  Time,
  State = NULL,
  Event = NULL,
  EventType = NULL,
  TimeScale = c("from_first", "observed", "from_event"),
  EventReference = NULL,
  StateInterval = c("forward", "point"),
  Format = c("state_path", "visit_points", "event_rug", "minimal"),
  SortBy = c("duration", "last_time", "first_time", "state", "id"),
  TimeUnit = c("auto", "days", "weeks", "months", "years", "visits"),
  Relabel = TRUE,
  codebook = NULL,
  LineWidth = 5,
  PointSize = 2.5,
  Alpha = 0.9,
  BaseSize = 13,
  Data = lifecycle::deprecated(),
  ID = lifecycle::deprecated(),
  Codebook = lifecycle::deprecated()
)
```

## Arguments

- data:

  A data frame in long format with one row per visit or observation.

- id_var:

  Character string naming the participant ID column.

- Time:

  Character string naming the time variable.

- State:

  Optional character string naming a state/group variable used for
  coloring timelines or visit points.

- Event:

  Optional character string naming a logical or binary event indicator
  variable.

- EventType:

  Optional character string naming an event category variable. Used for
  event point shapes/colors.

- TimeScale:

  Character. One of `"from_first"`, `"observed"`, or `"from_event"`.

  `"observed"` uses raw observed time values.

  `"from_first"` normalizes each subject relative to their first
  observed timepoint.

  `"from_event"` normalizes each subject relative to the first
  occurrence of `EventReference`.

- EventReference:

  Optional event value used when `TimeScale = "from_event"`.

- StateInterval:

  Character. One of `"forward"` or `"point"`.

  `"forward"` extends the current state forward until the next visit.

  `"point"` only colors visit points without extending intervals.

- Format:

  Character. One of `"state_path"`, `"visit_points"`, `"event_rug"`, or
  `"minimal"`.

- SortBy:

  Character. One of `"duration"`, `"last_time"`, `"first_time"`,
  `"state"`, or `"id"`.

- TimeUnit:

  Character. One of `"auto"`, `"days"`, `"weeks"`, `"months"`,
  `"years"`, or `"visits"`.

- Relabel:

  Logical. If `TRUE`, use labels from `Codebook` or variable attributes
  when available.

- codebook:

  Optional codebook data frame with columns `Variable` and `Label`.

- LineWidth:

  Numeric line width for swimmer segments. Default is `5`.

- PointSize:

  Numeric point size for visit/event points. Default is `2.5`.

- Alpha:

  Numeric alpha transparency for swimmer segments. Default is `0.9`.

- BaseSize:

  Base font size for the plot theme. Default is `13`.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- ID:

  **Deprecated** (since 19.15.0). Use `id_var` instead.

- Codebook:

  **Deprecated** (since 19.15.0). Use `codebook` instead.

## Value

A `ggplot` object.

## Details

Time can be represented using dates or numeric visit values. Timelines
can optionally be normalized relative to each participant's first visit
or the first occurrence of a specified event.

States can either be displayed as continuous intervals extending forward
until the next visit (`StateInterval = "forward"`) or shown only as
visit points (`StateInterval = "point"`).

## What the examples vary

State paths extend each state to the next visit, and ordering by
duration puts the longest observed timelines at the top. Visit points
instead retain only the observed visits, and `last_time` ranks
participants by their most recent observed follow-up.

Event rugs center time on each participant's first flare, so the input
includes only participants who experienced the selected event; ordering
by state then groups participants by their initial clinical state.

## Examples

``` r
df_swimmer <- tibble::tibble(
  ID = c("P01", "P01", "P01", "P02", "P02", "P03", "P03", "P03", "P04", "P04"),
  Day = c(0, 90, 365, 0, 270, 0, 180, 540, 0, 120),
  State = c("Stable", "Flare", "Recovered", "Stable", "Stable",
            "Flare", "Flare", "Recovered", "Stable", "Withdrawn"),
  Event = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE),
  EventType = c(NA, "Flare", NA, NA, NA, "Flare", NA, NA, NA, "Withdrawal")
)

# State paths, ordered by duration
PlotTimeSwimmer(
  data = df_swimmer,
  id_var = "ID",
  Time = "Day",
  State = "State",
  Event = "Event",
  EventType = "EventType",
  Format = "state_path",
  SortBy = "duration",
  TimeUnit = "months"
)
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).


# Visit points, ordered by last observed follow-up
PlotTimeSwimmer(
  data = df_swimmer,
  id_var = "ID",
  Time = "Day",
  State = "State",
  Format = "visit_points",
  SortBy = "last_time",
  TimeUnit = "months"
)


# Event rugs, centered on each participant's first flare
PlotTimeSwimmer(
  data = subset(df_swimmer, ID %in% c("P01", "P03")),
  id_var = "ID",
  Time = "Day",
  State = "State",
  Event = "Event",
  EventType = "EventType",
  TimeScale = "from_event",
  EventReference = TRUE,
  Format = "event_rug",
  SortBy = "state",
  TimeUnit = "months"
)

```
