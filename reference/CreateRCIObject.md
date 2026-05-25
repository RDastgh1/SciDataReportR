# Create a Reliable Change Index (RCI) object

Learn regression-based Reliable Change Index (RCI) models relative to a
user-defined reference visit and calculate projected RCI values.

## Usage

``` r
CreateRCIObject(
  Data,
  Variables,
  DataFormat = c("wide", "long"),
  ID,
  Method = "regression",
  BaselineSpecifier = NULL,
  FollowupSpecifier = NULL,
  SpecifierPosition = c("suffix", "prefix"),
  VisitColumn = NULL,
  VisitOrder = NULL,
  BaselineVisit = NULL,
  Confidence = 0.95,
  Relabel = TRUE
)
```

## Arguments

- Data:

  A data frame.

- Variables:

  Character vector of canonical variable names.

- DataFormat:

  Either "wide" or "long".

- ID:

  ID column.

- Method:

  Currently only "regression" is supported.

- BaselineSpecifier:

  Baseline visit identifier for wide data.

- FollowupSpecifier:

  Follow-up visit identifier for wide data.

- SpecifierPosition:

  Either "suffix" or "prefix".

- VisitColumn:

  Visit column for long data.

- VisitOrder:

  Optional ordering of visits.

- BaselineVisit:

  Reference visit used for RCI calculations.

- Confidence:

  Confidence interval threshold.

- Relabel:

  Logical; use variable labels when available.

  \#'

  ### Interpretation guide

  egression-based RCI values are interpreted similarly to z-scores.

  |            |                                 |
  |------------|---------------------------------|
  | RCI cutoff | Approximate confidence interval |
  | +/-0.50    | ~38%                            |
  | +/-1.00    | ~68%                            |
  | +/-1.645   | ~90%                            |
  | +/-1.96    | ~95%                            |
  | +/-2.58    | ~99%                            |

  Traditional Jacobson-Truax RCI thresholds typically use +/-1.96,
  corresponding to approximately 95% confidence.

## Value

A SciDataReportR_RCI object.

## Details

Supports both wide and long longitudinal data structures.

Long format is recommended for datasets with more than two visits.
