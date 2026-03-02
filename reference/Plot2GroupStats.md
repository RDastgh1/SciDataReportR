# Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE by p; COLOR by Category (vector or data frame); stable point size; palette via paletteer)

Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE
by p; COLOR by Category (vector or data frame); stable point size;
palette via paletteer)

## Usage

``` r
Plot2GroupStats(
  Data,
  Variables,
  VariableCategories = NULL,
  impClust,
  normalClust,
  GroupVar,
  missing_threshold = 0.8,
  max_levels = 10,
  label_q = 0.05,
  x_axis = c("signed_logp", "signed_effect", "effect", "logp"),
  sort_by = c("q", "p", "effect", "signed_logp", "signed_effect", "none"),
  mct_args = list(),
  palette = "pals::alphabet",
  point_size = 3.5
)
```

## Arguments

- Data:

  data.frame

- Variables:

  character vector of variables to analyze

- VariableCategories:

  optional: - data frame with columns Variable, Category; OR - vector of
  categories (named by variable OR unnamed aligned to `Variables`)

- impClust, normalClust:

  labels for the two groups (impClust plotted to the RIGHT for signed
  axes)

- GroupVar:

  column name in `Data` holding the group labels

- missing_threshold:

  drop vars with \> this fraction missing (default 0.80)

- max_levels:

  drop factors with \> this many levels (default 10)

- label_q:

  label threshold using q (default 0.05)

- x_axis:

  one of c("signed_logp","signed_effect","effect","logp")

- sort_by:

  one of c("q","p","effect","signed_logp","signed_effect","none")

- mct_args:

  list of extra args to SciDataReportR::MakeComparisonTable(); e.g.,
  AddEffectSize=TRUE

- palette:

  paletteer palette string for category colors (default
  "pals::alphabet")

- point_size:

  numeric constant for point size (default 3.5)

## Value

list(plot=ggplot, table=gtsummary, pvaltable=data.frame,
data_used=tibble)
