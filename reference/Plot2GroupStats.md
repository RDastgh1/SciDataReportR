# Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE by p; COLOR by Category (vector or data frame); stable point size; palette via paletteer)

Plot & Summarize Group Stats via MakeComparisonTable (BH q from p; SHAPE
by p; COLOR by Category (vector or data frame); stable point size;
palette via paletteer)

## Usage

``` r
Plot2GroupStats(
  data,
  variables,
  VariableCategories = NULL,
  impClust,
  normalClust,
  group_var,
  missing_threshold = 0.8,
  max_levels = 10,
  label_q = 0.05,
  x_axis = c("signed_logp", "signed_effect", "effect", "logp"),
  sort_by = c("q", "p", "effect", "signed_logp", "signed_effect", "none"),
  mct_args = list(),
  palette = NULL,
  point_size = 3.5,
  Data = lifecycle::deprecated(),
  Variables = lifecycle::deprecated(),
  GroupVar = lifecycle::deprecated()
)
```

## Arguments

- data:

  data.frame

- variables:

  character vector of variables to analyze

- VariableCategories:

  optional:

  - data frame with columns Variable, Category; OR

  - vector of categories (named by variable OR unnamed aligned to
    `Variables`)

- impClust, normalClust:

  labels for the two groups (impClust plotted to the RIGHT for signed
  axes)

- group_var:

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

  Optional paletteer palette string for category colors. When `NULL`
  (the default), the SciDataReportR palette is used. Passing a paletteer
  string such as `"pals::alphabet"` still works as before.

- point_size:

  numeric constant for point size (default 3.5)

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- Variables:

  **Deprecated** (since 19.15.0). Use `variables` instead.

- GroupVar:

  **Deprecated** (since 19.15.0). Use `group_var` instead.

## Value

list(plot=ggplot, table=gtsummary, pvaltable=data.frame,
data_used=tibble)

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable output
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A broad biomarker panel compared between Diagnosis groups
vars <- c(
  "age", "ACE_CD143_Angiotensin_Converti", "ACTH_Adrenocorticotropic_Hormon",
  "AXL", "Adiponectin", "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
  "Alpha_1_Microglobulin", "Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
  "Angiotensinogen", "Apolipoprotein_A_IV", "Apolipoprotein_A1",
  "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
  "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E",
  "Apolipoprotein_H", "B_Lymphocyte_Chemoattractant_BL", "BMP_6",
  "Beta_2_Microglobulin", "Betacellulin", "C_Reactive_Protein", "CD40",
  "CD5L", "Calbindin", "Calcitonin", "CgA", "GRO_alpha", "MMP10", "MMP7",
  "NT_proBNP", "PAI_1", "TRAIL_R3", "VEGF", "Ab_42", "p_tau", "tau"
)

result <- Plot2GroupStats(
  Labelled,
  variables = vars,
  group_var = "Diagnosis",
  impClust = "Impaired",
  normalClust = "Control",
  label_q = 0.0001
)

# Compact y-axis labels; full results stay in result$pvaltable
result$plot + ggplot2::theme(
  axis.text.y = ggplot2::element_text(size = 6),
  plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10)
)

# }
```
