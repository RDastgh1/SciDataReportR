# Plot significant interaction pairs

Select interaction rows from a matrix result and reuse
[`PlotInteractionEffectsContinuous()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotInteractionEffectsContinuous.md)
for each single-pair display.

## Usage

``` r
plotSigInteractions(
  data,
  InteractionObject,
  PVar = "PValue",
  Pthresh = 0.05,
  adjust = TRUE,
  standardize = TRUE,
  ncol = NULL,
  max_plots = 25
)
```

## Arguments

- data:

  Data used to fit the source matrix.

- InteractionObject:

  Output of
  [`PlotInteractionEffectsMatrix()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotInteractionEffectsMatrix.md)
  with Results and Metadata. Plot layer data is never used for
  selection.

- PVar:

  Selection column: "PValue" or "PValueFDR".

- Pthresh:

  Select p-values strictly below this cutoff.

- adjust:

  Display model-based covariate-adjusted partial residuals. The
  predictor is residualized on nuisance covariates and moderator terms;
  the outcome partial residual is the full-model residual plus the
  fitted predictor contribution on that residualized predictor.
  Independent covariate-only residualization of both axes would not
  preserve interaction slopes in general. FALSE displays raw values and
  an unadjusted model.

- standardize:

  Scale axes using full fitted-sample SDs, pooled over groups.

- ncol:

  Optional positive integer requesting a patchwork layout.

- max_plots:

  Maximum number of plots, ordered by selected p-value.

## Value

Named list of ggplots, with selected interaction rows in attribute
Results. When ncol is supplied, attribute Layout holds a patchwork.
Empty selections return an empty list with a message.

## Details

Annotations report displayed-model standardized simple slopes, intervals
and sample sizes. Interaction estimate, p and q retain the original
source inference even when adjust is FALSE. Metadata and source formulas
are retained as attributes. Continuous moderator lines are at the
model-sample mean minus/plus one SD; point colours remain continuous.
