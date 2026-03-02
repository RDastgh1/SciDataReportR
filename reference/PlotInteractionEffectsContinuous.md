# Plot Single Interaction Effect

Creates a scatter plot with regression lines showing the interaction
between a predictor and outcome variable, moderated by either a
continuous or categorical variable.

## Usage

``` r
PlotInteractionEffectsContinuous(
  Data,
  interVar = NULL,
  outcomeVar = NULL,
  predictorVar = NULL,
  covars = NULL,
  n_lines = 3,
  alpha = 0.6,
  point_size = 2
)
```

## Arguments

- Data:

  A data frame containing the variables to be analyzed

- interVar:

  Character string specifying the interaction variable (moderator)

- outcomeVar:

  Character string specifying the outcome variable

- predictorVar:

  Character string specifying the predictor variable

- covars:

  Character vector of covariate names to include in the model

- n_lines:

  For continuous moderators, number of lines to plot (default: 3 for
  low/med/high)

- alpha:

  Transparency level for points (default: 0.6)

- point_size:

  Size of points (default: 2)

## Value

A ggplot object showing the interaction effect

## Details

For categorical moderators, separate regression lines are plotted for
each category. For continuous moderators, regression lines are plotted
at the mean and ± 1 SD.

The subtitle shows the p-value for the interaction term. The caption
lists any covariates included in the model.

Variable labels are used if available in the data frame.
