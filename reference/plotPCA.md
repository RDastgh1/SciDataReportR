# Plot 3D PCA Scores

Creates an interactive 3-D scatter plot of PCA scores with Plotly. You
can colour points by a variable **and** customise the hover label.

## Usage

``` r
plotPCA(PCAObj, Var = NULL, t = "NULL", HoverVar = NULL)
```

## Arguments

- PCAObj:

  A list containing:

  - `Scores` – a matrix/data-frame with the PCA scores.

  - `CombinedData` – the original data used to calculate the PCA.

- Var:

  (string) Optional variable in `CombinedData` to colour by. *Default:*
  `NULL` (no colouring).

- t:

  (string) If colouring, set to `"Factor"` for categorical; anything
  else is treated as continuous. *Default:* `"NULL"`.

- HoverVar:

  (string) Optional variable in `CombinedData` whose value will be
  displayed on hover. *Default:* `NULL`, which shows the row number.

## Value

A Plotly htmlwidget.
