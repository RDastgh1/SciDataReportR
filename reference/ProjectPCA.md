# Project PCA scores onto new data

Use an existing PCA solution (either a PCA object from CreatePCATable or
a loading table) to compute principal component scores on a new dataset.

## Usage

``` r
ProjectPCA(
  Data,
  VarsToReduce = NULL,
  PCAInput,
  InputType = c("PCAObj", "LoadingTable"),
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- Data:

  Data frame on which to project PCA scores.

- VarsToReduce:

  Optional character vector of variable names to use. If NULL, uses all
  variables that appear in both Data and the PCA solution.

- PCAInput:

  Either:

  - the full object returned by CreatePCATable (when InputType =
    "PCAObj"), or

  - a loading table like CreatePCATable()\$LoadingTable (when InputType
    = "LoadingTable").

- InputType:

  One of "PCAObj" or "LoadingTable".

- center:

  Logical; only used when InputType is "LoadingTable". For "PCAObj", the
  centering choice is taken from PCAInput\$ScaleParams\$center and this
  argument is ignored (with a warning if it conflicts).

- scale:

  Logical; only used when InputType is "LoadingTable". For "PCAObj", the
  scaling choice is taken from PCAInput\$ScaleParams\$scale and this
  argument is ignored (with a warning if it conflicts).

## Value

A list with:

- Scores:

  Data frame of projected PCA scores.

- CombinedData:

  Original Data with scores appended as new columns.

- LoadingsUsed:

  Matrix of loadings used for projection.

- PCAObj:

  PCA object used (if InputType is "PCAObj"), otherwise NULL.

- VarsUsed:

  Variables used from Data for projection.

- Center:

  Logical flag indicating whether centering was applied for projection.

- Scale:

  Logical flag indicating whether scaling was applied for projection.
