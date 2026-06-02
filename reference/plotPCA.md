# Plot PCA scores

Create an interactive 2D or 3D Plotly scatter plot from a PCA object
created by
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).

## Usage

``` r
plotPCA(
  PCAObj,
  Var = NULL,
  t = "NULL",
  HoverVar = NULL,
  HoverVars = NULL,
  Components = NULL,
  Mode = c("auto", "3D", "2D"),
  ColorType = c("auto", "factor", "continuous"),
  Relabel = TRUE,
  Title = NULL
)
```

## Arguments

- PCAObj:

  A PCA object returned by
  [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).
  Must contain `Scores` and `CombinedData`.

- Var:

  Optional character string naming a variable in `PCAObj$CombinedData`
  used to color points. Default is `NULL`.

- t:

  Deprecated compatibility argument for color type. Use `ColorType`
  instead. If set to `"Factor"`, `Var` is treated as categorical. Any
  other value treats `Var` as continuous. Default is `"NULL"`.

- HoverVar:

  Optional character string naming one variable in `PCAObj$CombinedData`
  to display in hover text. Retained for backward compatibility. Prefer
  `HoverVars` for new code. Default is `NULL`.

- HoverVars:

  Optional character vector naming one or more variables in
  `PCAObj$CombinedData` to display in hover text. If `NULL`, row number
  is shown. Default is `NULL`.

- Components:

  Optional character or numeric vector specifying which score columns to
  plot. Supply two components for 2D or three components for 3D. If
  `NULL`, the first three score columns are used.

- Mode:

  Character. Either `"auto"`, `"3D"`, or `"2D"`. If `"auto"`, the plot
  dimension is inferred from `Components`. If `Components = NULL`,
  `"auto"` defaults to 3D. Default is `"auto"`.

- ColorType:

  Character. Either `"auto"`, `"factor"`, or `"continuous"`. `"auto"`
  treats character, factor, logical, and labelled variables as
  categorical and numeric variables as continuous. Default is `"auto"`.

- Relabel:

  Logical. If `TRUE`, labels attached to hover variables are used in
  hover text when available. If `FALSE`, raw variable names are used.
  Default is `TRUE`.

- Title:

  Optional plot title. Default is `NULL`, which produces no title.

## Value

A Plotly htmlwidget.

## Details

By default, the function plots the first three score columns in
`PCAObj$Scores` using a 3D plot. These are not assumed to be named
`RC1`, `RC2`, and `RC3`; the function uses the current column order in
`PCAObj$Scores`. Users may manually specify which components to plot
using `Components`.

If `Mode = "auto"`, the function chooses 2D when two components are
supplied and 3D when three components are supplied. If
`Components = NULL`, it defaults to the first three score columns and
creates a 3D plot.

The function can optionally color points by a variable in
`PCAObj$CombinedData` and customize hover text using one or more
variables. When `Relabel = TRUE`, labels attached to hover variables are
used in the hover display when available.

## Examples

``` r
if (FALSE) { # \dontrun{
plotPCA(PCAObj)

plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  Var = "SITE",
  ColorType = "factor"
)

plotPCA(
  PCAObj,
  Components = c("RC2", "RC1", "RC4"),
  Mode = "3D"
)

plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  HoverVars = c("SubjectID", "SITE", "Visit")
)

plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  HoverVars = c("SubjectID", "SITE", "Visit"),
  Relabel = FALSE
)
} # }
```
