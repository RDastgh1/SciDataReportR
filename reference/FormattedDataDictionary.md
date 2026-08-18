# Create a formatted data dictionary table

This function generates a formatted data dictionary table using the
specified data frame. The table includes variable names, labels, types,
and additional formatting based on variable types.

## Usage

``` r
FormattedDataDictionary(
  data,
  digits = 2,
  DataFrame = lifecycle::deprecated(),
  numdecimals = lifecycle::deprecated()
)
```

## Arguments

- data:

  The data frame for which the data dictionary is to be created.

- digits:

  Number of decimals to display for numeric variables (default: 2).

- DataFrame:

  **Deprecated** (since 19.15.0). Use `data` instead.

- numdecimals:

  **Deprecated** (since 19.15.0). Use `digits` instead.

## Value

A formatted data dictionary table (gt object).

## Details

This function requires the `gt` and `codebook` packages. If either is
not installed, the function will return an error.

## Giving the table room

A data dictionary is wide by nature - name, label, type, levels, and
summary statistics all belong on one row - so it reads badly squeezed
into a narrow text column. In a Quarto report the chunk option
`column: screen` handles this. The equivalent on a reference page is to
let the table spill into the page margins and scroll horizontally for
the rest, which the example's `ColumnScreen()` helper does.

The [`max()`](https://rdrr.io/r/base/Extremes.html) clamp in that helper
is what keeps it safe: the table breaks out as far as the margins allow
but never past the edge of the window, so no column is pushed off screen
and the page itself never scrolls sideways.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels first
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# Formatted data dictionary for a subset of variables
dictionary <- FormattedDataDictionary(
  Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "Adiponectin")]
)

# Break the table out to the full width of the window
ColumnScreen <- function(x) {
  htmltools::browsable(htmltools::div(
    style = paste(
      "position: relative;",
      "margin-left: max(-5rem, calc(50% - 50vw));",
      "margin-right: max(-5rem, calc(50% - 50vw));",
      "box-sizing: border-box;",
      "overflow-x: auto;"
    ),
    htmltools::HTML(gt::as_raw_html(x))
  ))
}

ColumnScreen(dictionary)


  

Variable1,2
```
