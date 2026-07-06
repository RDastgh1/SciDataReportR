# Create Statistics Table

Generate a table of statistics including means, standard deviations,
counts, and p-values.

## Usage

``` r
CreateStatisticsTable(data, TargetVar, Data = lifecycle::deprecated())
```

## Arguments

- data:

  The data frame containing the variables of interest.

- TargetVar:

  The target variable for which statistics will be calculated.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A formatted HTML table displaying statistics.
