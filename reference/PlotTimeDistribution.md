# Plot Time Distribution

This function plots the distribution of time-based data.

## Usage

``` r
PlotTimeDistribution(
  data,
  DateVariable = "Date",
  Data = lifecycle::deprecated()
)
```

## Arguments

- data:

  The data frame containing the time-based data.

- DateVariable:

  The name of the column in the data frame containing the date
  information. Default is "Date".

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A ggplot object displaying the distribution of time-based data.

## Examples

``` r
set.seed(1)
df <- data.frame(
  Date = as.Date("2024-01-01") + sample(0:364, 200, replace = TRUE)
)

PlotTimeDistribution(df, DateVariable = "Date")
```
