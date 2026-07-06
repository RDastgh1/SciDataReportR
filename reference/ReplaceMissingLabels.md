# Replace Missing Labels in Dataframe Columns

This function iterates through the columns of a dataframe and assigns
the column name as the label to any column that does not have a label.

## Usage

``` r
ReplaceMissingLabels(data, df = lifecycle::deprecated())
```

## Arguments

- data:

  A dataframe.

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

The input dataframe with missing labels replaced.

## Examples

``` r
# A data frame where only some columns carry a variable label
df <- data.frame(
  age    = c(52, 61, 77),
  bmi    = c(24.1, 29.7, 22.0),
  smoker = c(0, 1, 0)
)
labelled::var_label(df$age) <- "Age (years)"

# Before: bmi and smoker have no label
sapply(df, function(x) sjlabelled::get_label(x))
#> $age
#> [1] "Age (years)"
#> 
#> $bmi
#> NULL
#> 
#> $smoker
#> NULL
#> 

# After: unlabelled columns are labelled with their column name
filled <- ReplaceMissingLabels(df)
sapply(filled, function(x) sjlabelled::get_label(x))
#>           age           bmi        smoker 
#> "Age (years)"         "bmi"      "smoker" 
```
