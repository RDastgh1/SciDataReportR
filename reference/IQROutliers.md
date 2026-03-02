# Detect outliers using the Tukey IQR rule and visualize results

This function identifies potential outliers in a numeric variable using
the Tukey interquartile range (IQR) rule (Tukey, 1977). It returns a
tibble of the detected outlier rows and a ggplot visualization showing
the variable across groups with outlier points highlighted.

## Usage

``` r
IQROutliers(df, Variable, id = NULL, group = NULL)
```

## Arguments

- df:

  A data frame or tibble containing the variable to evaluate.

- Variable:

  A string specifying the name of the numeric variable to test.

- id:

  A string specifying the identifier column to include in the returned
  outlier table. If `NULL`, no ID column is included in the returned
  table. Defaults to `NULL`.

- group:

  A string specifying the grouping or batch column to use on the x-axis
  of the diagnostic plot. If `NULL`, the function will produce a single
  combined boxplot across all rows. Defaults to `NULL`.

## Value

A list with two elements:

- `outlierdf`: a tibble containing only rows flagged as outliers,
  including the ID (when requested), variable, group (when requested),
  and outlier flag.

- `p`: a ggplot2 object showing a boxplot and jittered points colored by
  outlier status.

## Details

The outlier rule is defined as: \$\$value \< Q1 - 1.5 \* IQR \\
\textrm{or} \\ value \> Q3 + 1.5 \* IQR\$\$

Outliers are visually highlighted using jittered points colored red,
while the boxplot remains uncolored to prevent creation of a separate
outlier-only box.

The color aesthetic is mapped only within `geom_jitter()`, ensuring the
boxplot is drawn once per group rather than once per outlier class.
Missing values are ignored when computing quartiles. If `id` or `group`
are provided they must be present in the input data frame.

## References

Tukey, J. W. (1977). *Exploratory Data Analysis*. Addison-Wesley.

## Examples

``` r
if (FALSE) { # \dontrun{
  IQROutliers(df_Revalued_Data, "PS",
              id = "record_id",
              group = "Cohort")

  # Without id or group
  IQROutliers(df_Revalued_Data, "PS", id = NULL, group = NULL)
} # }
```
