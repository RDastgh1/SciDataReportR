# Plot Bland-Altman Agreement Plot

Generates a Bland-Altman plot to visualize the agreement between two
variables.

## Usage

``` r
PlotBlandAltman(DataFrame, Variable1, Variable2)
```

## Arguments

- DataFrame:

  A data frame containing the variables to compare.

- Variable1:

  The name of the first variable (as a string) to compare.

- Variable2:

  The name of the second variable (as a string) to compare.

## Value

A list containing:

- plot:

  A ggplot2 object of the Bland-Altman plot.

- stats:

  A list of Bland-Altman statistics from the BlandAltmanLeh package.

## Note

This function is adapted from code written by Eran Shorer.
