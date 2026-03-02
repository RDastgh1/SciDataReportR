# Windsorize

This function performs windsorization on a numeric vector, replacing
values below or above a certain threshold with the corresponding
threshold value.

## Usage

``` r
windsorize(Data, sdlim = 2.5)
```

## Arguments

- Data:

  A numeric vector to be windsorized.

- sdlim:

  The number of standard deviations to use for the windsorization
  threshold. Values below (mean - sdlim \* sd) or above (mean + sdlim \*
  sd) will be replaced with the threshold values.

## Value

The windsorized numeric vector.
