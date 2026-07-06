# Winsorize a numeric vector using SD or IQR thresholds

This function performs winsorization on a numeric vector by capping
extreme values at calculated lower and upper thresholds. Thresholds can
be based on either standard deviation (assuming approximate normality)
or interquartile range (robust to skewed distributions).

## Usage

``` r
windsorize(
  data,
  sdlim = 2.5,
  iqrlim = 1.5,
  method = "sd",
  Data = lifecycle::deprecated()
)
```

## Arguments

- data:

  A numeric vector to be winsorized.

- sdlim:

  Numeric. Number of standard deviations for the "sd" method.

- iqrlim:

  Numeric. Multiplier for the IQR when method = "iqr" (default 1.5).

- method:

  Character string specifying the method: "sd" (default) or "iqr".

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A numeric vector with values winsorized to the specified thresholds.

## Details

Winsorization replaces extreme values with the nearest threshold value
rather than discarding them, limiting the influence of outliers while
retaining the sample size. The term originates with the statistician
Charles P. Winsor; see Dixon (1960) and Tukey (1962) for early
treatments of the estimator.

## References

Tukey, J. W. (1962). The future of data analysis. *The Annals of
Mathematical Statistics*, 33(1), 1-67.
[doi:10.1214/aoms/1177704711](https://doi.org/10.1214/aoms/1177704711)

Dixon, W. J. (1960). Simplified estimation from censored normal samples.
*The Annals of Mathematical Statistics*, 31(2), 385-391.
[doi:10.1214/aoms/1177705900](https://doi.org/10.1214/aoms/1177705900)

## Examples

``` r
x <- c(rnorm(100), 10, 15, -12)

# SD-based winsorization
windsorize(x, method = "sd", sdlim = 2.5)
#>   [1]  0.609996541 -1.095720589 -0.086880289 -0.467032116 -0.198470682
#>   [6]  0.748168528 -0.379906361 -1.542796961  0.544144501  0.443726297
#>  [11]  0.835608454 -0.289657352 -0.550534916  0.194222034  1.998376399
#>  [16] -0.137923175 -0.853814398 -2.591133289 -0.692715894 -0.315631545
#>  [21] -0.459876820  0.089198132 -0.080633540  0.340769115 -0.008257906
#>  [26]  1.571428502 -0.409596987 -0.301478813  0.395054457 -0.595236302
#>  [31]  1.893824910 -0.056627324 -0.680579015 -1.298496001 -0.189737210
#>  [36]  0.498595657  0.276169603 -0.175747774  0.639139078 -0.914278603
#>  [41] -0.759582570 -0.164117502  1.276754886 -0.131387834 -1.382509590
#>  [46] -0.494431142  0.971682365 -1.721934135 -0.032492757 -2.364609048
#>  [51] -0.749479370 -0.096079198  0.024409192 -0.791149391 -0.098692346
#>  [56] -0.170913302 -1.758703501 -1.645349576 -1.324402538 -0.993901371
#>  [61] -0.388928127  0.738527947 -0.076932442 -0.671862312 -0.908052285
#>  [66]  0.217315420  1.317197056 -0.283807239 -0.960399000  0.741116538
#>  [71]  0.031334319 -0.063206104 -0.524516678  0.380639054 -1.671511371
#>  [76] -1.501237627  0.727867477 -0.989464103  0.254849570  1.595667140
#>  [81] -1.927738017  0.276319937 -0.151877985 -0.010967372 -0.821258403
#>  [86]  0.691245925 -0.497975407  0.985494560  0.042463701 -0.380926864
#>  [91] -0.486110616 -0.864881658  1.576573561 -0.799783737  1.702113048
#>  [96]  0.423103056  1.471057078  0.024319341  0.297111046  0.798915206
#> [101]  5.800377607  5.800377607 -5.846774438

# IQR-based winsorization
windsorize(x, method = "iqr", iqrlim = 1.5)
#>   [1]  0.609996541 -1.095720589 -0.086880289 -0.467032116 -0.198470682
#>   [6]  0.748168528 -0.379906361 -1.542796961  0.544144501  0.443726297
#>  [11]  0.835608454 -0.289657352 -0.550534916  0.194222034  1.998376399
#>  [16] -0.137923175 -0.853814398 -2.416362216 -0.692715894 -0.315631545
#>  [21] -0.459876820  0.089198132 -0.080633540  0.340769115 -0.008257906
#>  [26]  1.571428502 -0.409596987 -0.301478813  0.395054457 -0.595236302
#>  [31]  1.893824910 -0.056627324 -0.680579015 -1.298496001 -0.189737210
#>  [36]  0.498595657  0.276169603 -0.175747774  0.639139078 -0.914278603
#>  [41] -0.759582570 -0.164117502  1.276754886 -0.131387834 -1.382509590
#>  [46] -0.494431142  0.971682365 -1.721934135 -0.032492757 -2.364609048
#>  [51] -0.749479370 -0.096079198  0.024409192 -0.791149391 -0.098692346
#>  [56] -0.170913302 -1.758703501 -1.645349576 -1.324402538 -0.993901371
#>  [61] -0.388928127  0.738527947 -0.076932442 -0.671862312 -0.908052285
#>  [66]  0.217315420  1.317197056 -0.283807239 -0.960399000  0.741116538
#>  [71]  0.031334319 -0.063206104 -0.524516678  0.380639054 -1.671511371
#>  [76] -1.501237627  0.727867477 -0.989464103  0.254849570  1.595667140
#>  [81] -1.927738017  0.276319937 -0.151877985 -0.010967372 -0.821258403
#>  [86]  0.691245925 -0.497975407  0.985494560  0.042463701 -0.380926864
#>  [91] -0.486110616 -0.864881658  1.576573561 -0.799783737  1.702113048
#>  [96]  0.423103056  1.471057078  0.024319341  0.297111046  0.798915206
#> [101]  2.104343341  2.104343341 -2.416362216

# Compare the distribution before and after winsorization
set.seed(42)
x <- c(rnorm(200, mean = 10, sd = 2), 30, 32, -8, -10)
compare_df <- data.frame(
  raw        = x,
  winsorized = windsorize(x, method = "iqr", iqrlim = 1.5)
)

PlotContinuousDistributions(compare_df, variables = c("raw", "winsorized"))

```
