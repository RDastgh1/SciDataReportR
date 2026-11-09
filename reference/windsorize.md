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
  side = "both",
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

- side:

  Character string specifying which tail(s) to winsorize: "both"
  (default), "right" for high values only, or "left" for low values
  only.

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
#>   [1]  0.022219088  1.002990836 -0.510837853 -0.561157279 -2.155919690
#>   [6] -1.299313837 -2.489035135  0.385699811 -0.375022105 -0.411129490
#>  [11]  1.150854208  0.466490107 -1.472774839  0.231793138 -0.353862971
#>  [16]  0.009962680 -1.057261963 -1.322741294  1.637855618 -1.466527144
#>  [21]  1.398309261 -1.357004256  3.263905401 -0.784925328  0.378145899
#>  [26]  1.904726206  0.412795172  0.875169986  0.485634836 -0.237414729
#>  [31] -0.766504204  1.856074490 -1.236186676  0.131550941  0.697692268
#>  [36] -0.852487709  1.047918098  0.509333319 -2.664207884  0.624910380
#>  [41]  1.184505035 -1.388541796  0.723454848  0.152894128  0.373553456
#>  [46] -0.511928978  1.695112620  0.556579110 -0.674325393 -1.715975932
#>  [51] -1.253239856  0.010011130 -1.117293663 -1.433744557 -1.765801527
#>  [56]  0.557843484  0.358130211 -0.663349451 -0.383905045  1.047631454
#>  [61]  1.391354759 -0.134404561 -0.972889822 -0.006030462  0.550473054
#>  [66] -4.052599213 -1.857064544  0.104120717 -0.401377603 -2.447339403
#>  [71] -0.293839109  0.011060185 -1.082572658 -1.892421521 -0.700328850
#>  [76]  0.167939654  0.107125304 -0.355187046 -0.107946668  0.274840230
#>  [81]  0.702565217 -1.320257533  0.333215788  0.350844472  0.658356711
#>  [86]  1.513082170 -0.564915910 -0.341629250  0.111237803 -0.151049754
#>  [91] -0.065553941 -1.324244689 -0.637254700 -1.202807676  0.066051856
#>  [96] -0.691412762 -0.217061106 -1.305980618 -0.170662622 -1.106203142
#> [101]  5.945721243  5.945721243 -6.162991935

# IQR-based winsorization
windsorize(x, method = "iqr", iqrlim = 1.5)
#>   [1]  0.022219088  1.002990836 -0.510837853 -0.561157279 -2.155919690
#>   [6] -1.299313837 -2.489035135  0.385699811 -0.375022105 -0.411129490
#>  [11]  1.150854208  0.466490107 -1.472774839  0.231793138 -0.353862971
#>  [16]  0.009962680 -1.057261963 -1.322741294  1.637855618 -1.466527144
#>  [21]  1.398309261 -1.357004256  2.795032144 -0.784925328  0.378145899
#>  [26]  1.904726206  0.412795172  0.875169986  0.485634836 -0.237414729
#>  [31] -0.766504204  1.856074490 -1.236186676  0.131550941  0.697692268
#>  [36] -0.852487709  1.047918098  0.509333319 -2.664207884  0.624910380
#>  [41]  1.184505035 -1.388541796  0.723454848  0.152894128  0.373553456
#>  [46] -0.511928978  1.695112620  0.556579110 -0.674325393 -1.715975932
#>  [51] -1.253239856  0.010011130 -1.117293663 -1.433744557 -1.765801527
#>  [56]  0.557843484  0.358130211 -0.663349451 -0.383905045  1.047631454
#>  [61]  1.391354759 -0.134404561 -0.972889822 -0.006030462  0.550473054
#>  [66] -3.388886984 -1.857064544  0.104120717 -0.401377603 -2.447339403
#>  [71] -0.293839109  0.011060185 -1.082572658 -1.892421521 -0.700328850
#>  [76]  0.167939654  0.107125304 -0.355187046 -0.107946668  0.274840230
#>  [81]  0.702565217 -1.320257533  0.333215788  0.350844472  0.658356711
#>  [86]  1.513082170 -0.564915910 -0.341629250  0.111237803 -0.151049754
#>  [91] -0.065553941 -1.324244689 -0.637254700 -1.202807676  0.066051856
#>  [96] -0.691412762 -0.217061106 -1.305980618 -0.170662622 -1.106203142
#> [101]  2.795032144  2.795032144 -3.388886984

# Winsorize high-value outliers only
windsorize(x, method = "iqr", iqrlim = 1.5, side = "right")
#>   [1]   0.022219088   1.002990836  -0.510837853  -0.561157279  -2.155919690
#>   [6]  -1.299313837  -2.489035135   0.385699811  -0.375022105  -0.411129490
#>  [11]   1.150854208   0.466490107  -1.472774839   0.231793138  -0.353862971
#>  [16]   0.009962680  -1.057261963  -1.322741294   1.637855618  -1.466527144
#>  [21]   1.398309261  -1.357004256   2.795032144  -0.784925328   0.378145899
#>  [26]   1.904726206   0.412795172   0.875169986   0.485634836  -0.237414729
#>  [31]  -0.766504204   1.856074490  -1.236186676   0.131550941   0.697692268
#>  [36]  -0.852487709   1.047918098   0.509333319  -2.664207884   0.624910380
#>  [41]   1.184505035  -1.388541796   0.723454848   0.152894128   0.373553456
#>  [46]  -0.511928978   1.695112620   0.556579110  -0.674325393  -1.715975932
#>  [51]  -1.253239856   0.010011130  -1.117293663  -1.433744557  -1.765801527
#>  [56]   0.557843484   0.358130211  -0.663349451  -0.383905045   1.047631454
#>  [61]   1.391354759  -0.134404561  -0.972889822  -0.006030462   0.550473054
#>  [66]  -4.052599213  -1.857064544   0.104120717  -0.401377603  -2.447339403
#>  [71]  -0.293839109   0.011060185  -1.082572658  -1.892421521  -0.700328850
#>  [76]   0.167939654   0.107125304  -0.355187046  -0.107946668   0.274840230
#>  [81]   0.702565217  -1.320257533   0.333215788   0.350844472   0.658356711
#>  [86]   1.513082170  -0.564915910  -0.341629250   0.111237803  -0.151049754
#>  [91]  -0.065553941  -1.324244689  -0.637254700  -1.202807676   0.066051856
#>  [96]  -0.691412762  -0.217061106  -1.305980618  -0.170662622  -1.106203142
#> [101]   2.795032144   2.795032144 -12.000000000

# Compare the distribution before and after winsorization. Both panels are
# drawn on the raw data's x range, because free scales would rescale the
# winsorized panel to its own narrower range and hide the very thing being
# demonstrated.
set.seed(42)
x <- c(rnorm(200, mean = 10, sd = 2), 30, 32, -8, -10)
compare_df <- data.frame(
  raw        = x,
  winsorized = windsorize(x, method = "iqr", iqrlim = 1.5)
)

df_Compare <- tidyr::pivot_longer(
  compare_df,
  cols = dplyr::everything(),
  names_to = "Version",
  values_to = "Value"
)
df_Compare$Version <- factor(
  df_Compare$Version,
  levels = c("raw", "winsorized"),
  labels = c("Raw", "Winsorized")
)

ggplot2::ggplot(df_Compare, ggplot2::aes(x = Value)) +
  ggplot2::geom_histogram(bins = 40, na.rm = TRUE) +
  ggplot2::facet_wrap(~ Version, ncol = 1) +
  ggplot2::coord_cartesian(xlim = range(compare_df$raw, na.rm = TRUE)) +
  ggplot2::labs(
    title = "Winsorization pulls outliers to the limits",
    subtitle = "Both panels share the raw data's x range",
    x = "Value", y = "Count"
  ) +
  ggplot2::theme_bw()


# The four extreme values are gone from the tails and have reappeared as
# taller bars at the winsorized limits; nothing was dropped.
range(compare_df$raw)
#> [1] -10  32
range(compare_df$winsorized)
#> [1]  4.985675 15.049195
length(compare_df$raw) == length(compare_df$winsorized)
#> [1] TRUE
```
