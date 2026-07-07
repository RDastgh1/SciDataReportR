# Derive Freesurfer bilateral totals and ICV-adjusted measures

Automatically derives bilateral Freesurfer volume measures from ASEG and
DKT outputs, then creates intracranial-volume-adjusted ratios.

## Usage

``` r
DeriveFreesurferVolumes(data, verbose = TRUE)
```

## Arguments

- data:

  A data frame containing Freesurfer ASEG and/or DKT variables.

- verbose:

  Logical. If `TRUE`, prints a short summary of derived variables.
  Default is `TRUE`.

## Value

A data frame containing only newly derived variables, with the same
number of rows as `data`.

## Details

This function is designed for Freesurfer data frames that already
contain cleaned variable names such as:

- `Left_Hippocampus`

- `Right_Hippocampus`

- `lh_fusiform_volume`

- `rh_fusiform_volume`

- `EstimatedTotalIntraCranialVol`

- `eTIV`

The function detects:

- ASEG-style bilateral pairs using `Left_` and `Right_`

- DKT-style bilateral cortical pairs using `lh_` and `rh_`

- selected global Freesurfer volume variables

Bilateral totals are computed as:

`left + right`

ICV-adjusted ratios are computed as:

`volume / intracranial volume`

If both `EstimatedTotalIntraCranialVol` and `eTIV` are present, the
function checks that they are equivalent before deriving ICV-adjusted
variables. If the two columns differ, the function stops and asks the
user to resolve which intracranial volume variable should be used.

The function returns only the newly derived variables, not the original
data. This makes it convenient to append the derived columns with
[`cbind()`](https://rdrr.io/r/base/cbind.html) or
[`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html).

## Examples

``` r
if (FALSE) { # \dontrun{
fs_derived <- DeriveFreesurferVolumes(df_Freesurfer)

df_Freesurfer <- cbind(
  df_Freesurfer,
  DeriveFreesurferVolumes(df_Freesurfer)
)
} # }
```
