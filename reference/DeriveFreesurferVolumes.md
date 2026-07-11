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
number of rows as `data`. A derivation log is stored in the attribute
`"Freesurfer_derivation_log"`.

## Details

This function is designed for Freesurfer data frames that contain ASEG,
DKT, and global Freesurfer volume variables. It supports both native
Freesurfer-style ASEG names using hyphens and cleaned names using
underscores.

Supported ASEG-style bilateral pairs include names such as:

- `Left-Hippocampus`

- `Right-Hippocampus`

- `Left_Caudate`

- `Right_Caudate`

Supported DKT-style bilateral cortical pairs include names such as:

- `lh_fusiform_volume`

- `rh_fusiform_volume`

Supported intracranial volume columns are:

- `EstimatedTotalIntraCranialVol`

- `eTIV`

If both `EstimatedTotalIntraCranialVol` and `eTIV` are present, the
function checks that they are equivalent before deriving ICV-adjusted
variables. If the two columns differ, the function stops and asks the
user to resolve which intracranial volume variable should be used.

Bilateral totals are computed as:

`left + right`

ICV-adjusted ratios are computed as:

`volume / intracranial volume`

ICV-adjusted ratios are produced for each individual left/right (or
`lh`/`rh`) measure, for the derived bilateral totals, and for the global
Freesurfer volumes. For example, `Left-Hippocampus` and
`Right-Hippocampus` yield `Left_Hippocampus_icv`,
`Right_Hippocampus_icv`, and `Hippocampus_total_icv`.

The function returns only the newly derived variables, not the original
data. This makes it convenient to append the derived columns with
[`cbind()`](https://rdrr.io/r/base/cbind.html) or
[`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html).

## Examples

``` r
if (FALSE) { # \dontrun{
fs_derived <- DeriveFreesurferVolumes(df_freesurfer)

df_freesurfer <- dplyr::bind_cols(
  df_freesurfer,
  DeriveFreesurferVolumes(df_freesurfer)
)

attr(fs_derived, "Freesurfer_derivation_log")
} # }
```
