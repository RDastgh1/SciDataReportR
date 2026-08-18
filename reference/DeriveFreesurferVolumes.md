# Derive Freesurfer bilateral measures and optional ICV-adjusted ratios

Automatically derives bilateral Freesurfer measures from ASEG and DKT
outputs, using either sums or means, and can create
intracranial-volume-adjusted ratios.

## Usage

``` r
DeriveFreesurferVolumes(
  data,
  icv_var = NULL,
  derive_icv_ratios = TRUE,
  bilateral_method = c("sum", "mean"),
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame containing Freesurfer ASEG and/or DKT variables.

- icv_var:

  Optional single character string naming the intracranial volume column
  to use for ratios. When `NULL` and ratios are requested,
  `EstimatedTotalIntraCranialVol` or `eTIV` is detected automatically.

- derive_icv_ratios:

  Logical. If `TRUE`, derive ICV-adjusted ratios. Default is `TRUE`.

- bilateral_method:

  Character string specifying whether matched left/right measures are
  combined using a `"sum"` or `"mean"`. Default is `"sum"`. Output names
  retain the `_total` suffix for compatibility.

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

- `lh_fusiform_thickness`

- `rh_fusiform_thickness`

When ICV ratios are requested and `icv_var` is `NULL`, supported
intracranial volume columns are:

- `EstimatedTotalIntraCranialVol`

- `eTIV`

If both `EstimatedTotalIntraCranialVol` and `eTIV` are present, the
function checks that they are equivalent before deriving ICV-adjusted
variables. If the two columns differ, the function stops and asks the
user to resolve which intracranial volume variable should be used.

Bilateral sums are computed as:

`left + right`

Bilateral means are computed as:

`(left + right) / 2`

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

thickness_derived <- DeriveFreesurferVolumes(
  df_Thickness,
  derive_icv_ratios = FALSE,
  bilateral_method = "mean"
)

df_freesurfer <- dplyr::bind_cols(
  df_freesurfer,
  DeriveFreesurferVolumes(df_freesurfer)
)

attr(fs_derived, "Freesurfer_derivation_log")
} # }
```
