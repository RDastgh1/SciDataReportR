# Prepare a data frame for SPSS export

`PrepSPSS()` converts column names to SPSS-safe variable names while
preserving the original column names as variable labels. It also returns
and optionally writes a name map so users can track original and
exported names.

## Usage

``` r
PrepSPSS(
  data,
  path = NULL,
  name_map_path = NULL,
  label_map_path = NULL,
  return = c("list", "data", "map"),
  quiet = FALSE,
  show_map = FALSE,
  max_length = 64,
  max_label_length = 120,
  compress = "byte",
  ...
)
```

## Arguments

- data:

  A data frame or tibble.

- path:

  Optional file path ending in `.sav`. If supplied, the prepared data
  are written using
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html).

- name_map_path:

  Optional file path for saving the name map as a CSV.

- label_map_path:

  Optional file path for saving the label map (all truncated labels
  alongside their original text) as a CSV.

- return:

  One of `"list"`, `"data"`, or `"map"`.

- quiet:

  Logical. If `FALSE`, prints a compact summary.

- show_map:

  Logical. If `TRUE`, prints the full original-to-SPSS name map. The
  default is `FALSE` because large scientific datasets can have
  thousands of renamed variables.

- max_length:

  Maximum SPSS variable-name length. Defaults to 64.

- max_label_length:

  Maximum value-label length in bytes. Defaults to 120, the SPSS limit
  enforced by
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html).
  Factor levels and value labels longer than this are truncated and
  recorded in the label map.

- compress:

  Compression type passed to
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html).

- ...:

  Additional arguments passed to
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html)
  if `path` is supplied.

## Value

Depending on `return`, either a list (with elements `data`, `name_map`,
and `label_map`), the prepared data frame, or the name map.

## Details

SPSS caps value labels (factor levels and
[`haven::labelled()`](https://haven.tidyverse.org/reference/labelled.html)
value labels) at 120 bytes and variable labels at 256 bytes;
[`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html)
refuses to write files that exceed these limits. `PrepSPSS()` truncates
over-long labels to fit (marking them with `"..."`), de-duplicates
factor levels that collide after truncation, and records every
truncation in a label map so no information is silently lost.
