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
  return = c("list", "data", "map"),
  quiet = FALSE,
  show_map = FALSE,
  max_length = 64,
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

- compress:

  Compression type passed to
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html).

- ...:

  Additional arguments passed to
  [`haven::write_sav()`](https://haven.tidyverse.org/reference/read_spss.html)
  if `path` is supplied.

## Value

Depending on `return`, either a list, data frame, or name map.
