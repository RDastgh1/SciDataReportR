# Read a scientific data file with optional inspection

`ReadSciData()` imports common scientific data file formats while
preserving original column names and labels as much as possible. It
optionally calls
[`InspectFile()`](https://rdastgh1.github.io/SciDataReportR/reference/InspectFile.md)
before import to flag common issues such as multiple sheets, metadata
rows, unnamed columns, duplicate column names, and formatting.

## Usage

``` r
ReadSciData(
  path,
  sheet = NULL,
  header_row = NULL,
  col_names = TRUE,
  range = NULL,
  inspect = TRUE,
  print_inspection = interactive(),
  strict = FALSE,
  guess_max = 10000,
  delim = NULL,
  repair_names = TRUE,
  inspect_styles = FALSE,
  fast_delimited = TRUE,
  ...
)
```

## Arguments

- path:

  Path to the file.

- sheet:

  Sheet name or index for Excel files.

- header_row:

  Row containing column names for Excel files. If `NULL`,
  `ReadSciData()` may use the probable header row detected by
  [`InspectFile()`](https://rdastgh1.github.io/SciDataReportR/reference/InspectFile.md).

- col_names:

  Passed to file readers where applicable. For Excel files, use `TRUE`,
  `FALSE`, or a character vector.

- range:

  Optional Excel range passed to
  [`readxl::read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html).

- inspect:

  Logical. If `TRUE`, inspect the file before importing.

- print_inspection:

  Logical. If `TRUE`, print compact inspection results when issues are
  detected. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- strict:

  Logical. If `TRUE`, stop when inspection detects potential issues.

- guess_max:

  Maximum rows used for type guessing where supported.

- delim:

  Delimiter for `.txt` files. Defaults to tab.

- repair_names:

  Logical. If `TRUE`, repair blank and duplicate column names after
  import.

- inspect_styles:

  Logical. If `TRUE`, inspect Excel workbook formatting. This can be
  slower and noisier for large workbooks.

- fast_delimited:

  Logical. If `TRUE`, use
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  for delimited text files when available and when no extra reader
  arguments are supplied through `...`. This is usually much faster than
  `readr` for large `.csv`, `.tsv`, and `.txt` files.

- ...:

  Additional arguments passed to the underlying reader.

## Value

Imported data object. Inspection metadata is attached as the
`scidata_inspection` attribute when `inspect = TRUE`.

## Details

After import, `ReadSciData()` repairs only names that make the data
frame difficult or impossible to use in tidyverse workflows. Blank or
`NA` names are renamed to `...unnamed_POSITION`, and duplicate names are
made unique only when needed. Existing valid names are preserved.
