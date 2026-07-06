# Inspect a scientific data file before import

`InspectFile()` checks common import risks before a file is read into an
analysis workflow. It is especially useful for Excel workbooks, where
multiple sheets, metadata rows, unnamed columns, duplicate column names,
and workbook formatting can change how the data should be interpreted.

## Usage

``` r
InspectFile(
  path,
  sheet = NULL,
  preview_rows = 20,
  check_styles = TRUE,
  check_sheets = TRUE,
  check_header = TRUE,
  quiet = FALSE
)
```

## Arguments

- path:

  Path to the file.

- sheet:

  Sheet name or index for Excel files. If `NULL`, the first sheet is
  inspected.

- preview_rows:

  Number of rows to preview when checking header and column name issues.

- check_styles:

  Logical. For `.xlsx` files, check whether workbook styles or
  formatting exist. This can be slower for large workbooks.

- check_sheets:

  Logical. For Excel files, check whether the workbook has multiple
  sheets.

- check_header:

  Logical. For Excel files, attempt to detect whether the header row is
  not row 1.

- quiet:

  Logical. If `FALSE`, print a compact inspection summary.

## Value

Invisibly returns a list containing file inspection metadata.

## Details

This function does not modify the file or the imported data. It only
reports what it finds.
