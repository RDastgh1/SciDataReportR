# Combine Two Codebooks with Conflict Detection

`CombineCodebooks` compares two codebook data frames (an old version and
a new version), identifies added or removed variables and columns,
detects cell-by-cell differences, and produces a combined codebook. For
records without any differences, it returns a single "Combined" row; for
records with differences, it returns both the "Old" and "New" rows,
flagged as conflicts.

## Usage

``` r
CombineCodebooks(OldCodebook, NewCodebook, key = "Variable")
```

## Arguments

- OldCodebook:

  A data.frame or tibble representing the old codebook. Each row must
  correspond to a single variable entry.

- NewCodebook:

  A data.frame or tibble representing the new codebook. Must have the
  same structure (column names) as `OldCodebook`, though extra or
  missing columns will be handled.

- key:

  A string giving the name of the key column to identify variables (e.g.
  "Variable"). Defaults to "Variable".

## Value

A list with elements:

- `added_variables`: character vector of keys present only in
  `NewCodebook`.

- `removed_variables`: character vector of keys present only in
  `OldCodebook`.

- `columns_added`: character vector of column names present only in
  `NewCodebook`.

- `columns_removed`: character vector of column names present only in
  `OldCodebook`.

- `value_differences`: tibble of cell-level differences (`RowID`,
  `Field`, `OldValue`, `NewValue`).

- `combined_df`: tibble containing the merged codebook with versions and
  conflict flag.

## Details

This function:

- Builds a `combined_df` that includes all columns from both inputs:

  - For non-conflicting records, a single row with
    `Version = "Combined"`.

  - For conflicting records, two rows: one with `Version = "Old"` and
    one with `Version = "New"`.

- Flags each row with `CHECKFORCONFLICTS` (0 for combined, 1 for old/new
  conflict rows).
