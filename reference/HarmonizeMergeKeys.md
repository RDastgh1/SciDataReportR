# Harmonize merge key columns before a safe merge

Internal helper used by safe_merge().

## Usage

``` r
HarmonizeMergeKeys(
  df_before,
  df_add,
  by,
  key_parser = NULL,
  stop_on_failed_numeric = TRUE,
  n_examples = 5
)
```

## Value

A list with elements `left` and `right` (the input data frames with
harmonized key columns) and `report` (a data frame describing the
harmonization applied to each key).
