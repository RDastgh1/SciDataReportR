# Insert Values into a String Array

This function inserts one or more strings into an array of strings at
specified locations.

## Usage

``` r
InsertValues(vec, values, after_value, location = "after")
```

## Arguments

- vec:

  A character vector in which values will be inserted.

- values:

  A character vector or a single string to insert into the original
  vector.

- after_value:

  A character string indicating the element in `vec` before or after
  which the values will be inserted. Should have a length of one

- location:

  A character string specifying whether to insert `values` "before" or
  "after" the `after_value`. Defaults to "after".

## Value

A character vector with the values inserted.
