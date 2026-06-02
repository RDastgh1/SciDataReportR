# Create PCA table and visualization

Compatibility alias for
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).
Prefer
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md)
in new code because this workflow returns a reusable PCA object, not
only a static table.

## Usage

``` r
CreatePCATable(...)
```

## Arguments

- ...:

  Arguments passed to
  [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).

## Value

The same object returned by
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).

## Examples

``` r
PCA <- CreatePCATable(
  Data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3,
  imputeMethod = "median"
)
```
