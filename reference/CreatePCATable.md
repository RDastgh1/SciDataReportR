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

## See also

[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md)
for the canonical function and full examples.

## Examples

``` r
PCA <- CreatePCATable(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

# Display the scree plot from the returned object
PCA$p_scree

```
