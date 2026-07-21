# Create MCA table and visualization

Compatibility alias for
[`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md).
Prefer
[`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md)
in new code because this workflow returns a reusable MCA object, not
only a static table.

## Usage

``` r
CreateMCATable(...)
```

## Arguments

- ...:

  Arguments passed to
  [`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md).

## Value

The same object returned by
[`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md).

## See also

[`CreateMCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateMCAObject.md)
for the canonical function and full examples.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

mca <- CreateMCATable(Labelled, VarsToReduce = c("Diagnosis", "Genotype"))
#> Warning: no non-missing arguments to min; returning Inf
mca$p_scree
```
