# SOM + latent profile clustering pipeline (with AHP and distance baselines)

Compatibility alias for
[`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md).
Prefer
[`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md)
in new code because this function fits a reusable SOM clustering model.

## Usage

``` r
Pipeline_SOMClust(...)
```

## Arguments

- ...:

  Arguments passed to
  [`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md).

## Value

The same `Pipeline_SOMClust` object returned by
[`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md).

## See also

[`CreateSOMClusterModel()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateSOMClusterModel.md)
for the canonical function and full examples.

## Examples

``` r
if (FALSE) { # \dontrun{
# NOTE: Not run - see CreateSOMClusterModel() for the tracked get_data() bug.
data(SampleData)

Pipeline_SOMClust(
  data = SampleData,
  variables = c("age", "AXL", "Adiponectin", "Alpha_1_Antitrypsin"),
  method = "exploratory",
  k_range = 2:4,
  models = 1
)
} # }
```
