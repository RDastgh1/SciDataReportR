# Project new data onto an existing SOM clinical phenotype space

Compatibility alias for
[`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md).
Prefer
[`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md)
in new code.

## Usage

``` r
Project_SOMClust(...)
```

## Arguments

- ...:

  Arguments passed to
  [`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md).

## Value

The same projection object returned by
[`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md).

## See also

[`ProjectSOMCluster()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectSOMCluster.md)
for the canonical function and full examples.

## Examples

``` r
if (FALSE) { # \dontrun{
# NOTE: Not run - see ProjectSOMCluster() and CreateSOMClusterModel() for the
# tracked get_data() bug that blocks the SOM workflow.
Project_SOMClust(object = model, new_df = SampleData)
} # }
```
