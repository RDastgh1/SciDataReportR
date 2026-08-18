# Neutral simulated clustering and phenotyping benchmark

A fixed 480-participant benchmark for evaluating train-once/project-many
clustering workflows. Twelve numeric variables support at least four
retained principal components and four balanced truth groups. Three
four-level categorical variables provide strongly separated,
class-dependent response patterns for categorical-method demonstrations,
while `DensityX` and `DensityY` provide two unequal-density groups plus
noise. `TruthCluster` and `TruthDensityGroup` are retained solely for
teaching and test evaluation; they must not be used as clustering
features.

## Usage

``` r
data(SimulatedPhenotypeData)
```

## Format

`SimulatedPhenotypeData` has 480 rows: 320 training participants and 160
projection participants. It contains `Var1` through `Var12`, `CatVar1`
through `CatVar3`, density coordinates, truth fields, and the fixed
cohort split. `SimulatedPhenotypeVariableTypes` is its complete labelled
codebook.

## Why the benchmark works

The twelve numeric variables are organized in blocks, and each cluster
is high on a different block. That is the structure a clustering method
is supposed to recover, and the reason a method that finds nothing here
is genuinely failing rather than facing an impossible problem.

`DensityX` and `DensityY` are a separate and deliberately harder
problem: two groups of unequal density plus a noise group, for methods
that are meant to detect outliers rather than partition everything.

`TruthCluster` and `TruthDensityGroup` exist only so results can be
scored. They must never be used as clustering features.

## See also

[SampleData](https://rdastgh1.github.io/SciDataReportR/reference/SampleData.md)
for the labelling and reporting examples, and
[`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md)
or
[`CreateClusterModel_PCA_KMeans()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_PCA_KMeans.md)
for the pipelines this benchmark was built for.

## Examples

``` r
data(SimulatedPhenotypeData)

# 480 participants, with the truth clusters balanced across both cohorts
dim(SimulatedPhenotypeData)
#> [1] 480  21
table(SimulatedPhenotypeData$Cohort, SimulatedPhenotypeData$TruthCluster)
#>             
#>              Cluster 1 Cluster 2 Cluster 3 Cluster 4
#>   Projection        40        40        40        40
#>   Training          80        80        80        80

# \donttest{
data(SimulatedPhenotypeVariableTypes)

ShowTable <- function(x, caption = NULL, height = NULL) {
  htmltools::browsable(htmltools::HTML(as.character(
    kableExtra::kable_styling(
      knitr::kable(x, format = "html", caption = caption, row.names = FALSE),
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    )
  )))
}

# A slice: identifiers, clustering variables, and the truth labels
ShowTable(
  utils::head(
    SimulatedPhenotypeData[, c(
      "ParticipantID", "Cohort", "TruthCluster",
      paste0("Var", 1:6), "CatVar1", "DensityX", "DensityY"
    )],
    8
  ),
  "First eight participants"
)

First eight participants
 
 ParticipantID 
```
