# SampleData for practicing SciDataReportR functions

A sample dataset adapted from the AppliedPredictive Modeling R package
with additional simulated variables, noise and missingness added. The
original dataset Craig-Schapiro et al. (2011) describes a clinical study
of 333 patients where laboratory measurements are used to predict which
subjects are most likely to develop cognitive impairment, such as
Alzheimer's disease.

## Usage

``` r
SampleData
```

## Format

### `SampleData`

A data frame with 333 rows and 131 columns:

- Diagnosis:

  Control or impaired

- age:

  Age in years

- sex:

  Coded 0 = Female, 1 = Male, left uncoded until
  [`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
  is applied

- Genotype:

  APOE genotype, one of E2E2, E2E3, E2E4, E3E3, E3E4, E4E4

- ...:

  127 further columns, almost all continuous protein and biomarker
  measurements

## Source

<http://appliedpredictivemodeling.com/data>

## Details

This is deliberately a *raw* extract rather than a tidy one, because
that is what the package is for. It arrives with no variable labels,
`sex` stored as bare 0/1, and missingness scattered through 13 of the
biomarker columns (up to 30% in some), so the examples throughout the
package can show what each function does to real, untidied data.

It is paired with
[SampleVariableTypes](https://rdastgh1.github.io/SciDataReportR/reference/SampleVariableTypes.md),
the codebook that describes it. Passing the two to
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
is the first step of nearly every example in this package.

## See also

[SampleVariableTypes](https://rdastgh1.github.io/SciDataReportR/reference/SampleVariableTypes.md)
for the companion codebook, and
[SimulatedPhenotypeData](https://rdastgh1.github.io/SciDataReportR/reference/SimulatedPhenotypeData.md)
for the clustering examples.

## Examples

``` r
data(SampleData)

# 333 participants, 131 columns, two diagnosis groups.
dim(SampleData)
#> [1] 333 131
table(SampleData$Diagnosis)
#> 
#>  Control Impaired 
#>      242       91 

# The first columns are demographics; the rest are biomarkers.
names(SampleData)[1:10]
#>  [1] "Diagnosis"                       "age"                            
#>  [3] "sex"                             "Genotype"                       
#>  [5] "ACE_CD143_Angiotensin_Converti"  "ACTH_Adrenocorticotropic_Hormon"
#>  [7] "AXL"                             "Adiponectin"                    
#>  [9] "Alpha_1_Antichymotrypsin"        "Alpha_1_Antitrypsin"            

# As shipped, it is unlabelled and `sex` is still a bare numeric code.
str(SampleData[, c("Diagnosis", "age", "sex", "Genotype", "AXL")])
#> 'data.frame':    333 obs. of  5 variables:
#>  $ Diagnosis: chr  "Control" "Control" "Control" "Control" ...
#>  $ age      : num  52 61 77 97 73 87 82 56 69 89 ...
#>  $ sex      : int  0 0 1 0 0 1 1 1 0 0 ...
#>  $ Genotype : chr  "E3E3" "E3E4" "E3E4" "E3E4" ...
#>  $ AXL      : num  1.098 0.683 -0.145 0.683 0.191 ...
sjlabelled::get_label(SampleData$age)
#> NULL

# \donttest{
data(SampleVariableTypes)

# The codebook turns it into the labelled frame the other examples use.
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(
    utils::head(
      Labelled[, c("Diagnosis", "age", "sex", "Genotype", "AXL", "tau", "p_tau")],
      8
    ),
    full_width = TRUE
  )
)))


 Diagnosis 
```
