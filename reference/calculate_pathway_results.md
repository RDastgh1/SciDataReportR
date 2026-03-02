# Calculate Pathway Results for Metabolite Comparisons

Calculate Pathway Results for Metabolite Comparisons

## Usage

``` r
calculate_pathway_results(
  df,
  comparison_var,
  covariates = NULL,
  metabolites,
  comparison_type = "auto",
  use_point_correlation = FALSE
)
```

## Arguments

- df:

  Data frame containing metabolite and comparison data

- comparison_var:

  Character string specifying the comparison variable name

- covariates:

  Character vector of covariate names (optional)

- metabolites:

  Character vector of metabolite names to analyze

- comparison_type:

  Character string: "auto", "binary", or "continuous"

- use_point_correlation:

  Logical, if TRUE uses point correlation for binary comparisons

## Value

Data frame with metabolite results including fold change or correlation
values

## Examples

``` r
if (FALSE) { # \dontrun{
results <- calculate_pathway_results(
  df = my_data,
  comparison_var = "poorSleep",
  covariates = c("Age", "BMI"),
  metabolites = c("Tryptophan", "Kynurenine"),
  comparison_type = "binary"
)
} # }
```
