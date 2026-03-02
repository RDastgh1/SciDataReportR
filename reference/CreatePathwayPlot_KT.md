# Create Kynurenine-Tryptophan Pathway Plot

Creates a pathway diagram for the kynurenine-tryptophan metabolic
pathway with color-coded fold changes or correlations and significance
indicators.

## Usage

``` r
CreatePathwayPlot_KT(
  results_table,
  title = "",
  value_type = "auto",
  metabolite_mapping = NULL,
  use_fdr = FALSE
)
```

## Arguments

- results_table:

  Data frame with columns: Metabolite, p_value, p_adj, and either "%
  Change" or "correlation"

- title:

  Character string for plot title

- value_type:

  Character string: "auto", "fold_change", or "correlation"

- metabolite_mapping:

  Named character vector mapping results table names to standard names.
  For example: c("N'-Formylkynurenine" = "N-Formylkynurenine",
  "Quinolinic Acid(log10)" = "Quinolinic Acid")

- use_fdr:

  Logical: if TRUE uses FDR-adjusted p-values (p_adj) for significance,
  if FALSE uses raw p-values. Default is FALSE.

## Value

A ggplot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with raw p-values
plot <- CreatePathwayPlot_KT(results, "My Pathway")

# Use FDR-adjusted p-values for significance
plot <- CreatePathwayPlot_KT(results, "My Pathway", use_fdr = TRUE)

# With custom metabolite name mapping
name_map <- c(
  "N'-Formylkynurenine" = "N-Formylkynurenine",
  "Quinolinic Acid(log10)" = "Quinolinic Acid",
  "3-OH-kynurenine" = "3-Hydroxykynurenine"
)
plot <- CreatePathwayPlot_KT(results, "My Pathway", metabolite_mapping = name_map)
} # }
```
