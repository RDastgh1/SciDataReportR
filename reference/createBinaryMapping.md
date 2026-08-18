# Create a Mapping Table for Binary Variables

Identifies binary variables and returns a deterministic mapping for 0/1
coding.

- Factors with *explicit order* (ordered = TRUE) do NOT use heuristics;
  the highest (last) level is Positive.

- Logicals map to Negative = "FALSE", Positive = "TRUE".

- Numeric 0/1 (or any 2-value numeric) maps Positive to the numeric
  maximum.

- Characters / unordered factors use minimal heuristics (no race/PWH/sex
  terms).

## Usage

``` r
createBinaryMapping(
  data,
  CatVars,
  prefer = NULL,
  Data = lifecycle::deprecated()
)
```

## Arguments

- data:

  A dataframe.

- CatVars:

  Character vector of candidate binary variables.

- prefer:

  Optional named character vector of explicit positive levels, e.g.,
  c(STATUS = "PWH", Smoker = "Yes"). This overrides other rules.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

## Value

A data.frame with columns: Variable, Label, PositiveLevel,
NegativeLevel.

## Details

Modelling a two-level variable means scoring one level against the
other, and which level counts as "positive" decides the sign of every
coefficient and odds ratio that follows. This resolves that choice once,
explicitly, and returns it as a table so it can be checked rather than
assumed.

The rules are applied in order: an explicit `prefer` entry always wins;
an ordered factor uses its highest level; otherwise a short list of
conventional affirmative labels ("Yes", "Present", "Case", and similar)
is consulted, then the larger of two numbers, and finally the second
level in sorted order. The heuristics deliberately exclude race, sex,
and serostatus terms, because there is no defensible default "positive"
level for those - name them through `prefer`.

## See also

[`getBinaryVars()`](https://rdastgh1.github.io/SciDataReportR/reference/getBinaryVars.md)
to find the candidates.

## Examples

``` r
# \donttest{
data(SampleData)
data(SampleVariableTypes)

Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
vars_Binary <- getBinaryVars(Labelled)

# One row per variable, naming the level scored as positive
mapping <- createBinaryMapping(Labelled, vars_Binary)

htmltools::browsable(htmltools::HTML(as.character(
  FreezeTableHeader(mapping, full_width = TRUE)
)))


 Variable 
```
