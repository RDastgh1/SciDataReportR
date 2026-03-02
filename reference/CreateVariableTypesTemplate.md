# Create a Template for Variable Types

This function generates a data frame that summarizes the types and
labels of the variables in a given data frame. It optionally saves this
summary to a CSV file.

## Usage

``` r
CreateVariableTypesTemplate(
  DataFrame,
  CSVFileName = NULL,
  GuessCategorical = TRUE
)
```

## Arguments

- DataFrame:

  A data frame containing the variables to be summarized.

- CSVFileName:

  A string specifying the path and name of the CSV file to save the
  summary. If NULL (the default), the CSV file will not be created.

- GuessCategorical:

  A logical variable specifying if the function should guess what
  variables are categorical based on having \<= 5 unique values

## Value

A data frame with the following columns:

- Variable:

  The names of the variables in the input data frame.

- Label:

  The labels of the variables, if available; otherwise, the variable
  names.

- Type:

  The data types of the variables, converted to more user-friendly
  descriptions.

- Category:

  A placeholder column for categorizing variables (default is NA).

- Recode:

  A placeholder column for recoding information (default is NA).

- Code:

  A placeholder column for code information (default is NA).

- Notes:

  A placeholder column for any additional notes (default is an empty
  string).

- Exclude:

  A placeholder column for exclusion flags (default is NA).

## Examples

``` r
df <- data.frame(
  num = c(1.1, 2.2),
  int = c(1L, 2L),
  fact = factor(c("A", "B")),
  char = c("a", "b"),
  date = as.Date(c("2021-01-01", "2021-01-02"))
)
CreateVariableTypesTemplate(df)
#>      Variable Label        Type Category Recode Code Notes Exclude MissingCode
#> num       num   num Categorical       NA     NA   NA            NA            
#> int       int   int Categorical       NA     NA   NA            NA            
#> fact     fact  fact Categorical       NA     NA   NA            NA            
#> char     char  char Categorical       NA     NA   NA            NA            
#> date     date  date Categorical       NA     NA   NA            NA            
CreateVariableTypesTemplate(df, "variable_types.csv")
#>      Variable Label        Type Category Recode Code Notes Exclude MissingCode
#> num       num   num Categorical       NA     NA   NA            NA            
#> int       int   int Categorical       NA     NA   NA            NA            
#> fact     fact  fact Categorical       NA     NA   NA            NA            
#> char     char  char Categorical       NA     NA   NA            NA            
#> date     date  date Categorical       NA     NA   NA            NA            
```
