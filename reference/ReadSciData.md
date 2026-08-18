# Read a scientific data file with optional inspection

`ReadSciData()` imports common scientific data file formats while
preserving original column names and labels as much as possible. It
optionally calls
[`InspectFile()`](https://rdastgh1.github.io/SciDataReportR/reference/InspectFile.md)
before import to flag common issues such as multiple sheets, metadata
rows, unnamed columns, duplicate column names, and formatting.

## Usage

``` r
ReadSciData(
  path,
  sheet = NULL,
  header_row = NULL,
  col_names = TRUE,
  range = NULL,
  inspect = TRUE,
  print_inspection = interactive(),
  strict = FALSE,
  guess_max = 10000,
  delim = NULL,
  repair_names = TRUE,
  inspect_styles = FALSE,
  fast_delimited = TRUE,
  ...
)
```

## Arguments

- path:

  Path to the file.

- sheet:

  Sheet name or index for Excel files.

- header_row:

  Row containing column names for Excel files. If `NULL`,
  `ReadSciData()` may use the probable header row detected by
  [`InspectFile()`](https://rdastgh1.github.io/SciDataReportR/reference/InspectFile.md).

- col_names:

  Passed to file readers where applicable. For Excel files, use `TRUE`,
  `FALSE`, or a character vector.

- range:

  Optional Excel range passed to
  [`readxl::read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html).

- inspect:

  Logical. If `TRUE`, inspect the file before importing.

- print_inspection:

  Logical. If `TRUE`, print compact inspection results when issues are
  detected. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- strict:

  Logical. If `TRUE`, stop when inspection detects potential issues.

- guess_max:

  Maximum rows used for type guessing where supported.

- delim:

  Delimiter for `.txt` files. Defaults to tab.

- repair_names:

  Logical. If `TRUE`, repair blank and duplicate column names after
  import.

- inspect_styles:

  Logical. If `TRUE`, inspect Excel workbook formatting. This can be
  slower and noisier for large workbooks.

- fast_delimited:

  Logical. If `TRUE`, use
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  for delimited text files when available and when no extra reader
  arguments are supplied through `...`. This is usually much faster than
  `readr` for large `.csv`, `.tsv`, and `.txt` files.

- ...:

  Additional arguments passed to the underlying reader.

## Value

Imported data object. Inspection metadata is attached as the
`scidata_inspection` attribute when `inspect = TRUE`.

## Details

After import, `ReadSciData()` repairs only names that make the data
frame difficult or impossible to use in tidyverse workflows. Blank or
`NA` names are renamed to `...unnamed_POSITION`, and duplicate names are
made unique only when needed. Existing valid names are preserved.

## Formats, and what the examples show

The file extension decides the reader, so one call covers CSV, Excel,
SPSS, Stata, SAS, Parquet, Feather, JSON, and R's own formats.

**Excel.** Real workbooks have more than one sheet, and the data is
rarely on the first one, so `sheet` is usually needed. Exported sheets
also tend to carry a title and an export stamp above the real column
names; read naively those become the header and every column arrives as
text. With `inspect = TRUE` (the default) the real header row is
detected and used, and what inspection found is attached as the
`scidata_inspection` attribute so the import stays auditable.
`header_row` overrides the guess when detection gets it wrong.

**SPSS.** `.sav` files carry variable labels and value labels, and those
are the reason to import them directly rather than converting to CSV
first - a CSV export drops both. Because the labels survive the round
trip, the codebook template arrives pre-filled: `Label` comes from the
variable labels, and `Recode`/`Code` from the value labels, so coded
SPSS variables need no manual transcription.

## See also

[`InspectFile()`](https://rdastgh1.github.io/SciDataReportR/reference/InspectFile.md)
to inspect a file without importing it,
[`CreateVariableTypesTemplate()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateVariableTypesTemplate.md)
to build a codebook from the imported frame, and
[`RevalueData()`](https://rdastgh1.github.io/SciDataReportR/reference/RevalueData.md)
to apply it.

## Examples

``` r
# \donttest{
# Written to a temporary directory so the example is self-contained
dir_Example <- tempdir()

df_Study <- data.frame(
  subj_id = 1:8,
  age = c(58, 61, 47, 72, 66, 54, 69, 43),
  sex = c(0, 1, 1, 0, 1, 0, 1, 0),
  mmse = c(29, 24, 30, 18, 26, 27, 21, 28)
)

# CSV
path_Csv <- file.path(dir_Example, "study_data.csv")
utils::write.csv(df_Study, path_Csv, row.names = FALSE)

df_FromCsv <- ReadSciData(path_Csv)
df_FromCsv
#> # A tibble: 8 × 4
#>   subj_id   age   sex  mmse
#>     <int> <int> <int> <int>
#> 1       1    58     0    29
#> 2       2    61     1    24
#> 3       3    47     1    30
#> 4       4    72     0    18
#> 5       5    66     1    26
#> 6       6    54     0    27
#> 7       7    69     1    21
#> 8       8    43     0    28

# Excel: a workbook with more than one sheet
path_Excel <- file.path(dir_Example, "study_data.xlsx")
openxlsx::write.xlsx(
  list(
    ReadMe = data.frame(Note = "Exported by the lab core"),
    Measurements = df_Study
  ),
  path_Excel
)

df_FromExcel <- ReadSciData(path_Excel, sheet = "Measurements")
df_FromExcel
#> # A tibble: 8 × 4
#>   subj_id   age   sex  mmse
#>     <dbl> <dbl> <dbl> <dbl>
#> 1       1    58     0    29
#> 2       2    61     1    24
#> 3       3    47     1    30
#> 4       4    72     0    18
#> 5       5    66     1    26
#> 6       6    54     0    27
#> 7       7    69     1    21
#> 8       8    43     0    28

# A sheet with a title and export stamp above the real column names
path_Messy <- file.path(dir_Example, "messy_export.xlsx")
openxlsx::write.xlsx(
  as.data.frame(rbind(
    c("Study XYZ biomarker export", NA, NA, NA),
    c("Exported 2026-01-05 by lab core", NA, NA, NA),
    c("subj_id", "age", "sex", "mmse"),
    c("1", "58", "0", "29"),
    c("2", "61", "1", "24"),
    c("3", "47", "1", "30")
  )),
  path_Messy,
  colNames = FALSE
)

# Without inspection
ReadSciData(path_Messy, inspect = FALSE)
#> ReadSciData(): 3 unnamed columns renamed to ...unnamed_2, ...unnamed_3, ...unnamed_4.
#> # A tibble: 5 × 4
#>   `Study XYZ biomarker export`    ...unnamed_2 ...unnamed_3 ...unnamed_4
#>   <chr>                           <chr>        <chr>        <chr>       
#> 1 Exported 2026-01-05 by lab core NA           NA           NA          
#> 2 subj_id                         age          sex          mmse        
#> 3 1                               58           0            29          
#> 4 2                               61           1            24          
#> 5 3                               47           1            30          

# With inspection (the default)
df_Fixed <- ReadSciData(path_Messy)
df_Fixed
#> # A tibble: 3 × 4
#>   subj_id age   sex   mmse 
#>   <chr>   <chr> <chr> <chr>
#> 1 1       58    0     29   
#> 2 2       61    1     24   
#> 3 3       47    1     30   

# What inspection found, kept with the data
inspection <- attr(df_Fixed, "scidata_inspection")
inspection$probable_header_row
#> [1] 3
inspection$issues
#> [1] "Probable header row appears to be row 3, not row 1."

# Overriding the detected header row
ReadSciData(path_Messy, header_row = 3)
#> # A tibble: 3 × 4
#>   subj_id age   sex   mmse 
#>   <chr>   <chr> <chr> <chr>
#> 1 1       58    0     29   
#> 2 2       61    1     24   
#> 3 3       47    1     30   

# SPSS, with variable labels and value labels
df_Labelled <- df_Study
df_Labelled$sex <- haven::labelled(
  df_Labelled$sex,
  labels = c(Female = 0, Male = 1),
  label = "Sex assigned at birth"
)
attr(df_Labelled$age, "label") <- "Age at visit (years)"
attr(df_Labelled$mmse, "label") <- "MMSE total score"

path_Spss <- file.path(dir_Example, "study_data.sav")
haven::write_sav(df_Labelled, path_Spss)

df_FromSpss <- ReadSciData(path_Spss)
df_FromSpss
#> # A tibble: 8 × 4
#>   subj_id   age sex         mmse
#>     <dbl> <dbl> <dbl+lbl>  <dbl>
#> 1       1    58 0 [Female]    29
#> 2       2    61 1 [Male]      24
#> 3       3    47 1 [Male]      30
#> 4       4    72 0 [Female]    18
#> 5       5    66 1 [Male]      26
#> 6       6    54 0 [Female]    27
#> 7       7    69 1 [Male]      21
#> 8       8    43 0 [Female]    28

# The labels survive the round trip
sjlabelled::get_label(df_FromSpss)
#>                 subj_id                     age                     sex 
#>                      ""  "Age at visit (years)" "Sex assigned at birth" 
#>                    mmse 
#>      "MMSE total score" 
sjlabelled::get_labels(df_FromSpss$sex)
#> [1] "Female" "Male"  

# So the codebook template arrives pre-filled
CreateVariableTypesTemplate(df_FromSpss)
#>         Variable                 Label        Type Category Recode
#> subj_id  subj_id               subj_id      Double       NA     NA
#> age          age  Age at visit (years)      Double       NA     NA
#> sex          sex Sex assigned at birth Categorical       NA      1
#> mmse        mmse      MMSE total score      Double       NA     NA
#>                     Code Notes Exclude MissingCode
#> subj_id             <NA>            NA            
#> age                 <NA>            NA            
#> sex     0=Female; 1=Male            NA            
#> mmse                <NA>            NA            
# }
```
