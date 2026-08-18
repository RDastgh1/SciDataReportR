#' Create a Template for Variable Types
#'
#' Generates the variable-types table for a data frame, already in the format
#' [RevalueData()] expects: one row per column, with the variable's name, its
#' label, a guessed type, and empty columns waiting to be filled in. Optionally
#' writes it straight to CSV.
#'
#' @details
#' This is the starting point of the codebook workflow, and its value is that
#' you never build the file by hand. Point it at your data frame and it
#' produces every row, in the right order, with the right column names.
#'
#' @section Editing the template:
#' The generated file is meant to be edited - by you in a spreadsheet, or by
#' the person who actually knows the data. Sending the CSV to a collaborator
#' and asking them to fill it in is usually faster and more accurate than
#' interviewing them about their variables, because the questions are already
#' laid out one per row:
#'
#' * **Label** - what the variable should be called in tables and on plot axes.
#'   `Weight_kg_v2` becomes "Weight (kg)". Every SciDataReportR function that
#'   relabels output reads this.
#' * **Type** - the guess is only a guess. It comes from the stored class, plus
#'   a rule that treats anything with five or fewer distinct values as
#'   categorical. A 0/1 diagnosis stored as a number is the case that always
#'   needs correcting, and an ordinal scale stored as a factor is the other.
#' * **Recode** and **Code** - set `Recode` to `yes` and write the mapping in
#'   `Code` (`0 = Female; 1 = Male`) to turn coded numbers into labelled
#'   factors. Variables read from SPSS arrive with this already filled in from
#'   their value labels.
#' * **MissingCode** - the sentinel values that mean "missing" in this dataset,
#'   such as `999` or `-9`, so they are converted to `NA` instead of being
#'   analyzed as real measurements.
#' * **Notes** - free text, for the reasoning that would otherwise be lost.
#'
#' Editing the file changes nothing on its own. The edits take effect when the
#' table is passed to [RevalueData()], which applies the labels, types,
#' recodings, and missing codes to the data frame and reports anything it could
#' not act on.
#'
#' @section Selecting variables with it:
#' The template is also a convenient place to keep the groupings you analyze
#' by. Extra columns you add are carried along untouched, so a `Category`
#' column ("Demographics", "Lab Measures", "Cognition") or an `Include` flag
#' turns the codebook into the single place variable sets are defined:
#'
#' ```r
#' vars_Labs <- VarTypes$Variable[VarTypes$Category == "Lab Measures"]
#' MakeTable1(data, variables = vars_Labs)
#' ```
#'
#' That way a variable set lives in the codebook rather than being retyped into
#' every function call, and updating the set is a spreadsheet edit rather than
#' a code change. The bundled [SampleVariableTypes] uses `Category`,
#' `Subcategory`, and `Include` this way.
#'
#' @section The worked example:
#' The example starts from a raw extract of the kind that arrives from a study
#' database: cryptic names, a diagnosis coded 0/1, and 999 standing in for
#' missing.
#'
#' Generating the template gives every column a row, with `Type` filled in as a
#' guess - `sex_cd` and `dx_grp` are correctly caught as categorical because
#' they have few distinct values, but nothing yet knows that 999 means missing.
#' The edits are then made in code so the example runs, but in practice this is
#' the CSV coming back from a collaborator and being read with `read.csv()`.
#'
#' The edits take effect only when the table reaches [RevalueData()]: codes
#' become labelled factors, 999 becomes `NA`, and the labels follow the data
#' into every downstream table. The `Category` column then defines variable
#' sets, so later analyses refer to the codebook instead of repeating
#' hard-coded name vectors.
#'
#' @seealso [RevalueData()] to apply an edited template, [UpdateDataDictionary()]
#'   to add rows for new variables without losing existing edits, and
#'   [FormattedDataDictionary()] to render the finished codebook.
#'
#' @param data A data frame containing the variables to be summarized.
#' @param CSVFileName A string specifying the path and name of the CSV file to save the summary.
#'                    If NULL (the default), the CSV file will not be created.
#' @param GuessCategorical A logical variable specifying if the function should guess what variables are categorical based on having <= 5 unique values
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Variable}{The names of the variables in the input data frame.}
#'   \item{Label}{The labels of the variables, if available; otherwise, the variable names.}
#'   \item{Type}{The data types of the variables, converted to more user-friendly descriptions.}
#'   \item{Category}{A placeholder column for categorizing variables (default is NA).}
#'   \item{Recode}{A placeholder column for recoding information (default is NA).}
#'   \item{Code}{A placeholder column for code information (default is NA).}
#'   \item{Notes}{A placeholder column for any additional notes (default is an empty string).}
#'   \item{Exclude}{A placeholder column for exclusion flags (default is NA).}
#' }
#'
#' @examples
#' \donttest{
#' # A raw extract: cryptic names, a 0/1 diagnosis, 999 for missing
#' df_Raw <- data.frame(
#'   subj_id = 1:6,
#'   age_yrs = c(58, 61, 999, 47, 72, 66),
#'   sex_cd = c(0, 1, 1, 0, 1, 0),
#'   dx_grp = c(0, 1, 1, 0, 1, 1),
#'   mmse_tot = c(29, 24, 21, 30, 18, 26),
#'   visit_dt = as.Date(c(
#'     "2024-01-05", "2024-01-11", "2024-02-02",
#'     "2024-02-14", "2024-03-01", "2024-03-19"
#'   ))
#' )
#'
#' # Generate and inspect the template before editing it.
#' VarTypes <- CreateVariableTypesTemplate(df_Raw)
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(VarTypes, full_width = TRUE)
#' )))
#'
#' # Write it out for a collaborator to edit in a spreadsheet
#' path_Template <- file.path(tempdir(), "variable_types.csv")
#' CreateVariableTypesTemplate(df_Raw, path_Template)
#'
#' # The edits they would make, done here in code
#' VarTypes$Label <- c(
#'   "Participant ID", "Age at visit (years)", "Sex",
#'   "Diagnostic group", "MMSE total score", "Visit date"
#' )
#' VarTypes$Recode <- c(NA, NA, "yes", "yes", NA, NA)
#' VarTypes$Code <- c(
#'   NA, NA, "0 = Female; 1 = Male", "0 = Control; 1 = Impaired", NA, NA
#' )
#' VarTypes$MissingCode <- c(NA, "999", NA, NA, NA, NA)
#' VarTypes$Category <- c(
#'   "Identifier", "Demographics", "Demographics",
#'   "Clinical", "Cognition", "Design"
#' )
#'
#' # The edits take effect here
#' revalued <- RevalueData(df_Raw, VarTypes)
#' df_Labelled <- revalued$RevaluedData
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_Labelled, full_width = TRUE)
#' )))
#'
#' # The Category column now defines variable sets
#' vars_Demographics <- VarTypes$Variable[VarTypes$Category == "Demographics"]
#' vars_Demographics
#' }
#'
#' @importFrom sjlabelled get_label
#' @importFrom utils write.csv
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
CreateVariableTypesTemplate <- function(data,
    CSVFileName = NULL,
    GuessCategorical = TRUE,
    DataFrame = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "CreateVariableTypesTemplate(DataFrame)", "CreateVariableTypesTemplate(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data


  # Get the classes of the variables in the DataFrame
  Types <- sapply(DataFrame, class)
  Types <- factor(Types, levels = c("numeric", "integer", "factor",
                                    "character", "Date", "logical"), labels = c("Double",
                                                                                "Double", "Categorical", "String", "Date", "Categorical"))

  # Decide if categorical Values need to be calculated based on whether there were any categorical variables or not.
  if(GuessCategorical == T || (is.null(GuessCategorical) & sum(Types == "Categorical", na.rm = T) == 0)){

    # Get the number of unique values for each. If it's <=5, decide that it is categorical
    unique_counts <- sapply(DataFrame, function(x) length(unique(x)))

    Types[unique_counts <= 5] <- "Categorical"
  }

  # Get variable labels (if available)
  DataLabels <- sjlabelled::get_label(DataFrame, def.value = colnames(DataFrame))

  # Create the template dataframe
  VariableTypes <- data.frame(Variable = colnames(DataFrame),
                              Label = DataLabels, Type = Types, Category = NA, Recode = NA,
                              Code = NA, Notes = "", Exclude = NA, MissingCode = "")

  # Handle labelled factors (e.g., variables loaded from SPSS using `haven`)
  labelled_factors <- sapply(DataFrame, sjlabelled::is_labelled)

  for (i in which(labelled_factors)) {
    var_name <- colnames(DataFrame)[i]

    # Set Recode to 1 for labelled factors
    VariableTypes$Recode[i] <- 1

    # Get the labels for the factor and create a recoding string
    labels <- levels(as.factor(DataFrame[[var_name]]))
    codes <- sjlabelled::get_labels(DataFrame[[var_name]])

    recode_str <- paste0(paste(labels, codes, sep = "=",
                               collapse = "; "))
    # Set the Code column to show recoding
    VariableTypes$Code[i] <- recode_str
  }

  # Optionally save the output to a CSV file
  if (!is.null(CSVFileName)) {
    write.csv(VariableTypes, CSVFileName, row.names = FALSE)
  }

  return(VariableTypes)
}
