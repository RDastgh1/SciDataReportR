#' Create a data dictionary for a data frame
#'
#' This function generates a data dictionary for a given data frame, including variable names, labels,
#' types, missing value statistics, summary statistics, and other relevant information.
#'
#' @param DataFrame The data frame for which the data dictionary is to be created.
#' @param numdecimals Number of decimals to display for numeric variables (default: 2).
#'
#' @return A data frame representing the data dictionary.
#'
#'
#' @importFrom sjlabelled get_label get_labels
#' @importFrom codebook skim_codebook
#' @importFrom dplyr select rename left_join
#'
#' @export
Make_DataDictionary <- function(DataFrame, numdecimals = 2) {

  # Extract variable labels and values
  df.var <- sjlabelled::get_label(DataFrame)
  df.val <- sjlabelled::get_labels(DataFrame)

  # Generate summary statistics
  c_skim <- codebook::skim_codebook(DataFrame) %>% as_tibble() %>%
    rename(Variable = skim_variable)

  # Rename the 'factor.n_unique' column to 'n_unique' and handle missing values
  c_skim$n_unique <- c_skim$factor.n_unique
  if("character.n_unique" %in% colnames(c_skim)){
    c_skim$n_unique <- replace(c_skim$factor.n_unique, is.na(c_skim$factor.n_unique), c_skim$character.n_unique)
  }

  # Merge dataframes
  CB <- data.frame(Variable = colnames(DataFrame), Label = df.var)
  CB <- left_join(CB, c_skim %>% select(Variable, skim_type, factor.ordered, n_missing, complete_rate,
                                        n_unique, factor.top_counts, numeric.mean, numeric.median,
                                        numeric.sd, numeric.min, numeric.max, numeric.hist), by = "Variable")

  # Calculate unique values for numeric variables
  numVars <- CB$Variable[CB$skim_type == "numeric"]
  for (var in numVars){
    CB$n_unique[CB$Variable == var] <- length(unique(DataFrame[[var]]))
  }

  # Round numeric variables
  is.num <- sapply(CB, is.numeric)
  CB[is.num] <- lapply(CB[is.num], round, numdecimals)

  # Convert to character and handle missing values
  CB <- lapply(CB, as.character) %>% as.data.frame()
  CB[is.na(CB)] <- " "

  # Rename columns
  CB <- CB %>% rename(Type = skim_type, 'Ordered Factor' = factor.ordered, 'Top Counts' = factor.top_counts,
                      Mean = numeric.mean, Median = numeric.median, SD = numeric.sd,
                      Min = numeric.min, Max = numeric.max, Histogram = numeric.hist)

  return(CB)
}
