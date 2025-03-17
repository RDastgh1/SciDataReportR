
#' Update an existing codebook based on a given dataframe
#'
#' This function updates a codebook by adding missing variables, removing variables
#' that no longer exist in the dataframe (if specified), and optionally replacing
#' outdated labels.
#'
#' @param Dataframe A dataframe for which the codebook needs to be updated.
#' @param Codebook A dataframe representing the existing codebook with at least a 'Variable' column.
#' @param RemoveMissing Logical; if TRUE, removes variables from the codebook that are not in the dataframe.
#' @param ReplaceLabels Logical; if TRUE, replaces outdated labels in the codebook with new ones from the dataframe.
#'
#' @return A list containing:
#'   - `UpdatedCodebook`: The updated codebook dataframe.
#'   - `NewVariables`: Variables present in the dataframe but missing from the original codebook.
#'   - `NotExistingVariables`: Variables in the codebook that are not present in the dataframe.
#'   - `MismatchedLabels`: A dataframe of variables with mismatched labels between the codebook and the dataframe.

#' @export
UpdateCodebook <- function(Dataframe, Codebook, RemoveMissing = TRUE, ReplaceLabels = FALSE) {

  # Ensure the 'MissingCode' column is character type
  Codebook$MissingCode <- as.character(Codebook$MissingCode)

  # Identify variables that are in the Dataframe but missing from the Codebook
  vars_MissingFromCodebook <- setdiff(colnames(Dataframe), Codebook$Variable)

  # If there are missing variables, create a new codebook template for them
  if (length(vars_MissingFromCodebook) > 0) {
    NewCodebook <- CreateVariableTypesTemplate(Dataframe %>% select(all_of(vars_MissingFromCodebook)))
    Codebook <- bind_rows(Codebook, NewCodebook)
  }

  # Identify variables that exist in the Codebook but are not in the Dataframe
  vars_MissingFromDataFrame <- setdiff(Codebook$Variable, colnames(Dataframe))

  # Remove missing variables from Codebook if RemoveMissing is set to TRUE
  if (RemoveMissing) {
    Codebook <- Codebook %>% filter(Variable %in% names(Dataframe))
  }

  # Create a reference codebook for the full Dataframe
  FullCodebook <- CreateVariableTypesTemplate(Dataframe)

  # Identify variables with mismatched labels between the Codebook and FullCodebook
  mismatched_labels <- Codebook %>%
    inner_join(FullCodebook, by = "Variable", suffix = c("_old", "_new")) %>%
    filter(Label_old != Variable & Label_old != Label_new) %>%
    select(Variable, Label_old, Label_new)

  # Update labels if ReplaceLabels is TRUE
  if (ReplaceLabels) {
    Codebook <- Codebook %>%
      left_join(FullCodebook %>% select(Variable, Label), by = "Variable", suffix = c("_old", "_new")) %>%
      mutate(Label = ifelse(Label_new != Variable, Label_new, Label_old)) %>%
      select(Variable, Label)
  }

  # Return a list containing the updated codebook and relevant metadata
  list(
    UpdatedCodebook = Codebook,
    NewVariables = vars_MissingFromCodebook,
    NotExistingVariables = vars_MissingFromDataFrame,
    MismatchedLabels = mismatched_labels
  )
}
