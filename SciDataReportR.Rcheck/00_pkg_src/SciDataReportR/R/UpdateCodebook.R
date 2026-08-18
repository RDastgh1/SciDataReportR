
#' Update an existing codebook based on a given dataframe
#'
#' This function updates a codebook by adding missing variables, removing variables
#' that no longer exist in the dataframe (if specified), and optionally replacing
#' outdated labels.
#'
#' @param data A dataframe for which the codebook needs to be updated.
#' @param codebook A dataframe representing the existing codebook with at least a 'Variable' column.
#' @param RemoveMissing Logical; if TRUE, removes variables from the codebook that are not in the dataframe.
#' @param ReplaceLabels Logical; if TRUE, replaces outdated labels in the codebook with new ones from the dataframe.
#'
#' @return A list containing:
#'   - `UpdatedCodebook`: The updated codebook dataframe.
#'   - `NewVariables`: Variables present in the dataframe but missing from the original codebook.
#'   - `NotExistingVariables`: Variables in the codebook that are not present in the dataframe.
#'   - `MismatchedLabels`: A dataframe of variables with mismatched labels between the codebook and the dataframe.

#' @param Dataframe \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Codebook \strong{Deprecated} (since 19.15.0). Use \code{codebook} instead.
#' @export
UpdateCodebook <- function(data,
    codebook,
    RemoveMissing = TRUE,
    ReplaceLabels = FALSE,
    Dataframe = lifecycle::deprecated(),
    Codebook = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Dataframe)) {
    lifecycle::deprecate_warn("19.15.0", "UpdateCodebook(Dataframe)", "UpdateCodebook(data)")
    data <- Dataframe
  }
  if (!missing(data)) Dataframe <- data
  if (lifecycle::is_present(Codebook)) {
    lifecycle::deprecate_warn("19.15.0", "UpdateCodebook(Codebook)", "UpdateCodebook(codebook)")
    codebook <- Codebook
  }
  if (!missing(codebook)) Codebook <- codebook


  # Ensure the 'MissingCode' column is character type
  Codebook$MissingCode <- as.character(Codebook$MissingCode)

  # Identify variables that are in the Dataframe but missing from the Codebook
  vars_MissingFromCodebook <- setdiff(colnames(Dataframe), Codebook$Variable)

  # If there are missing variables, create a new codebook template for them
  if (length(vars_MissingFromCodebook) > 0) {
    NewCodebook <- CreateVariableTypesTemplate(Dataframe %>% dplyr::select(dplyr::all_of(vars_MissingFromCodebook)))
    Codebook <- bind_rows(Codebook, NewCodebook)
  }

  # Identify variables that exist in the Codebook but are not in the Dataframe
  vars_MissingFromDataFrame <- setdiff(Codebook$Variable, colnames(Dataframe))

  # Remove missing variables from Codebook if RemoveMissing is set to TRUE
  if (RemoveMissing) {
    Codebook <- Codebook %>% dplyr::filter(Variable %in% names(Dataframe))
  }

  # Create a reference codebook for the full Dataframe
  FullCodebook <- CreateVariableTypesTemplate(Dataframe)

  # Identify variables with mismatched labels between the Codebook and FullCodebook
  mismatched_labels <- Codebook %>%
    inner_join(FullCodebook, by = "Variable", suffix = c("_old", "_new")) %>%
    dplyr::filter(Label_old != Variable & Label_old != Label_new) %>%
    dplyr::select(Variable, Label_old, Label_new)

  # Update labels if ReplaceLabels is TRUE
  if (ReplaceLabels) {
    Codebook <- Codebook %>%
      dplyr::left_join(FullCodebook %>% dplyr::select(Variable, Label), by = "Variable", suffix = c("_old", "_new")) %>%
      dplyr::mutate(Label = ifelse(Label_new != Variable, Label_new, Label_old)) %>%
      dplyr::select(Variable, Label)
  }

  # Return a list containing the updated codebook and relevant metadata
  list(
    UpdatedCodebook = Codebook,
    NewVariables = vars_MissingFromCodebook,
    NotExistingVariables = vars_MissingFromDataFrame,
    MismatchedLabels = mismatched_labels
  )
}
