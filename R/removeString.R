#To do: add warnings if the strings were not in the original

#' Remove Strings from a Vector
#'
#' This function removes strings from a vector that are present in another vector.
#'
#' @param Orig The original vector containing strings.
#' @param Remove The vector of strings to be removed from the original vector.
#' @return A vector containing the strings from the original vector that were not present in the removal vector.
#' @export
removeString <- function(Orig, Remove) {
  leftovers <- Orig[Orig %in% Remove == FALSE]
  return(leftovers)
}
