#' Negated "in" Operator
#'
#' This custom operator returns TRUE if the element on the left-hand side of the operator is not found
#' in the vector on the right-hand side when using the %in% operator. Otherwise, it returns FALSE.
#'
#' @name not_in_operator
#' @aliases %!in%
#' @aliases %notin%
#' @param x The element to be tested for absence in the vector.
#' @param y The vector in which to search for the element.
#' @return TRUE if the element is not found in the vector, otherwise FALSE.
#' @export
`%!in%` <- function(x, y) {
  # Return the negation of the result of the %in% operator
  !('%in%'(x, y))
}

`%notin%` <- function(x, y) {
  # Return the negation of the result of the %in% operator
  !('%in%'(x, y))
}
